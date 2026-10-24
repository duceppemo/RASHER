# Tutorial: plastomes of 24 Arabidopsis accessions

This tutorial runs RASHER on public data: whole-genome Illumina reads of 24 *Arabidopsis thaliana*
accessions (ENA project [PRJEB62038](https://www.ebi.ac.uk/ena/browser/view/PRJEB62038); HiSeq 2x100 bp
pairs, 9–16 Gbases per run). The reference is the Col-0 plastome
([NC_000932.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_000932.1), 154,478 bp) — and Col-0 itself is among
the 24 accessions, so one sample should come out identical to the reference. Let's check, and see how the
other 23 plastomes differ. This is the dataset of the 0.1.0 [benchmark](Validation#real-data-benchmark-010);
the same project holds matched Nanopore runs of the same accessions, which is why it was chosen there.

## 1. Download the reads and the reference

About 160 GB of fastq (this is whole-genome data; RASHER does the skimming). The project holds over a
hundred accessions; this tutorial uses the 24 that also have a matched Nanopore run (the
[benchmark](Validation#real-data-benchmark-010) set, so the two platforms can be compared on the same
samples). The script downloads those 24 paired runs from ENA, named after their accession, and checks
their MD5:

```bash
mkdir arabidopsis && cd arabidopsis
ACCESSIONS="Are-1 Are-10 Are-2 Are-6 Bay-0 Bla-1 Blh-1 Bur-0 Can-0 Cant-1 Col-0 Ct-1 Cvi-0 Ge-0 Ita-0 Jea Kz-9 Nok-1 Oy-0 Ri-0 Shahdara Shigu-2 Tsu-0 Yo-0"
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB62038&result=read_run&fields=sample_alias,fastq_ftp,fastq_md5&format=tsv&limit=0" \
    | awk -F'\t' -v acc="$ACCESSIONS" 'BEGIN{split(acc,a," "); for(i in a) keep[a[i]]} \
        NR>1 && ($2 in keep) && $3 ~ /_1.fastq.gz/' > runs.tsv
mkdir -p illumina
while IFS=$'\t' read -r run alias urls md5s; do
    name=$(echo "$alias" | tr ' /' '__')_${run}
    url1=${urls%%;*}; url2=${urls##*;}
    curl -s -C - --retry 5 -o "illumina/${name}_R1.fastq.gz" "https://$url1"
    curl -s -C - --retry 5 -o "illumina/${name}_R2.fastq.gz" "https://$url2"
    echo "${md5s%%;*}  illumina/${name}_R1.fastq.gz"
    echo "${md5s##*;}  illumina/${name}_R2.fastq.gz"
done < runs.tsv > md5.txt
md5sum -c md5.txt

# The reference: the Col-0 plastome, as fasta and as the annotated GenBank record
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000932.1&rettype=fasta" \
    > NC_000932.1.fasta
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000932.1&rettype=gbwithparts&retmode=text" \
    > NC_000932.1.gb
```

The `awk` keeps the 24 accessions' paired Illumina runs (the project also holds their Nanopore runs, which
have no `_1`/`_2` files). If an MD5 check fails, run the loop again: `curl -C -` resumes interrupted
downloads. The `_R1`/`_R2` files pair into one sample per accession by name ([Usage](Usage#inputs)).

## 2. Run RASHER

```bash
rasher -r NC_000932.1.fasta -i illumina/ -o rasher_plastome -t 24 -p 3 --annotation NC_000932.1.gb
```

The defaults: baiting with minimap2 `-x sr` (a pair is kept when either mate matches), fastp filtering and
an even subsampling to 100x, templated assembly with `samtools consensus`, SNPs read along the reference
with minimap2, and a FastTree tree. With 24 threads it takes about 4.5 hours, nearly all of it baiting the
9–16 Gbases of each run; everything downstream works on the baited ~20%. `--annotation` gives the report
the genes and regions of the plastome and the effect of each SNP; it changes nothing else and can be added
to a finished run (nothing is redone).

Open `rasher_plastome/report.html` in a browser for everything below on one page. The report of this run is
online:
[view it](https://duceppemo.github.io/RASHER/reports/tutorial_arabidopsis_report.html)
or [download it](https://raw.githubusercontent.com/duceppemo/RASHER/main/docs/reports/tutorial_arabidopsis_report.html)
(local paths removed).

## 3. Reads and assemblies: `summary.tsv`

```
Sample              Status  Raw_reads  Raw_bases    Baited_reads  Baited_bases  Baited_pct  ...
Are-10_ERR11436050  ok      73723854   11058578100  14630030      2194504500    19.844
Are-1_ERR11436047   ok      109224206  16383630900  20190292      3028543800    18.485
Are-2_ERR11436048   ok      86115886   12917382900  17521862      2628279300    20.347
...
```

- 9–21% of the bases of these whole-genome runs are plastid reads — leaf DNA, as in BACoN's potato
  tutorial. Nuclear copies of plastome segments (NUPTs) are baited too; they are a small fraction and the
  consensus turns their conflicting bases into `N` rather than SNPs ([Methods](Methods#3-assembly)).
- Every sample is capped at 100x (`Est_depth`), so the assemblies cost the same whatever the run size.
- Every consensus is full length (154,421–154,518 bp, one contig, both inverted repeats), with 11–146 `N`
  bases (median 19).

## 4. SNPs and tree: `4_compared/minimap2/`

The SNP alignment has 503 sites that every plastome covers and where they differ. What they say:

- **Col-0 is 0 SNPs from the reference** — the expected result, since NC_000932.1 *is* the Col-0 plastome;
  a nice end-to-end check of the whole pipeline on real data.
- The other accessions are 39–97 SNPs from the reference, and up to 109 from each other.
- Jea and Ri-0 carry identical plastomes (distance 0); every other pair is distinct.

`snps.vcf` holds every position where the genomes differ (more than the 503 sites every genome covers:
[Outputs](Outputs#comparison)), with the gene, codon and effect of each SNP from the annotation, and the
report's genome map places them on the LSC, inverted repeats and SSC.

## 5. A de novo check: `-a skesa`

```bash
rasher -r NC_000932.1.fasta -i illumina/ -o rasher_plastome -a skesa -t 24 -p 3
```

The same output folder reuses the baited and filtered reads (nothing upstream is redone), and SKESA — the
recommended de novo assembler ([Methods](Methods#3-assembly)) — assembles each sample in a few minutes: 1–3
contigs totalling about 128 kb, which is the plastome with its inverted repeat assembled once, as short
reads must ([FAQ](FAQ)). De novo assemblies are compared with SKA2 by default (`4_compared/ska/`). The
0.1.0 [benchmark](Validation#real-data-benchmark-010) ran this dataset through every assembler, Snippy,
GetOrganelle and NOVOPlasty: the templated distances match Snippy's (r = 0.976) and SKESA matches SPAdes
pair for pair, so the de novo run is a structure check more than a SNP improvement here.

The same accessions work against the *mitochondrial* genome
([NC_037304.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_037304.1), 367,808 bp) by changing `-r` (and `-o`):
the templated assembly is again a single full-length sequence, while the mitogenome's repeats fragment any
short-read de novo assembly — the benchmark record shows that contrast in detail.

## What this shows

- Whole-genome Illumina runs need no prior organelle enrichment: RASHER extracts the ~20% of plastid reads
  and caps the depth, so 24 runs of 9–16 Gbases compare in an afternoon, unattended.
- The templated consensus reproduces the reference exactly where the sample *is* the reference (Col-0), and
  its SNP distances match an independent caller (Snippy) on the same reads.
- A de novo assembly (SKESA) confirms the picture independently and shows the structure; expect the
  inverted repeat once, not twice.
