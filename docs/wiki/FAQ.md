# FAQ

**A sample failed. Where do I look?**
`summary.tsv` gives the step (`failed (assemble)`) and the reason in `Note`; `logs/<step>/<sample>.log` has
the commands and the programs' messages (fastp's own report is in `logs/2_filter/<sample>.fastp.html`). A failed
sample does not stop the others, and is retried when RASHER runs again in the same output folder.

**"no reads matched the reference"**
The sample has no read aligning to the reference: wrong reference, a failed library, or an empty sample. With a
distant reference, try `-b bbduk`.

**My R1 and R2 files became two samples.**
They are paired only when they sit in the same folder and their names differ only by the mate marker before the
extension: `X_R1`/`X_R2`, `X_1`/`X_2` or `X.R1`/`X.R2` (an `_001` tail allowed), with the same marker on both
(`X_1` does not pair with `X_R2`). Names such as `X-R1.fastq.gz` or `X_read1.fastq.gz` are not recognized, in a
folder or in a sample sheet (which pairs a sample's files by the same rules): rename them, or make symbolic links
named `X_R1.fastq.gz` and `X_R2.fastq.gz`. See [Usage](Usage#inputs).

**"desynchronise at read ..." or "do not hold the same number of reads"**
The R1 and R2 files of a pair must list the mates in the same order. Files that were filtered or sorted
separately, or truncated by an interrupted download, break that. Take the files as the sequencer or the SRA wrote
them (`fasterq-dump` keeps the order).

**A warning says a file was skipped next to the pairs of a sample.**
A sample (a subfolder, or a sample's rows of a sheet) made of pairs cannot also take single reads: the file
without a mate is left out. Give it as a sample of its own to use it.

**How do I merge the lanes of a sample?**
In a flat folder, each lane of a bcl2fastq run (`X_S1_L001_R1_001`, `X_S1_L002_R1_001`) is a sample of its own.
Put a sample's files in a subfolder named after it, or list them in a sample sheet: each lane's pair stays a pair
within the sample.

**"No SNP site is shared by all the genomes: no tree"**
With `--snp-method ska`, core SNPs must have their context in every genome. One fragmented or incomplete
assembly is enough to remove every site. Look at `Contigs` and `Length_vs_reference` in `summary.tsv`, remove
the poor samples, or use `--ska-min-freq 0.9` to keep SNPs missing from a few genomes.

**"No tree: no genome differs from the others at any position compared"**
The VCF of the SNPs (all of them, whatever `--ska-min-freq`) is empty: the genomes are identical at every
position compared, so there is nothing to build a tree on, and `--ska-min-freq` would not change that. The
distances are all 0.

**Which `--snp-method`?**
The default, `auto`, compares templated assemblies with minimap2 and de novo assemblies (SPAdes, SKESA) with SKA2.
SKA2 needs the same 15 bases on either side of a SNP in every genome: it loses the SNPs that lie next to an indel
or another variant, and those near the ends of contigs and runs of `N`. `--snp-method minimap2` reads each assembly
position by position along the reference and needs no such context, but it counts every difference an assembly
carries, including the errors of a de novo assembly, which SKA2 skips: hence SKA2 for de novo assemblies. These
defaults come from BACoN and held up on the simulated short-read [validation](Validation). See
[Methods](Methods#4-comparison).

**What does `--ska-min-freq` below 1 do?**
SKA2 then keeps the SNP sites that some genomes lack (they get a gap in the alignment), as long as at least that
fraction of the genomes has them: a pan-genome SNP alignment instead of core SNPs. It helps when core SNPs are few
or none: fragmented assemblies, gaps, indels or other SNPs within 15 bases of a site, genomes of another
structure. The cost is missing data: each pair of genomes is compared on the sites both have, and a genome
missing many sites looks closer to the others than it is. `snps.vcf` is the same whatever the value. Details in
[Methods](Methods#4-comparison).

**Which assembler?**
The templated assembly (the default) to compare samples of one species against a close reference. For a de novo
check or a distant reference, SKESA (the recommended de novo assembler: the same SNP distances as SPAdes on the
real-data benchmark, in far fewer contigs); SPAdes for the assembly graph, below roughly 30x depth, or where
heavy NUMT/NUPT contamination is expected. See [Methods](Methods#3-assembly).

**The de novo assembly is much longer than the reference.**
Reads from other genomes were baited with the target: mitochondrial and nuclear sequences similar to the
plastome (for example plastid DNA inserted in the mitochondrion), or contamination. The extra contigs are
usually at low coverage, and `Off_reference_pct` gives the share of the assembly that aligns nowhere on the
reference (a note above 10%). The templated assembly is not affected; for SPAdes assemblies, check the assembly
graph in `3_assembled/assembly_graphs/`.

**I am baiting a mitogenome.**
Plant mitogenomes contain plastid-derived stretches, and nuclear genomes copies of mitochondrial DNA: reads of
both are baited with the mitogenome's, and count towards `--target-depth` (see [Methods](Methods#2-filtering)).
With whole-genome reads of a large nuclear genome, check the real depth: the mitogenome may be covered only a few
times. Look at the `N` bases of templated assemblies, and at the off-reference notes of de novo ones.

**My chloroplast assembly is in several contigs, not one circle.**
Short reads, and their inserts, are much shorter than the inverted repeat (about 25 kb): no assembler can tell
from them how the single-copy regions connect through the repeat, so the large single copy, the repeat (once)
and the small single copy come out separately (SKESA breaks at other repeats too). SPAdes' assembly graph shows
the connections. This does not affect SNP calling with SKA2; the templated assembly follows the reference's
structure.

**Why are there `N` in the templated assembly?**
Where fewer than three reads cover the reference, or no base reaches a call quality of 20 (the reads disagree,
often inside an insertion), samtools writes `N`. `N_bases` in `summary.tsv` counts them. A de novo assembly
resolves insertions.

**Which reads can I use?**
Illumina short reads (or other short reads of similar accuracy), paired or single-end, fastq or fasta. fastp
trims the adapters and filters by quality and length; fasta reads have no qualities and are filtered by length
only. Long reads belong to [BACoN](https://github.com/duceppemo/BACoN).

**Do I need to trim the reads first?**
No: fastp trims adapters (and polyG tails of NextSeq and NovaSeq data) with its defaults. Reads already trimmed
work too; reads shorter than `--min-read-length` (50 bp) are dropped.

**Can I compare samples sequenced on different instruments or with different read lengths?**
Yes. Systematic errors differ between runs and can add to the distances between samples from different runs.
The templated assembly guards against them by needing three reads and a call quality of 20 at every base (`N`
otherwise); compare `N_bases` between the runs.

**How do I compare with a public genome?**
Give its fasta file with `--add-genomes`: finished genomes (published plastomes, earlier assemblies) join
the comparison as they are, named after their file, without baiting or assembly.

```bash
rasher -r ref.fasta -i reads/ -o out/ --add-genomes published/*.fasta
```

**How do I resume, or change one parameter?**
Run the same command again; see [Usage](Usage#resuming-and-changing-parameters).

**Can RASHER resume a BACoN output folder?**
No: every step runs again (the reads, filters and assemblers differ). Use a new output folder.
