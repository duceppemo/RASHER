# RASHER

RASHER takes Illumina short reads, paired or single-end, from many samples and a reference sequence, typically an
organelle genome (chloroplast, mitochondrion) or another region of interest, and for each sample:

1. **baits** the reads that match the reference (minimap2, or BBDuk k-mers); a pair is kept whole when either
   mate matches;
2. **filters** them with fastp (adapters, base qualities, length), and caps the depth by an even, deterministic
   subsampling;
3. **assembles** them: a templated consensus of the reads aligned to the reference (samtools, the default), or a
   de novo assembly (SPAdes, or SKESA);

then **compares** the samples: SNPs (each assembly read along the reference with minimap2, the default for
templated assemblies; SKA2 split k-mers, the default for de novo ones; or a Parsnp core-genome alignment), a
matrix of pairwise SNP distances, a VCF of the SNPs relative to the reference, and a tree (FastTree or IQ-TREE).
Everything is summarized in an HTML report, and in MultiQC sections (next to fastp's own).

```
reads ──► 1. bait ──────► 2. filter ──────► 3. assemble ──────────────────► 4. compare
          minimap2 -x sr   fastp, then       templated: samtools consensus    SNPs: minimap2, SKA2 or Parsnp
          or BBDuk         even subsampling  de novo:   SPAdes or SKESA       tree: FastTree or IQ-TREE
```

RASHER was designed for genome skimming: low-coverage sequencing of total DNA, where the organelle reads are
abundant enough to assemble. Enriching the target (targeted sequencing, hybridization capture) gives more depth
and better assemblies.

A rasher is a thin slice of bacon: RASHER is the short-read fork of [BACoN](https://github.com/duceppemo/BACoN),
which does the same for Nanopore reads. The comparison, the report, the annotation of SNP effects, the checkpoints
and most of the code are BACoN's; the reads (pairing, baiting, filtering) and the assemblers are RASHER's own.
[Methods](Methods#changes-from-bacon) lists the differences.

## Quick start

```bash
git clone https://github.com/duceppemo/RASHER
cd RASHER
conda env create -f environment.yml && conda activate RASHER
pip install .

rasher -r chloroplast.fasta -i reads/ -o results/
```

`reads/` holds the fastq files: `sample_R1.fastq.gz` and `sample_R2.fastq.gz` (or `_1`/`_2`, `.R1`/`.R2`) make
one paired sample, `sample`; a file without a mate is a single-end sample; a subfolder is a sample named after the
subfolder ([Usage](Usage#inputs)). Open `results/report.html`.

| Page | Contents |
|---|---|
| [Installation](Installation) | conda environment, pip, checking the installation |
| [Usage](Usage) | inputs and read pairing, all options, resuming, performance |
| [Methods](Methods) | what each step does, choosing an assembler and a SNP method, limits, references |
| [Example](Example) | a small simulated dataset with known SNPs, and its report |
| [Tutorial](Tutorial) | a complete analysis of public data |
| [Outputs](Outputs) | every file and column |
| [Validation](Validation) | how the tools were chosen and tested |
| [FAQ](FAQ) | troubleshooting |
| [Development](Development) | tests, continuous integration, releases |

[Releases](https://github.com/duceppemo/RASHER/releases) ·
[Changelog](https://github.com/duceppemo/RASHER/blob/main/CHANGELOG.md) ·
[Issues](https://github.com/duceppemo/RASHER/issues) · [BACoN](https://github.com/duceppemo/BACoN)
