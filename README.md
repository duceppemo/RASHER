# RASHER

<p align="center">
  <a href="https://github.com/duceppemo/RASHER/actions/workflows/ci.yml"><img src="https://github.com/duceppemo/RASHER/actions/workflows/ci.yml/badge.svg" alt="CI"></a>
  <a href="https://codecov.io/gh/duceppemo/RASHER"><img src="https://codecov.io/gh/duceppemo/RASHER/graph/badge.svg" alt="Coverage"></a>
  <img src="https://img.shields.io/badge/python-3.10%2B-blue" alt="Python 3.10+">
  <a href="LICENSE"><img src="https://img.shields.io/github/license/duceppemo/RASHER" alt="License: MIT"></a>
  <a href="https://github.com/duceppemo/RASHER/wiki"><img src="https://img.shields.io/badge/docs-wiki-informational" alt="Documentation"></a>
</p>

RASHER extracts the Illumina reads that match a reference sequence (typically an organelle genome) from many
samples, assembles them, and compares the samples: SNP distances and a tree. It was designed for genome
skimming, and works best with targeted sequencing that enriches the region of interest.

A rasher is a thin slice of bacon: RASHER is the short-read fork of
[BACoN](https://github.com/duceppemo/BACoN), which does the same for Nanopore reads.

```
reads ──► 1. bait (minimap2) ──► 2. filter (fastp) ──► 3. assemble ──────────────► 4. compare
                                                       templated: samtools consensus    SNPs: minimap2, SKA2 or Parsnp
                                                       de novo:   SPAdes or SKESA       tree: FastTree or IQ-TREE
```

## Quick start

```bash
git clone https://github.com/duceppemo/RASHER
cd RASHER
mamba env create -f environment.yml && conda activate RASHER
pip install .

rasher -r chloroplast.fasta -i reads/ -o results/
```

`reads/` holds the fastq files: `sample_R1.fastq.gz` and `sample_R2.fastq.gz` (or `_1`/`_2`) are paired
automatically into one sample; a file without a mate is a single-end sample; a subfolder is a sample named
after the subfolder with all the files it contains. Gzipped or not.

Open `results/report.html` for an overview: samples with their read counts, depth and assembly, the tree, a
heatmap of the SNP distances, and a methods paragraph. The same results are in `summary.tsv`,
`4_compared/minimap2/snp_distances.tsv`, `4_compared/minimap2/tree.nwk` and `4_compared/minimap2/snps.vcf`,
and MultiQC picks up RASHER's `*_mqc.json` files. Rerunning the same command resumes where it stopped;
changing a parameter reruns only the steps it affects. With `--annotation chloroplast.gb` (or a GenBank file
as `-r`), the report also shows the genes and the effect of each SNP; with `--metadata samples.tsv` (a table
with a `sample` column and any others, such as a group or a year), it shows the metadata in the samples table
and colours every genome in the report (tree, heatmap, bar charts, tables) by one column (`--color-by`).

To check an installation, run the example (simulated paired reads with known SNPs; about ten seconds). It is
in the repository:

```bash
curl -sL https://github.com/duceppemo/RASHER/archive/refs/tags/v0.1.0.tar.gz | tar -xz --strip-components=1 RASHER-0.1.0/example
bash example/run_example.sh
```

## Documentation

Everything else is in the [wiki](https://github.com/duceppemo/RASHER/wiki), whose sources are maintained in
[`docs/wiki`](docs/wiki):

- [Installation](https://github.com/duceppemo/RASHER/wiki/Installation) — conda, pip, troubleshooting
- [Usage](https://github.com/duceppemo/RASHER/wiki/Usage) — inputs, read pairing, options, resuming, performance
- [Methods](https://github.com/duceppemo/RASHER/wiki/Methods) — what each step does, and which assembler or SNP method to choose
- [Example](https://github.com/duceppemo/RASHER/wiki/Example) — the bundled simulated dataset
- [Outputs](https://github.com/duceppemo/RASHER/wiki/Outputs) — files and columns
- [Validation](https://github.com/duceppemo/RASHER/wiki/Validation) — how the methods were chosen and tested
- [FAQ](https://github.com/duceppemo/RASHER/wiki/FAQ) — troubleshooting
- [Development](https://github.com/duceppemo/RASHER/wiki/Development) — tests, CI, releases

## Citation

If RASHER is useful in your work, please cite it (see [`CITATION.cff`](CITATION.cff)):

> Duceppe, M.-O. (2026). RASHER: Bait, Assemble and Compare Illumina short reads (v0.1.0).

and the programs it runs: minimap2, fastp, samtools or SPAdes or SKESA, SKA2 or Parsnp, and FastTree or
IQ-TREE
([references](https://github.com/duceppemo/RASHER/wiki/Methods#references)). RASHER is a fork of
[BACoN](https://github.com/duceppemo/BACoN) (Duceppe, 2026, doi:10.5281/zenodo.22970412).

## License

[MIT](LICENSE)
