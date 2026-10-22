# Installation

RASHER is a Python package (standard library only) that runs command-line programs, which are installed with
conda.

## From the source code (recommended for now)

```bash
git clone https://github.com/duceppemo/RASHER
cd RASHER
conda env create -f environment.yml
conda activate RASHER
pip install .
rasher --version
```

[`environment.yml`](https://github.com/duceppemo/RASHER/blob/main/environment.yml) lists the programs and their
minimum versions:

| Program | Used for | Minimum |
|---|---|---|
| minimap2 | baiting (`-x sr`); templated assembly; off-reference share of de novo assemblies; SNPs (`--snp-method minimap2`, the default for templated assemblies) | 2.24 |
| samtools | templated assembly (`samtools consensus -X hiseq`); `--keep-bam` | 1.21 |
| fastp | read filtering and adapter trimming | 1.0 |
| SPAdes | de novo assembly (`-a spades`) | 4.0 |
| SKA2 | SNPs (`--snp-method ska`; the default for de novo assemblies) | 0.5 |
| Parsnp, HarvestTools | SNPs (`--snp-method parsnp`) | 2.1.2 (with 2.1.1, every SNP is one base off) |
| FastTree | tree (default) | 2.1.11 |
| IQ-TREE | tree (`--tree iqtree`) | 2.2 |
| BBMap (BBDuk) | baiting (`-b bbduk`) | 39 |
| Bandage | optional: pictures of the SPAdes assembly graphs | |
| python-isal | optional: faster reading of gzipped reads | |

SKESA (`-a skesa`) is not in the environment: it is an optional assembler, still being evaluated
([Validation](Validation)). Install it into the environment to use it:

```bash
conda install -c conda-forge -c bioconda skesa
```

RASHER checks at start-up that the programs needed by the chosen options are on the `PATH`, and lists the
missing ones with the command to install them.

## bioconda (planned)

There is no bioconda package yet: it will be submitted with the first release
([`recipe/meta.yaml`](https://github.com/duceppemo/RASHER/blob/main/recipe/meta.yaml)). Once it is published,
this section will give the command, which is expected to be:

```bash
conda create -n rasher -c conda-forge -c bioconda rasher
conda activate rasher
rasher --version
```

It will install RASHER with every program above except Bandage and SKESA. Until then, install from the source
code.

## Checking the installation

The example is in the repository (`example/`). From a clone, or after downloading the example folder of the
release:

```bash
curl -sL https://github.com/duceppemo/RASHER/archive/refs/tags/v0.1.0.tar.gz | tar -xz --strip-components=1 RASHER-0.1.0/example
bash example/run_example.sh
```

generates a small simulated dataset (an annotated reference, a few samples of short reads with known SNPs, and
their metadata), runs RASHER with the default settings and checks the results against the truth: every pairwise
SNP distance, the region, gene, context and effect of every SNP, and the report's annotation and metadata
([Example](Example)); it prints one `OK` line per check.

## pip only

If the programs are already installed (for example in an HPC module system), RASHER itself installs with
`pip install https://github.com/duceppemo/RASHER/archive/refs/tags/v0.1.0.tar.gz` (or
`pip install git+https://github.com/duceppemo/RASHER` for the latest code).

## Troubleshooting

- **The environment takes long to solve, or mamba reports "nothing provides ..." for packages that exist**:
  some mamba 2.x versions misreport conflicts; `conda env create` (with the libmamba solver) solves it.
- **`samtools consensus: unrecognised option -X`**: samtools is older than 1.17 (which added `-X`); update
  it. RASHER is tested with samtools 1.21 and later.
- **Next to BACoN**: install RASHER in its own environment (BACoN's has Filtlong, Flye and myloasm, not fastp or
  SPAdes). The two commands do not share output folders: a BACoN output folder does not resume under RASHER
  (every step runs again).
