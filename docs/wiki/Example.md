# Example

The repository bundles a small simulated example: an annotated 30 kb circular reference and four samples of
paired 2x150 bp reads with known SNPs. It runs in about ten seconds and checks its own results, so it doubles
as an installation test.

## Running it

From a clone of the repository (the example is not part of a package install):

```bash
bash example/run_example.sh
```

or from a release download:

```bash
curl -sL https://github.com/duceppemo/RASHER/archive/refs/tags/v0.1.0.tar.gz | tar -xz --strip-components=1 RASHER-0.1.0/example
bash example/run_example.sh
```

`run_example.sh [OUTPUT_FOLDER] [THREADS]` does three things:

1. `make_example.py` simulates the data into `example_output/data/`: `reference.fasta`, its GenBank
   annotation `reference.gb` (synthetic genes, an inverted repeat), `metadata.tsv`, and
   `reads/<sample>_R1.fastq.gz` + `_R2.fastq.gz` for the samples `alpha`, `beta`, `gamma` and `delta`
   (about 50x each, 20% off-target pairs, 0.2% substitution errors). The simulation is deterministic: the
   same files every time (standard library only).
2. `rasher` runs with the default settings (minimap2 baiting, fastp, templated assembly, minimap2 SNPs,
   FastTree), plus `--annotation` and `--metadata`.
3. `check_example.py` compares the results with the planted truth and fails loudly on any difference.

## What is planted, and what is checked

- `alpha` matches the reference; `beta` carries 6 SNPs; `gamma` carries beta's 6 plus 4 of its own;
  `delta` carries 10 of its own. The pairwise SNP distances must come out exactly:
  alpha–beta 6, alpha–gamma 10, alpha–delta 10, beta–gamma 4, beta–delta 16, gamma–delta 20.
- The VCF must hold the 20 planted sites, and the report's SNP table their region, gene, codon and effect
  (among the 12 coding changes: 5 synonymous, 6 missense and 1 nonsense).
- The report must show the metadata columns and colour the genomes by group.

Open `example_output/rasher/report.html` to see the result (the report of this run is also
[online](https://duceppemo.github.io/RASHER/reports/example_report.html)): the samples table (read counts, depth,
assembly), the tree, the SNP distance heatmap, the genome map with the annotated genes and every SNP, and
the methods paragraph. The same numbers are in `summary.tsv` and `4_compared/minimap2/` — see
[Outputs](Outputs).

The reads are simulated the same way as the [validation](Validation) scenarios (`validation/simulate.py`),
with 20% off-target pairs instead of the validation's 15%; the example is just small enough to run anywhere.
