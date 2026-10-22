# Development

## Setup

```bash
git clone https://github.com/duceppemo/RASHER && cd RASHER
conda env create -f environment.yml && conda activate RASHER
pip install -e ".[test]"
pre-commit install
```

## Code

| Module | Contents |
|---|---|
| `rasher/cli.py` | command line |
| `rasher/pipeline.py` | the steps, checkpoints, summary and provenance |
| `rasher/samples.py` | finding samples (file, folder, sample sheet) |
| `rasher/metadata.py` | sample metadata: TSV/CSV reading (shared with the sample sheet), merging, the column that colours the report |
| `rasher/steps.py` | baiting, filtering, assembly |
| `rasher/compare.py` | SKA2, Parsnp, distances, trees |
| `rasher/newick.py` | Newick parsing, midpoint rooting, SVG drawing |
| `rasher/report.py` | the HTML report (`python -m rasher.report OUTPUT` rebuilds it) |
| `rasher/annotation.py` | GenBank and GFF3 annotations: genes, plastome regions, the effects of SNPs on coding sequences |
| `rasher/multiqc.py` | the MultiQC custom-content files |
| `rasher/seqio.py` | fasta/fastq reading and writing |
| `rasher/tools.py` | running external programs and logging their output |

RASHER uses the Python standard library only (python-isal is used when installed). External programs are run
with `rasher.tools.run`, which logs the command and the program's messages and raises an error with the end of
the log when a program fails.

## Tests

```bash
ruff check .
pytest --cov=rasher
```

The tests replace the external programs by small stub scripts (`tests/test_pipeline.py`), so they run in
seconds, without conda, and test the pipeline logic: resuming, parameter changes, failures, provenance. The
real programs are exercised by the example and the validation:

```bash
bash example/run_example.sh                   # seconds; checks SNP distances, effects and the report against the truth
bash validation/run_validation.sh /tmp/val    # a few minutes; scores every assembler and SNP method
```

Continuous integration (`.github/workflows/ci.yml`) runs ruff, the tests on Python 3.10 to 3.13 (and macOS),
and the example in the conda environment.

## Adding or replacing a program

A program is added or replaces another only if it is maintained, installs with the others from bioconda,
and does at least as well in the validation. Add it to `environment.yml`, `recipe/meta.yaml` and
`rasher/tools.py` (`PACKAGES`), with a stub in `tests/test_pipeline.py`, then run the validation and record the
results in a new dated folder under `validation/results/` (never rewrite an old record).

## Documentation

The wiki is maintained in `docs/wiki/` and published to the GitHub wiki by `.github/workflows/wiki.yml` on
every push to `main` that changes it. Do not edit the wiki on GitHub: the next publication overwrites it.

## Web site

GitHub Pages serves `docs/` at https://duceppemo.github.io/RASHER/ (`docs/.nojekyll`: files are served as they
are). It holds the example reports of `docs/reports/`; rebuild them after changes to the report, with local
paths removed from `run_info.json` first.

## Releases

1. Update `version` in `pyproject.toml`, `__version__` in `rasher/__init__.py`, `CITATION.cff` (version and
   date), `CHANGELOG.md`, the citation and the example's download line in `README.md`, and the download and
   `pip install` lines of `docs/wiki/Installation.md` (`tests/test_cli.py` checks that they agree). Then
   rebuild the reports of `docs/reports/`: they show the RASHER version of `run_info.json`, so run the example
   (`bash example/run_example.sh`) and the tutorial's command again with the bumped version (a resume reruns
   nothing; the report is rebuilt) before copying their `report.html`, with the local paths removed from
   `run_info.json` first (`tests/test_cli.py` checks that the published reports name the current version).
   Then take the wiki's images of the report again from the rebuilt reports (`docs/images/example_report_*.png`
   and `report_*.png`: the same parts of the page, cropped alike), so that they show the release. In
   `CHANGELOG.md`, rename `## Unreleased` to `## X.Y.Z (YYYY-MM-DD)` (`tests/test_cli.py` checks that the
   section of the current version exists; the release workflow publishes it).
2. Commit, tag `vX.Y.Z` and push the tag: `.github/workflows/release.yml` checks the versions, builds the
   package and creates the GitHub release with the changelog section.
3. Zenodo archives the release and mints a version DOI: add it to the `identifiers` of `CITATION.cff`. The
   top-level `doi` stays the concept DOI (10.5281/zenodo.22970412), which always points to the latest version:
   GitHub's "Cite this repository" at any tag then gives a DOI that resolves (a version DOI is only known after
   the tag). The README badge uses the concept DOI too.
4. Update `recipe/meta.yaml` (version, sha256 of the tag's tarball, build number 0). Bioconda's autobump bot
   usually opens the update pull request in bioconda-recipes by itself; otherwise open one with the new recipe.
