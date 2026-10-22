# Usage

```bash
rasher -r REFERENCE.fasta -i READS -o OUTPUT [options]
```

## Inputs

**Reference** (`-r`): a fasta file, gzipped or not, with one or more sequences, for example a complete
chloroplast genome. It is used to bait the reads, as the template of the templated assembly, and as the
anchor of the SNP comparison. A reference of the same species works best; a closely related species works for
baiting and de novo assembly. RASHER uses an upper-case copy (soft-masked, lower-case regions are compared like
the rest) in which IUPAC ambiguity codes are `N`; the same applies to the assemblies and to `--add-genomes`.

The reference can also be a GenBank file (`.gb`, `.gbk`, `.gbff`, `.genbank`, gzipped or not), such as NCBI's
record of a plastome: its sequence (`ORIGIN`) is the reference and its features annotate the report. Each
record is named after its `VERSION` (`NC_008096.2`), or its `LOCUS` name without one, followed by its
`DEFINITION`, as in NCBI's fasta of the same record: `reference.fasta` in the output folder is then identical to
the one made from that fasta, so switching between the two reruns nothing. `run_info.json` records the MD5 of
the file given, fasta or GenBank.

**Annotation** (`--annotation`): a GenBank or GFF3 file (`.gff`, `.gff3`; gzipped or not) describing the
reference, used in the report only: the genome map gets the genes, the regions of a plastome and the effect of
each SNP ([Outputs](Outputs#reporthtml)). It is not part of any checkpoint: adding, changing or removing it
reruns nothing. The LSC/IRb/SSC/IRa band of a plastome comes from the annotated inverted repeats or, when the
annotation does not mark them (older RefSeq plastomes), from the large inverted repeat RASHER detects in the
reference sequence itself, when the regions have the layout of a plastome (the repeats at least 5% of the
sequence, the larger single-copy region at most 200 kb, the record not a mitochondrion or a chromosome); the
band is drawn with or without an annotation, the gene track only with one
([Methods](Methods#5-snp-effects)). The default is the reference itself when it is a GenBank file. The annotated sequences must
have the names of the reference's sequences (a single annotated sequence of the same length as a single
reference sequence is accepted whatever its name); features on other sequences or beyond the end of a sequence
are ignored with a warning, and an annotated sequence whose length (`LOCUS`, `##sequence-region`) is not the
reference's gets a warning asking whether it is the annotation of this reference. RASHER copies the file,
uncompressed, to `OUTPUT/annotation.gb` or `OUTPUT/annotation.gff3`, so that the report can be rebuilt after the
folder is moved. What is read from the annotation, and how SNP effects are derived, is in
[Methods](Methods#5-snp-effects).

**Sample metadata** (`--metadata`): a TSV or CSV file with a `sample` column (any case) naming the samples and
any other columns, such as a group, a cultivar, a site or a year, used in the report only
([Outputs](Outputs#reporthtml)): the columns are added to the samples table, and one of them colours the bar
charts, the tree and the heatmap (each sample has the same colour and shape everywhere in the report), with a
table of the groups of identical genomes against its values. Like the annotation, it is
not part of any checkpoint: adding, changing or removing it reruns nothing, and the report is rebuilt. The
columns of a sample sheet other than `sample` and `file` are metadata too (a sample on several rows must have the
same values; otherwise the first is kept, with a warning); when both are given, a column of `--metadata` replaces
the sheet's column of the same name (whatever its case), and its columns come first. Names and values are
stripped, and runs of whitespace inside a metadata value become one space (a sample sheet's file paths are kept
as written); empty cells, `NA`, `na` and `-` are missing values. Two columns of the same name, whatever their
case, are an error; a column named like one of the samples table's own columns (`Status`, `Note`, `Depth`...)
is shown as `Status (metadata)`. A sample listed twice keeps its first row (warning); rows naming no sample of
the run are reported; samples without a row get blank cells. Genomes given with `--add-genomes` may have a row
too, under their file name. RASHER writes the merged table, limited to the run's samples and added genomes, to
`OUTPUT/metadata.tsv`, so that the report can be rebuilt after the folder is moved; a later run without metadata
removes that copy (a `metadata.tsv` that RASHER did not write is kept, and the report uses it, saying that it was
found in the folder and not given to the run; `run_info.json` and the MultiQC table do not have it).

The file format is the same for a metadata file and a sample sheet, with one difference: blank lines are
ignored, and lines starting with `#` before the header are comments; after the header, a line starting with `#`
is a comment in a sample sheet too (a sample is left out by commenting its line) but a row in a metadata file
(so a value such as `#FF0000` may come first). A tab in the header makes the file a TSV: a tab always separates
cells, a cell entirely in quotes loses them (`"abc"` is `abc`, and `""` inside such a cell is one quote, as a
spreadsheet exports them) and other quotes are ordinary characters (`5" tube`); otherwise it is a CSV, whose
quoted values may contain commas but not line breaks (an unclosed quote is an error naming the line, not a value
that silently swallows the following rows).

```
sample      group      year  comment
alpha       A          2021  Identical to the reference by design
beta        A          2021  Carries six planted SNPs
gamma       B          2022  Shares beta's six SNPs, plus four of its own
delta       NA         2023  Ten SNPs of its own; group unknown
```

`--color-by COLUMN` chooses the column that colours every genome in the report (tree, heatmap, bar charts,
tables); `--color-by none` leaves them uncoloured (and needs no metadata). By
default it is the first column that can be coloured: a column with at most 48 distinct values that does not
look like free text, that is, whose distinct values are at most 30 characters long on average and, once 10 or
more samples have a value, do not outnumber half of them (so 28 samples can have up to 14 values, and 48 values
need 96 samples). Numbers with few distinct values (years, batches) are categories like any other. A requested
column that does not exist stops RASHER before any step; one that cannot be coloured gives a warning and
uncoloured figures. The values are marked in sorted order (numerically when they are all numbers), each with a
colour and a shape of its own: 12 colours from the colourblind-safe palettes of Okabe & Ito and of Paul Tol,
chosen and ordered to stay as distinct as possible under deuteranopia, protanopia and tritanopia (the first ones,
which a column with few values gets, the most distinct), and seven shapes (circle, triangle, square, diamond,
inverted triangle, plus, cross). The first 12 values take the colours in turn with the first four shapes in turn,
so that the values sharing a shape are four colours apart and stay distinct to a colourblind reader; the next 36
were chosen by a search to keep the values sharing a shape as distinct as possible. Up to 25 values, values
whose colours look alike to a colourblind reader therefore differ in shape; with more, a few pairs sharing a
shape may look alike (more so towards 48 values), and the legend and the value written after each name in the
tree and in the tables tell them apart. Consecutive values differ in both colour and shape. The markers follow
the sorted values, so adding or removing a value (a new sample, an edited sheet) can change the colour and shape
of the values after it. A genome without a value gets a hollow circle.

**Reads**: Illumina short reads, paired or single-end, as fastq or fasta, gzipped or not (`.fastq`, `.fq`,
`.fasta`, `.fa`, `.fna`, `.fas`, with or without `.gz`). Three ways to give them:

| Input | Samples |
|---|---|
| `-i sample.fastq.gz` | one single-end sample, named after the file (`sample`) |
| `-i folder/` | two files of the folder that differ only by a mate marker are one paired sample, named by what comes before the marker (below); every other sequence file directly in the folder is a single-end sample named after the file; each subfolder is a sample named after the subfolder, made of all the sequence files it contains, in any depth (linked subfolders included; a subfolder without any is skipped, with a message), its pairs found the same way. Subfolders named `unclassified` and `mixed` (reads without a barcode, as demultiplexers name them) are skipped |
| `--sample-sheet samples.tsv` | a TSV or CSV file with the columns `sample` and `file`; several files per sample separated by `;` or on several rows (the R1 and R2 files of a sample, several lanes); relative paths start from the sheet's folder; blank lines and lines starting with `#` are ignored; any other column is sample metadata (above) |

**Read pairing.** Two files make a pair when they sit in the same folder and their names differ only by the mate
marker just before the extension:

| Mates | Sample |
|---|---|
| `X_R1.fastq.gz`, `X_R2.fastq.gz` (or `_r1`/`_r2`) | `X` |
| `X_1.fastq.gz`, `X_2.fastq.gz` | `X` |
| `X.R1.fastq.gz`, `X.R2.fastq.gz` | `X` |
| `X_S1_L001_R1_001.fastq.gz`, `X_S1_L001_R2_001.fastq.gz` (bcl2fastq) | `X_S1_L001` |

The `R` may be either case, the `_001` tail of bcl2fastq and BCL Convert is allowed, and the two mates must use
the same marker (`X_1` does not pair with `X_R2`). A file whose mate is missing, or a marker claimed by two files,
is not paired: in a folder, each such file is a single-end sample named after its whole stem (`X_R1`). Within one
sample (a subfolder, or a sample's files in a sheet), pairs and unpaired files do not mix: the files without a
mate next to the sample's pairs are skipped, with a warning naming them (list them as a sample of their own to
use them). A sample with several pairs (lanes, runs) keeps each pair as a pair. In a flat folder, each lane of a
bcl2fastq run is its own sample (`X_S1_L001`, `X_S1_L002`): put a sample's lanes in a subfolder, or list them in
a sample sheet, to make one sample of them.

The R1 and R2 files of a pair must list the mates in the same order, as sequencers and the SRA write them: the
mates may have identical names or end in `/1` and `/2`; a pair whose files desynchronise (different names at the
same rank, or a different number of reads) is an error. Each pair travels whole through the pipeline
(`<sample>_R1` and `<sample>_R2` files in `1_extracted/` and `2_filtered/`).

Sample names may contain letters, digits and `.` `_` `+` `-`. Two inputs giving the same sample name, a sample
mixing fasta and fastq files, and a paired sample `X` next to a sample named `X_R1` or `X_R2` (whose files would
have the same names) are errors. Symbolic links are a quick way to rename samples.

## Options

| Option | Default | Description |
|---|---|---|
| `-r`, `--reference` | | Reference fasta, or GenBank (its sequence is used and its features annotate the report); required |
| `-i`, `--input` | | Reads file or folder (this or `--sample-sheet`) |
| `--sample-sheet` | | TSV/CSV with `sample` and `file` columns; other columns are metadata |
| `-o`, `--output` | | Output folder (required) |
| `--annotation` | the reference, if GenBank | GenBank or GFF3 annotation of the reference, for the report only (genes, regions, SNP effects); never reruns a step |
| `--metadata` | | TSV/CSV with a `sample` column and any other columns, for the report only (samples table, colours of every genome in the report); never reruns a step |
| `--color-by` | first usable column | Colours every genome in the report (tree, heatmap, bar charts, tables) by this metadata column (at most 48 distinct values, each with a colour and a shape; not free text), or `none` |
| `-b`, `--baiting-method` | `minimap2` | `minimap2`: reads with an alignment to the reference; `bbduk`: reads sharing a k-mer (with `--hdist` mismatches) |
| `-k`, `--kmer-size` | 31 | BBDuk k-mer size (at most 31) |
| `--hdist` | 1 | BBDuk: mismatches allowed in a k-mer (0, 1 or 2). Each one multiplies BBDuk's memory: with 2, a 155 kb plastome needs about 14 GB per sample (ignored, with a warning, with `-b minimap2`) |
| `--keep-bam` | off | Keep the sorted BAM of the baited reads aligned to the reference (`1_extracted/<sample>.bam`; for paired samples, the baited pairs aligned as pairs). minimap2 baiting only: an error with `-b bbduk` |
| `--min-read-length` | 50 | Reads shorter than this after fastp's trimming are discarded; of a pair, both mates go when either is shorter |
| `--target-depth` | 100 | After fastp, the reads (pairs kept whole) are subsampled evenly to at most this depth (depth = bases / genome size, `-s` or the reference length); see [Methods](Methods#2-filtering) |
| `-a`, `--assembly-method` | `samtools` | `samtools` (templated), `spades` (de novo, SPAdes `--isolate`) or `skesa` (de novo, SKESA; not in the conda environment: `conda install -c bioconda skesa`); see [Methods](Methods#3-assembly) |
| `--template-gaps` | `n` | Templated assembly: what to put where fewer than three reads cover the reference. `n`: an `N` (no base is guessed); `reference`: the reference's bases, at the ends of each sequence only (gaps inside stay `N`) |
| `-s`, `--size` | reference length | Expected genome size, for the estimated depth and `--target-depth` |
| `--snp-method`, `-snp` | `auto` | `auto`: `minimap2` for templated assemblies, `ska` for SPAdes and SKESA. `minimap2` (each assembly aligned to the reference, read position by position), `ska` (SKA2 split k-mers), `parsnp` (core-genome alignment), `none` (stop after the assembly); see [Methods](Methods#4-comparison) |
| `--ska-min-freq` | 1.0 | SKA2: fraction of the genomes that must contain a variant's context; minimap2: fraction that must have a base at its position; 1 = core SNPs, lower = pan-genome SNPs (like kSNP); see [Methods](Methods#4-comparison) |
| `--add-genomes` | | Finished genomes (fasta) to include in the comparison, such as published plastomes; named after their file |
| `--tree` | `fasttree` | `fasttree` (GTR, SH-like supports from 100 resamples) or `iqtree` (model selection, 1000 ultrafast bootstraps) |
| `--redo` | | Rerun this step and the following ones: `bait`, `filter`, `assemble`, `compare` |
| `-t`, `--threads` | all | Total threads, shared between the samples processed in parallel |
| `-p`, `--parallel` | 2 | Samples processed at the same time |
| `-m`, `--memory` | 85% of RAM | Total memory for BBDuk, SPAdes and SKESA, in GB, divided between the samples processed in parallel (`-p`) |
| `--debug` | | Verbose log |
| `-v`, `--version` | | Print RASHER's version and exit |
| `-h`, `--help` | | Print the options and exit |

## Resuming and changing parameters

Each step records the parameters it ran with in `OUTPUT/.checkpoints/`. Running the same command again skips
every finished step; changing a parameter reruns that step and the ones after it, and only them
(`--template-gaps` with a de novo assembler reruns nothing, and `--keep-bam` only adds BAMs: it reruns only the
samples without one). For example, after a first run with the defaults:

```bash
rasher -r ref.fasta -i reads/ -o out/ -a spades                        # reuses baiting and filtering, assembles with SPAdes
rasher -r ref.fasta -i reads/ -o out/ -a spades --snp-method parsnp    # reuses the SPAdes assemblies
```

Each comparison is written to its own folder (`4_compared/minimap2/`, `4_compared/ska/`, `4_compared/parsnp/`,
`4_compared/ska_0.5/`), so several can be kept side by side. A sample that failed is retried on the next run;
samples added to the input are processed without redoing the others' baiting and filtering (their
assemblies are compared again). The files of a sample removed from the input (or whose line of the sample
sheet is commented out) stay in the output folder, unused, and are reused if it comes back unchanged; after
switching from a de novo to the templated assembly, the de novo assembly graphs of each sample are
removed.

`--redo STEP` forces a step to run again, for example after installing a newer assembler.

A run stopped with Ctrl-C (or SIGTERM, SIGHUP: `kill`, a closed terminal, a job scheduler's time limit) stops
the programs it started; each sample finished before the interruption is recorded, so the next run redoes only
the samples that were still running. A run started with `nohup` keeps running when the terminal is closed. Two
runs cannot use the same output folder at the same time: the second one stops with an error.

A run killed with SIGKILL (`kill -9`, the kernel's out-of-memory killer) cannot stop the programs it started:
each runs in its own process group, and they keep running (and writing to the output folder) until they end.
Their command lines name the output folder: list them with `pgrep -af /full/path/to/out` (or
`ps -ef | grep /full/path/to/out`) and, once the list holds only them, stop them with `pkill -f /full/path/to/out`
before running RASHER again on that folder.

An input file in the output folder that RASHER would overwrite (`-r`, `--annotation`, `--metadata`,
`--sample-sheet` named as one of the files RASHER writes there, e.g. `reference.fasta` or `metadata.tsv`) is an
error, unless it is already exactly the copy RASHER would write (`-r out/reference.fasta` of an earlier run).

An output folder can be moved or copied and resumed from its new place, whether or not the copy kept the files'
times (`cp -r`, `scp`, `rsync` without `-t`, an archive): RASHER recognizes each step's input by its size and
time in the folder where it was made, and by its size alone in a moved or copied folder. An output folder of
BACoN does not resume under RASHER: every step runs again.

## Performance

Baiting reads every input read once (minimap2 aligns each file once); the other steps work on the baited reads
only, capped at `--target-depth`. With whole-genome reads, baiting therefore takes most of the time of a
templated run; a de novo assembly adds SPAdes' or SKESA's time per sample. `-t` threads are shared between the
`-p` samples run at the same time, and `-m` memory between the samples that BBDuk, SPAdes or SKESA process at the
same time. Timings and memory on real data will be given with the short-read [validation](Validation) and
[tutorial](Tutorial).
