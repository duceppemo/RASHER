# Validation

RASHER's methods are chosen and tested against simulated datasets with a known truth and against public real
data. The validation suite lives in
[`validation/`](https://github.com/duceppemo/RASHER/blob/main/validation/README.md): the simulator
(plastid, linear, rdna, heteroplasmy and NUMT/NUPT scenarios as paired 2x150 bp reads), the runner and the
scoring scripts, with a dated record of every validation under `validation/results/`.

## Simulated scenarios (0.1.0)

Record:
[`validation/results/2026-10-22_v0.1.0/`](https://github.com/duceppemo/RASHER/blob/main/validation/results/2026-10-22_v0.1.0/SUMMARY.md).
Pairs of genomes with exactly the true SNP distance:

| Scenario | samtools + minimap2 (default) | SPAdes + SKA2 (de novo default) | SKESA + SKA2 |
|---|---|---|---|
| plastid | 36/45 | 36/45 | 36/45 |
| linear | 10/10 | 10/10 | 10/10 |
| rdna (6x tandem array) | **10/10** | 10/10 | 10/10 |
| heteroplasmy | 36/45 | **45/45** | 36/45 |
| NUMT/NUPT | 3/21 (all too low) | **21/21** | 7/21 |

What the record establishes:

- The templated consensus is error-free where 150 bp pairs map uniquely, and says **N, never a wrong base**,
  where they do not (the one exception is a 50/50 mixture, where either allele may be picked): over a long inverted repeat at low depth (reads there map with MAPQ 0), across a
  deletion longer than the insert, at mixed (heteroplasmic) sites from a 5% minor fraction, and at NUMT/NUPT
  sites. The minor allele of a mixture is never reported below a 50/50 tie.
- NUPT reads are baited with the organelle reads (up to 16% of the filtered bases at 25x nuclear depth), and
  their alleles never reach the consensus — but the resulting Ns hide the sample's own SNPs inside copied
  segments, biasing templated distances **low**. Where NUMTs/NUPTs are expected, use `-a spades`: de novo
  assembly separates the copies into their own contigs and SKA2's distances were exact.
- SPAdes `--isolate` assembled every clean sample without errors; it collapses the inverted repeat and
  tandem arrays to one copy (the templated assembly keeps them apart).
- SKESA behaved as its design suggests — it breaks at repeats: the plastome IR collapses, about 20 bp are
  lost at each contig end, and 20 NUPT sites go missing in one sample. Its contigs are error-free on clean
  samples but not under contamination or strong mixtures (up to 257 mismatches at 25x nuclear depth). With
  its default SKA2 pairing it beat SPAdes on no scenario. It stays available as `-a skesa` while the
  real-data benchmark runs; expect it to be dropped unless that overturns the simulations.

## Real-data benchmark (0.1.0)

Record:
[`validation/results/2026-10-23_benchmark_short/`](https://github.com/duceppemo/RASHER/blob/main/validation/results/2026-10-23_benchmark_short/SUMMARY.md).
Public Illumina data — 24 whole-genome *Arabidopsis thaliana* accessions (plastome and the repeat-rich
367 kb mitogenome) and 32 human 1000 Genomes mitochondria — against GetOrganelle, NOVOPlasty and Snippy on
the same reads, and against BACoN on matched Nanopore runs of the same accessions. 336/336 RASHER sample
runs succeeded. What the record establishes:

- **The default (templated + minimap2) reproduces Snippy**, an independent caller on the same reads:
  r = 0.998 on human mt (mean distance difference 1.9 SNPs), r = 0.976 on the plastome. The templated
  consensus was also the only full-length assembly on all three genomes (both plastome IR copies, the
  whole mitogenome, single contig each).
- **Illumina sees SNPs that a Nanopore consensus masks.** Against the matched BACoN-ONT runs, RASHER's
  distances are consistently higher (never lower, 274/276 pairs), and Snippy sides with RASHER — the ONT
  runs average 40 SNPs low per pair (r = 0.84 to Snippy).
- **SKESA, not SPAdes, is the cleaner de novo assembler on real data** — the simulations' verdict is
  overturned. Same SNP distances (r ≥ 0.999) and the same reference coverage, but a median of 2 contigs vs
  SPAdes' 22 on the plastome, and 32 contigs at reference length vs SPAdes' 205 contigs at 24% over on the
  mitogenome. The one place SPAdes wins is very low coverage: below ~20–30x SKESA's assemblies and
  distances degrade first.
- **Against the competition**: GetOrganelle is the stronger pure assembler for one complete circular
  organelle (55/56 complete, ~7 min per WGS sample; it extends beyond the bait). NOVOPlasty handled the
  small human mt (32/32 circular) but failed or stopped short on 9/24 plastome samples at ~62 min each.
  Snippy matches RASHER's SNPs but assembles nothing. RASHER (~11 min per WGS sample) is the only tool of
  the four that goes from raw reads to assemblies, SNP distances, a tree and a report in one run.

The parent tool's long-read validation is documented in
[BACoN's wiki](https://github.com/duceppemo/BACoN/wiki/Validation); its records remain under
`validation/results/` with their original BACoN names. Those results are about Nanopore reads and do not
carry over to RASHER.
