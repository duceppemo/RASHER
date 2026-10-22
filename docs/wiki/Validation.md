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

## Real-data benchmark

In progress: public Illumina data (whole-genome *Arabidopsis thaliana* accessions, human 1000 Genomes
mitochondria) against GetOrganelle, NOVOPlasty and Snippy, and against BACoN on matched Nanopore runs of the
same samples. It will decide the final defaults and add timings; its record will be summarised here.

The parent tool's long-read validation is documented in
[BACoN's wiki](https://github.com/duceppemo/BACoN/wiki/Validation); its records remain under
`validation/results/` with their original BACoN names. Those results are about Nanopore reads and do not
carry over to RASHER.
