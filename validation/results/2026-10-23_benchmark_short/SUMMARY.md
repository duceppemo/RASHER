# Real-data benchmark of RASHER 0.1.0 (2026-10-23/24, pre-release)

Real-data benchmark of the short-read fork: RASHER's three assemblers against each other, against the
matched BACoN-ONT runs, and against GetOrganelle, NOVOPlasty and Snippy on the same reads. Stage 4
(human-mt SPAdes, the Arabidopsis *mitogenome*, a depth sweep) was added mid-benchmark before concluding
on SKESA vs SPAdes. The run scripts are in [`scripts/`](scripts/), the matrix scorer is
[`score_acc.py`](score_acc.py).

Run from RASHER 0.1.0 (commit `decfcce`, unreleased), env `RASHER_dev`, one machine: 2 Intel Xeon
Gold 6226R (64 threads, runs used `-t 24`), 503 GB RAM, Linux. No RASHER run failed: 336/336 sample runs
OK across 13 legs. Stage 4 ran alongside the (single-threaded) NOVOPlasty leg of stage 3 on the otherwise
idle machine, so its wall times carry a small contention penalty; stages 1–2 had the machine to
themselves. Helper scripts: [`aln2dist.py`](aln2dist.py) (snippy-core alignment → distance matrix),
[`comp_score.py`](comp_score.py) (competitor assemblies → reference coverage).

## Datasets (public, no simulation)

- **Arabidopsis accessions** (`acc_cp_*`): 24 PRJEB62038 accessions, HiSeq 2x100 PE WGS (0.46–1.1 Gbases
  baited of 11–16 Gbases raw per sample), reference plastome NC_000932.1 (154,478 bp). Matched ONT runs of
  the same accessions exist as BACoN 0.4.0 results (templated and Flye, SKA2 distances) for concordance.
- **Human mitochondria** (`mt_*`): 32 1KGP chrM Illumina BAM slices converted to paired fastq
  (`samtools collate | samtools fastq`), reference rCRS NC_012920.1 (16,569 bp).
- **Arabidopsis mitogenome** (`acc_mt_*`, stage 4): the same 24 accessions against NC_037304.1
  (367,808 bp) — large and repeat-rich, the hard de novo case the plastome is not.
- **Depth sweep** (`acc_cp_*_d30/_d15`, stage 4): the plastome runs repeated at `--target-depth` 30 and
  15 (the main runs cap at the default 100x), with the already-baited reads of the templated run as
  input.

## Per-leg results

| Leg | Samples OK | Contigs (median, range) | Length (median, range) | SNP sites | Wall | Peak RAM |
|---|---|---|---|---|---|---|
| mt templated (minimap2) | 32/32 | 1 | 16,569 (16,561–16,573) | 443 | 14:08 | 1.3 GB |
| mt SKESA (SKA2) | 32/32 | 1 | 16,568 (16,526–16,573) | 174 | 14:38 | 1.3 GB |
| acc plastome templated (minimap2) | 24/24 | 1 | 154,466 (154,421–154,518) | 503 | 4:23:10 | 4.3 GB |
| acc plastome SPAdes (SKA2) | 24/24 | 22 (13–31) | 133,419 (130,652–136,584) | 321 | 4:33:00 | 4.3 GB |
| acc plastome SKESA (SKA2) | 24/24 | 2 (1–3) | 128,178 (127,964–128,245) | 321 | 4:29:40 | 4.9 GB |
| mt SPAdes (SKA2) | 32/32 | 1 (1–3) | 16,646 (16,638–16,878) | 176 | 14:21 | 1.3 GB |
| acc mitogenome templated (minimap2) | 24/24 | 1 | 367,796 (367,759–367,813) | 216 | 4:16:22 | 1.3 GB |
| acc mitogenome SPAdes (SKA2) | 24/24 | 205 (176–243) | 455,928 (426,730–471,003) | 184 | 4:18:59 | 1.3 GB |
| acc mitogenome SKESA (SKA2) | 24/24 | 32 (25–40) | 366,010 (346,137–376,961) | 172 | 4:20:36 | 1.3 GB |
| acc plastome SPAdes 30x | 24/24 | 8 (5–11) | 129,288 | 321 | 1:53:11 | 7.5 GB |
| acc plastome SKESA 30x | 24/24 | 3 (1–4) | 128,122 | 320 | 1:51:56 | 7.5 GB |
| acc plastome SPAdes 15x | 24/24 | 5 (2–9) | 128,704 | 322 | 1:44:55 | 7.5 GB |
| acc plastome SKESA 15x | 24/24 | 6 (3–9) | 127,892 | 307 | 1:51:26 | 7.5 GB |

Wall time is dominated by baiting the full WGS runs (each leg re-baits in its own folder so the timings are
honest per-tool figures); the assembler choice changes the total by minutes. The depth-sweep legs bait the
already-baited reads (a ~20x smaller input where every read hits), hence the shorter walls and a different
memory profile.

## SNP-distance concordance, Arabidopsis plastome (276 sample pairs, `score_acc.py`)

| Comparison | Identical | Mean \|diff\| | Max \|diff\| | r |
|---|---|---|---|---|
| BACoN-ONT templated vs BACoN-ONT Flye (both SKA2) | 56/276 | 1.8 | 5 | 0.985 |
| RASHER SPAdes vs RASHER SKESA (both SKA2) | 233/276 | 0.2 | 2 | **0.9996** |
| RASHER SPAdes vs BACoN-ONT templated (SKA2 vs SKA2) | 2/276 | 14.1 | 26 | 0.932 |
| RASHER SKESA vs BACoN-ONT templated (SKA2 vs SKA2) | 2/276 | 14.3 | 26 | 0.931 |
| RASHER templated (minimap2) vs BACoN-ONT templated | 1/276 | 42.2 | 66 | 0.907 |
| RASHER templated (minimap2) vs RASHER SPAdes (SKA2) | 1/276 | 28.1 | 43 | 0.953 |

Median pairwise distances: BACoN-ONT 37–38, RASHER SPAdes/SKESA 53, RASHER templated (minimap2) 82.

## Interim findings

1. **RASHER never sees fewer SNPs than BACoN on the matched ONT data.** In the like-for-like SKA2
   comparison, RASHER SPAdes is higher in 274/276 pairs (median +14) and lower in none. The direction is
   what N-masking predicts — the ONT consensus masks ambiguous sites that accurate Illumina bases resolve —
   rather than Illumina noise, but the Snippy leg (stage 3) is the arbiter, as an independent caller on the
   same reads. The templated + minimap2 count is higher again (median +29 vs SPAdes + SKA2): whole-genome
   alignment also counts sites in regions SKA2's split k-mers skip (clustered variants, indel-adjacent).
2. **Real data reverses the simulated SKESA contiguity verdict.** The simulations (previous record) had
   SKESA fragmenting at repeats and never beating SPAdes; on the real plastome skims SKESA produced a
   median of 2 contigs against SPAdes' 22, with essentially identical distances (233/276 identical,
   max diff 2, r = 0.9996) and the same cost. Mapping both assemblies to the reference shows **identical
   unique content**: both cover 82.96% of the plastome — exactly the genome minus one inverted-repeat copy
   (both collapse the IR, as expected de novo) — and SPAdes' extra ~5 kb of total length is small
   off-target fragments (only 5–6 of its ~22 contigs align), not extra plastome. On the human mt dataset
   SKESA is single-contig full-length on 32/32 samples.
   **Decision deferred**: stage 4 adds the missing human-mt SPAdes leg, the repeat-rich Arabidopsis
   mitogenome (NC_037304.1, ~367 kb), and a 30x/15x depth sweep before the recommendation is fixed.
3. **The templated leg is the only full-length assembly** (154.5 kb with both IR copies, single contig,
   34 N bases typical) and the cheapest, consistent with the simulations.

## Stage 4: the SKESA-vs-SPAdes question on harder ground

4. **Human mt confirms the plastome pattern.** SPAdes vs SKESA distances are interchangeable (435/496
   pairs identical, max diff 2, r = 0.9994), but SKESA's contigs stay at or under the rCRS length
   (16,526–16,573) while SPAdes inflates some samples to 16,878 and splits up to 3 contigs.
5. **The repeat-rich mitogenome separates them decisively.** The templated consensus is single-contig
   full-length (367.8 kb). De novo, both fragment at the repeats — but SKESA yields a median of
   32 contigs totalling ~366 kb (≈ the reference), while SPAdes yields a median of **205 contigs
   totalling ~456 kb (24% over the reference)**. Mapped to the reference, both cover the same 94–96%;
   SPAdes' extra 60–90 kb per sample is redundant or off-target sequence. Distances again agree
   (r = 0.993, mean diff 1.3). The simulations' "SKESA breaks at repeats" holds — but SPAdes breaks
   *more*, and pads the pieces with junk.
6. **Depth sweep: 30x is free, 15x is where SKESA starts to crack.** At 30x both match the 100x
   distances essentially perfectly (SKESA 100x vs SPAdes 30x: 276/276 identical) and SKESA stays more
   contiguous (3 vs 8 contigs). At 15x the roles flip mildly: SPAdes is slightly more contiguous
   (5 vs 6 contigs) and more stable (vs 100x: mean diff 0.7, max 2, 322 SNP sites) than SKESA
   (mean diff 2.2–2.4, max 5, 307 sites, shortest assemblies). The conservative extender gives up
   first when coverage thins.
7. **Where this leaves the de novo choice**: SKESA beat SPAdes
   on contiguity and off-target restraint on all three genomes at normal depth, with identical distances
   and identical cost — the simulation-based "drop SKESA" recommendation of the previous record is
   overturned, and SKESA now looks like the better de novo *recommendation*, with a documented caveat to
   prefer SPAdes (or the templated default) below ~20–30x. Both remain available either way; the
   templated consensus remains the overall default and the only full-length assembly everywhere.

## Stage 3: Snippy arbitration and the assembly competitors

8. **Snippy sides with RASHER's templated + minimap2 default, strongly.** Pairwise distances from
   `snippy-core`'s core SNP alignment (same reads, independent caller — bwa + freebayes):

   | Snippy vs | Identical | Mean \|diff\| | Max \|diff\| | r |
   |---|---|---|---|---|
   | RASHER templated (minimap2), human mt | 74/496 | 1.9 | 5 | **0.998** |
   | RASHER templated (minimap2), plastome | 23/276 | 3.8 | 13 | **0.976** |
   | RASHER SPAdes/SKESA (SKA2), both datasets | 1–2 | 26–28 | 48–62 | 0.92–0.93 |
   | BACoN-ONT templated (SKA2), plastome | 1/276 | 40.2 | 66 | 0.84 |

   Medians line up as RASHER templated ≈ Snippy (82 vs 81 on the plastome, 42 vs 40 on human mt), the
   SKA2 de novo legs sit lower (split k-mers skip clustered/indel-adjacent sites), and BACoN-ONT lowest
   (37). **The finding-1 offset is resolved: the SNPs RASHER sees and the ONT runs don't are real** —
   the ONT consensus masks them; Illumina reads and two independent Illumina callers agree on them.
9. **GetOrganelle is the assembly reference point, and the only competitor that matches RASHER's
   reliability.** 23/24 complete circular plastomes (median 1 contig, 154.5 kb, ≥ 99.9% reference
   coverage; the one failure, Ct-1, stops at 128 kb/83%) and 32/32 single-contig human mt at 99.5–99.9%
   coverage. Cost: 7.3 min per plastome sample, 42 s per mt sample.
10. **NOVOPlasty is fast only on small genomes and unreliable on WGS skims.** Human mt: 32/32
    circularized, ~5 min per sample. Plastome: only 15/24 assemblies reach 99.9% reference coverage —
    3 total failures (~300 bp output), 2 stop near 17%, and 4 "full-length" 154 kb assemblies align over
    only 83% of the reference — at 62.5 min per sample (single-threaded, up to 106 min).
11. **Runtime, per Arabidopsis WGS sample, like for like** (each tool reads the full ~12–16 Gbase run):
    GetOrganelle 7.3 min (assembly only), Snippy 9.1 min (SNPs only), RASHER templated 11.0 min
    (assembly + SNPs + tree + report), NOVOPlasty 62.5 min (assembly only, often incomplete). On the
    small human-mt inputs RASHER is ~26 s per sample, Snippy ~26 s, GetOrganelle ~42 s, NOVOPlasty ~5 min.
    RASHER's per-sample wall is dominated by baiting; everything downstream runs on ~20% of the data.

## Where this leaves the defaults (decision with the user)

- **Templated (`-a samtools`) + minimap2 stays the overall default**: the only full-length assemblies on
  all three genomes, distances near-identical to Snippy on both datasets (findings 3, 8), and the
  cheapest leg.
- **SKESA looks like the better de novo recommendation** (finding 7), with the ~20–30x caveat; SPAdes
  stays available (assembly graph, low-depth robustness).
- **Against the competition**: GetOrganelle remains the stronger pure assembler for a single complete
  circular organelle (it extends beyond the bait, where RASHER de novo assembles only baited reads);
  RASHER is the only tool of the four that goes from raw WGS to SNP distances, a tree and a report, and
  its SNP calls match Snippy without needing a second pipeline. NOVOPlasty is not a concern at WGS-skim
  depths.
- Whether broader evidence (other taxa/instruments) is needed before fixing the defaults is an open
  user decision.
