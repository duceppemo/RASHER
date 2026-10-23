# Real-data benchmark of RASHER 0.1.0 — INTERIM (stages 1–2 of 4)

First real-data results for the short-read fork. **Interim record**: the competitor legs (stage 3:
GetOrganelle, NOVOPlasty, Snippy) and the extra assembler legs (stage 4: human-mt SPAdes, the Arabidopsis
*mitogenome*, a depth sweep) were still running or queued when this was written; **no defaults or
SKESA-vs-SPAdes conclusion is drawn here** — more benchmarks were explicitly requested before concluding.
The run scripts are in [`scripts/`](scripts/), the matrix scorer is [`score_acc.py`](score_acc.py).

Run from RASHER 0.1.0 (commit `decfcce`, unreleased), env `RASHER_dev`, one machine: 2 Intel Xeon
Gold 6226R (64 threads, runs used `-t 24`), 503 GB RAM, Linux. No run failed: 136/136 sample runs OK.

## Datasets (public, no simulation)

- **Arabidopsis accessions** (`acc_cp_*`): 24 PRJEB62038 accessions, HiSeq 2x100 PE WGS (0.46–1.1 Gbases
  baited of 11–16 Gbases raw per sample), reference plastome NC_000932.1 (154,478 bp). Matched ONT runs of
  the same accessions exist as BACoN 0.4.0 results (templated and Flye, SKA2 distances) for concordance.
- **Human mitochondria** (`mt_*`): 32 1KGP chrM Illumina BAM slices converted to paired fastq
  (`samtools collate | samtools fastq`), reference rCRS NC_012920.1 (16,569 bp).

## Per-leg results

| Leg | Samples OK | Contigs (median, range) | Length (median, range) | SNP sites | Wall | Peak RAM |
|---|---|---|---|---|---|---|
| mt templated (minimap2) | 32/32 | 1 | 16,569 (16,561–16,573) | 443 | 14:08 | 1.3 GB |
| mt SKESA (SKA2) | 32/32 | 1 | 16,568 (16,526–16,573) | 174 | 14:38 | 1.3 GB |
| acc plastome templated (minimap2) | 24/24 | 1 | 154,466 (154,421–154,518) | 503 | 4:23:10 | 4.3 GB |
| acc plastome SPAdes (SKA2) | 24/24 | 22 (13–31) | 133,419 (130,652–136,584) | 321 | 4:33:00 | 4.3 GB |
| acc plastome SKESA (SKA2) | 24/24 | 2 (1–3) | 128,178 (127,964–128,245) | 321 | 4:29:40 | 4.9 GB |

Wall time is dominated by baiting the full WGS runs (each leg re-baits in its own folder so the timings are
honest per-tool figures); the assembler choice changes the total by minutes.

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

Stage 3 (GetOrganelle, NOVOPlasty, Snippy on both datasets) and stage 4 results, the assembly-completeness
scoring, and the defaults discussion will extend this record.
