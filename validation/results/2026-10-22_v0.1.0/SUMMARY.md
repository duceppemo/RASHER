# Validation of RASHER 0.1.0 (2026-10-22, pre-release)

The first validation of the short-read fork: every scenario of `validation/simulate.py`, now simulating
Illumina 2x150 bp pairs (inserts 350 ± 50 bp, 0.2% substitutions — see `validation/README.md`), run with the
three assemblers (`samtools` templated, SPAdes, SKESA) and the SNP methods. **The read model changed with the
fork, so none of these numbers is comparable to the BACoN records in the sibling folders**, which simulated
R10.4.1 Nanopore reads; what carries over is the truth (the genomes and planted variants are unchanged).

Run from the working tree at commit `4b9074b` (RASHER 0.1.0 unreleased), env `RASHER_dev`, one machine:
2 Intel Xeon Gold 6226R (64 threads), 503 GB RAM, Linux.

```bash
bash validation/run_validation.sh <folder> 16       # plastid, linear, rdna
bash validation/run_validation.sh <folder> 16 new   # heteroplasmy, numt
```

Program versions: [`versions.txt`](versions.txt) (minimap2 2.31, fastp 1.1.0, samtools 1.24, SPAdes 4.2.0,
SKESA 2.5.1, SKA2 0.5.1, Parsnp 2.1.5, FastTree 2.2.0). Scores: `SCORES_<scenario>.md` (and `.json`).

## Exact SNP distances (pairs with exactly the true distance, best convention)

| Scenario | samtools + minimap2 (default) | samtools + SKA2 | SPAdes + SKA2 (default) | SPAdes + minimap2 | SKESA + SKA2 |
|---|---|---|---|---|---|
| plastid | 36/45 (max err 3) | 28/45 | 36/45 (max err 1) | **45/45** | 36/45 |
| linear | **10/10** | 10/10 | 10/10 | 10/10 | 10/10 |
| rdna | **10/10** | 10/10 | 10/10 | 10/10 | 10/10 (2/10 with minimap2) |
| heteroplasmy | 36/45 (max err 4) | 36/45 | **45/45** | 28/45 | 36/45 |
| numt | 3/21 (all too LOW) | 3/21 | **21/21** | 10/21 | 7/21 |

Wall times of a full run (plastid, 9 samples, `-t 16`): samtools 7.8 s, SPAdes 20.6 s, SKESA 21.8 s.

## Findings

1. **The templated consensus is error-free where reads map uniquely.** Plastid s01–s07, all linear and rdna
   samples, numt n01: 0 errors/100 kb, full length. The rdna tandem array — which de novo assemblers collapse
   to a third of its length — is perfectly reproduced, since the template keeps the copies apart.
2. **150 bp pairs cannot see over long repeats or deletions; the consensus says N, never a wrong base.**
   - The 10x plastid sample (s09) calls N over most of the inverted repeat: 11,362 N bases, 97% of them in
     the two IR copies (reads mapping equally to both copies get MAPQ 0, and at low depth the consensus has
     nothing left). At the scenario's 40x the IR is called correctly.
   - The 1.5 kb deletion sample (s08) keeps an N run the size of the deletion instead of removing it (no
     pair spans it; BACoN's long reads did). The scorer counts s08's N run as inserted bases and part of
     s09's as mismatches, which is why those two carry all the "errors" of the templated plastid table.
   - Knock-on effect: under the core rule (`--ska-min-freq 1`), one low-depth sample's N can drop another
     sample's SNP sites. Every wrong templated plastid pair is an s08 pair 3 too low, through s09's Ns.
3. **Heteroplasmy: the consensus never reports the minor allele below a 50% tie.** Minor sites go N from a
   5% minor fraction (1/30 sites) through 20–40% (20/30 sites N); distances stay exact up to 40%. At the
   50/50 tie, 4 minor alleles are picked (truth has no right answer there). More sensitive to mixtures than
   BACoN's long-read consensus (clean Q36 bases make a mixture visible), in the same conservative direction.
4. **NUMT/NUPT contamination biases templated distances LOW, not high.** Baiting keeps nearly all
   plastome-like nuclear pairs (NUPT share of filtered bases: 0.4% at 0.5x nuclear depth to 16% at 25x)
   while removing the rest of the nuclear reads (119,788 → 192). The NUPT allele never reaches the
   consensus (0 sites across all samples) — the mixed pileup is called N — but those N also hide the
   sample's own SNPs inside copied segments: distances come out 2–4 where the truth is 3–7. De novo
   assembly separates the NUPT copies into their own contigs instead: **SPAdes + SKA2 scores 21/21 exact**.
   (The extra NUPT contigs are 0.3–2.1% off-reference — below the 10% note threshold, so they are not
   flagged; whether the threshold should drop is a question for the real-data benchmark.)
5. **SPAdes (`--isolate`) assembles without errors on clean samples** and its distances with SKA2 are exact
   on heteroplasmy and numt; with minimap2 they suffer where mixtures fragment the assembly (63–76 contigs
   at 40–50% minor). On the plastid scenario SPAdes + minimap2 is the only 45/45 — the collapsed IR aligned
   to the reference recovers the IR SNPs at both positions.
6. **SKESA confirms the expectation that motivated testing it** (conservative, breaks at repeats): it
   collapses the plastome IR (59,954 of 70,000 bp; 9 contigs at 10x) and the rdna array, loses about 20 bp
   at contig ends (126 bp over the three linear molecules), and 20 NUPT sites go missing in n06. Its
   contigs are error-free on the clean samples but not under contamination or strong mixtures (numt n04 at
   25x nuclear: 257 mismatches, 3,221 duplicated bp; heteroplasmy h050: 6 mismatches). With its default
   SKA2 pairing it beats SPAdes on no scenario (with minimap2 its errors are smaller than SPAdes' on the
   mixture scenarios, but both are behind the SKA2 pairing there). **Recommendation: drop SKESA before
   release unless the real-data benchmark overturns this** (the user decides; `-a skesa` costs little to
   keep while the benchmark runs).
7. **Defaults look right so far**: templated + minimap2 as the overall default (fastest, exact or
   conservative everywhere except numt-heavy genomes), SPAdes + SKA2 as the de novo pairing (and the better
   choice where NUPTs/NUMTs are expected). To be confirmed on real data (Arabidopsis accessions HiSeq,
   human 1KGP chrM) in the benchmark record.

The real-data benchmark (vs GetOrganelle, NOVOPlasty, Snippy, and BACoN on matched ONT runs) is the next
record.
