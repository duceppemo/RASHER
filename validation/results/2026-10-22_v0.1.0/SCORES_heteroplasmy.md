## Scenario `heteroplasmy`

### Assemblies

| Run | Sample | Contigs | Circular (truth) | Length / truth | Covered | Duplicated bp | Mismatches | Ins bp | Del bp | Errors /100 kb |
|---|---|---|---|---|---|---|---|---|---|---|
| samtools | h000 | 1 | NA (1) | 70,000 / 70,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| samtools | h002 | 1 | NA (1) | 70,000 / 70,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| samtools | h005 | 1 | NA (1) | 70,000 / 70,000 | 1.0000 | 0 | 1 | 0 | 0 | 1.4 |
| samtools | h010 | 1 | NA (1) | 70,000 / 70,000 | 1.0000 | 0 | 8 | 0 | 0 | 11.4 |
| samtools | h020 | 1 | NA (1) | 70,000 / 70,000 | 1.0000 | 0 | 20 | 0 | 0 | 28.6 |
| samtools | h030 | 1 | NA (1) | 70,000 / 70,000 | 1.0000 | 0 | 20 | 0 | 0 | 28.6 |
| samtools | h040 | 1 | NA (1) | 70,000 / 70,000 | 1.0000 | 0 | 22 | 0 | 0 | 31.4 |
| samtools | h050 | 1 | NA (1) | 70,002 / 70,000 | 1.0000 | 0 | 29 | 3 | 1 | 47.1 |
| samtools | h020s | 1 | NA (1) | 70,000 / 70,000 | 1.0000 | 0 | 20 | 0 | 0 | 28.6 |
| skesa | h000 | 3 | 0 (1) | 59,954 / 70,000 | 0.9992 | 0 | 0 | 0 | 0 | 0.0 |
| skesa | h002 | 3 | 0 (1) | 59,954 / 70,000 | 0.9992 | 0 | 0 | 0 | 0 | 0.0 |
| skesa | h005 | 3 | 0 (1) | 59,954 / 70,000 | 0.9992 | 0 | 0 | 0 | 0 | 0.0 |
| skesa | h010 | 3 | 0 (1) | 59,954 / 70,000 | 0.9992 | 0 | 0 | 0 | 0 | 0.0 |
| skesa | h020 | 3 | 0 (1) | 59,954 / 70,000 | 0.9992 | 0 | 0 | 0 | 0 | 0.0 |
| skesa | h030 | 11 | 0 (1) | 59,621 / 70,000 | 0.9937 | 0 | 0 | 0 | 0 | 0.0 |
| skesa | h040 | 21 | 0 (1) | 59,336 / 70,000 | 0.9889 | 0 | 1 | 0 | 0 | 1.7 |
| skesa | h050 | 19 | 0 (1) | 59,071 / 70,000 | 0.9846 | 0 | 6 | 0 | 5 | 18.6 |
| skesa | h020s | 4 | 0 (1) | 59,911 / 70,000 | 0.9985 | 0 | 0 | 0 | 0 | 0.0 |
| spades | h000 | 3 | 0 (1) | 60,000 / 70,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| spades | h002 | 3 | 0 (1) | 60,000 / 70,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| spades | h005 | 3 | 0 (1) | 60,000 / 70,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| spades | h010 | 3 | 0 (1) | 60,000 / 70,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| spades | h020 | 10 | 0 (1) | 60,685 / 70,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| spades | h030 | 6 | 0 (1) | 60,309 / 70,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| spades | h040 | 63 | 0 (1) | 66,643 / 70,000 | 0.9997 | 0 | 3 | 0 | 0 | 4.9 |
| spades | h050 | 76 | 0 (1) | 68,034 / 70,000 | 0.9997 | 0 | 1 | 0 | 5 | 9.9 |
| spades | h020s | 3 | 0 (1) | 60,000 / 70,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |

### SNP distances

| Run | Method | Samples | Best convention | Exact pairs | Max error | Wrong pairs |
|---|---|---|---|---|---|---|
| samtools | minimap2 | 10/10 | events | 36/45 | 4 | Reference-h050: 4 (truth 0); h000-h050: 4 (truth 0); h002-h050: 4 (truth 0); h005-h050: 4 (truth 0); h010-h050: 4 (truth 0); h020-h050: 4 (truth 0) … |
| samtools | ska | 10/10 | events | 36/45 | 2 | Reference-h050: 2 (truth 0); h000-h050: 2 (truth 0); h002-h050: 2 (truth 0); h005-h050: 2 (truth 0); h010-h050: 2 (truth 0); h020-h050: 2 (truth 0) … |
| skesa | minimap2 | 10/10 | events | 36/45 | 1 | Reference-h050: 1 (truth 0); h000-h050: 1 (truth 0); h002-h050: 1 (truth 0); h005-h050: 1 (truth 0); h010-h050: 1 (truth 0); h020-h050: 1 (truth 0) … |
| skesa | ska | 10/10 | events | 36/45 | 1 | Reference-h050: 1 (truth 0); h000-h050: 1 (truth 0); h002-h050: 1 (truth 0); h005-h050: 1 (truth 0); h010-h050: 1 (truth 0); h020-h050: 1 (truth 0) … |
| spades | minimap2 | 10/10 | events | 28/45 | 2 | Reference-h040: 2 (truth 0); Reference-h050: 1 (truth 0); h000-h040: 2 (truth 0); h000-h050: 1 (truth 0); h002-h040: 2 (truth 0); h002-h050: 1 (truth 0) … |
| spades | ska | 10/10 | events | 45/45 | 0 |  |

### Minor-haplotype sites

At the 30 SNP sites where the minor haplotype differs (inverted-repeat SNPs count twice), what the assembly has: the major allele (correct), the minor allele, N, another base, both alleles in different contigs or repeat copies (split) or nothing (missing); then the 3 indels (mdel1, mins1, mdel2): M major, m minor, N, o other, s split, - missing. Distance: SNPs to the reference reported by each method (truth).

| Run | Sample | Minor fraction | Major | Minor | N | Other | Split | Missing | Indels | N bases | Distance to reference (truth) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| samtools | h000 | 0.00 | 30 (100%) | 0 (0%) | 0 (0%) | 0 | 0 | 0 | M/M/M | 0 | minimap2 0; ska 0 (0) |
| samtools | h002 | 0.02 | 30 (100%) | 0 (0%) | 0 (0%) | 0 | 0 | 0 | M/M/M | 0 | minimap2 0; ska 0 (0) |
| samtools | h005 | 0.05 | 29 (97%) | 0 (0%) | 1 (3%) | 0 | 0 | 0 | M/M/M | 1 | minimap2 0; ska 0 (0) |
| samtools | h010 | 0.10 | 22 (73%) | 0 (0%) | 8 (27%) | 0 | 0 | 0 | M/M/M | 8 | minimap2 0; ska 0 (0) |
| samtools | h020 | 0.20 | 10 (33%) | 0 (0%) | 20 (67%) | 0 | 0 | 0 | M/M/M | 20 | minimap2 0; ska 0 (0) |
| samtools | h020s | 0.20 | 10 (33%) | 0 (0%) | 20 (67%) | 0 | 0 | 0 | M/M/M | 20 | minimap2 12; ska 12 (12) |
| samtools | h030 | 0.30 | 10 (33%) | 0 (0%) | 20 (67%) | 0 | 0 | 0 | M/M/M | 20 | minimap2 0; ska 0 (0) |
| samtools | h040 | 0.40 | 6 (20%) | 0 (0%) | 20 (67%) | 0 | 4 | 0 | M/M/M | 22 | minimap2 0; ska 0 (0) |
| samtools | h050 | 0.50 | 0 (0%) | 4 (13%) | 24 (80%) | 0 | 2 | 0 | M/m/o | 25 | minimap2 4; ska 2 (0) |
| skesa | h000 | 0.00 | 30 (100%) | 0 (0%) | 0 (0%) | 0 | 0 | 0 | M/M/M | 0 | minimap2 0; ska 0 (0) |
| skesa | h002 | 0.02 | 30 (100%) | 0 (0%) | 0 (0%) | 0 | 0 | 0 | M/M/M | 0 | minimap2 0; ska 0 (0) |
| skesa | h005 | 0.05 | 30 (100%) | 0 (0%) | 0 (0%) | 0 | 0 | 0 | M/M/M | 0 | minimap2 0; ska 0 (0) |
| skesa | h010 | 0.10 | 30 (100%) | 0 (0%) | 0 (0%) | 0 | 0 | 0 | M/M/M | 0 | minimap2 0; ska 0 (0) |
| skesa | h020 | 0.20 | 30 (100%) | 0 (0%) | 0 (0%) | 0 | 0 | 0 | M/M/M | 0 | minimap2 0; ska 0 (0) |
| skesa | h020s | 0.20 | 28 (93%) | 0 (0%) | 0 (0%) | 0 | 0 | 2 | M/M/M | 0 | minimap2 12; ska 12 (12) |
| skesa | h030 | 0.30 | 17 (57%) | 0 (0%) | 0 (0%) | 0 | 0 | 13 | M/M/M | 0 | minimap2 0; ska 0 (0) |
| skesa | h040 | 0.40 | 8 (27%) | 1 (3%) | 0 (0%) | 0 | 0 | 21 | -/M/- | 0 | minimap2 0; ska 0 (0) |
| skesa | h050 | 0.50 | 2 (7%) | 6 (20%) | 0 (0%) | 0 | 0 | 22 | M/-/m | 0 | minimap2 1; ska 1 (0) |
| spades | h000 | 0.00 | 30 (100%) | 0 (0%) | 0 (0%) | 0 | 0 | 0 | M/M/M | 0 | minimap2 0; ska 0 (0) |
| spades | h002 | 0.02 | 30 (100%) | 0 (0%) | 0 (0%) | 0 | 0 | 0 | M/M/M | 0 | minimap2 0; ska 0 (0) |
| spades | h005 | 0.05 | 30 (100%) | 0 (0%) | 0 (0%) | 0 | 0 | 0 | M/M/M | 0 | minimap2 0; ska 0 (0) |
| spades | h010 | 0.10 | 30 (100%) | 0 (0%) | 0 (0%) | 0 | 0 | 0 | M/M/M | 0 | minimap2 0; ska 0 (0) |
| spades | h020 | 0.20 | 28 (93%) | 0 (0%) | 0 (0%) | 0 | 0 | 2 | M/M/M | 0 | minimap2 0; ska 0 (0) |
| spades | h020s | 0.20 | 30 (100%) | 0 (0%) | 0 (0%) | 0 | 0 | 0 | M/M/M | 0 | minimap2 12; ska 12 (12) |
| spades | h030 | 0.30 | 28 (93%) | 0 (0%) | 0 (0%) | 0 | 0 | 2 | M/M/M | 0 | minimap2 0; ska 0 (0) |
| spades | h040 | 0.40 | 7 (23%) | 3 (10%) | 0 (0%) | 0 | 0 | 20 | -/M/M | 0 | minimap2 2; ska 0 (0) |
| spades | h050 | 0.50 | 4 (13%) | 1 (3%) | 0 (0%) | 0 | 0 | 25 | M/M/m | 0 | minimap2 1; ska 0 (0) |

| Run | Wall time of the full run (s) |
|---|---|
| samtools | 14.96 |
| skesa | 44.74 |
| spades | 29.35 |
