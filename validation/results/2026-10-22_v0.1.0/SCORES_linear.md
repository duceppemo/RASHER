## Scenario `linear`

### Assemblies

| Run | Sample | Contigs | Circular (truth) | Length / truth | Covered | Duplicated bp | Mismatches | Ins bp | Del bp | Errors /100 kb |
|---|---|---|---|---|---|---|---|---|---|---|
| samtools | l01 | 3 | NA (0) | 30,000 / 30,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| samtools | l02 | 3 | NA (0) | 30,000 / 30,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| samtools | l03 | 3 | NA (0) | 30,000 / 30,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| samtools | l04 | 3 | NA (0) | 30,000 / 30,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| skesa | l01 | 3 | 0 (0) | 29,874 / 30,000 | 0.9958 | 0 | 0 | 0 | 0 | 0.0 |
| skesa | l02 | 3 | 0 (0) | 29,874 / 30,000 | 0.9958 | 0 | 0 | 0 | 0 | 0.0 |
| skesa | l03 | 3 | 0 (0) | 29,874 / 30,000 | 0.9958 | 0 | 0 | 0 | 0 | 0.0 |
| skesa | l04 | 3 | 0 (0) | 29,874 / 30,000 | 0.9958 | 0 | 0 | 0 | 0 | 0.0 |
| spades | l01 | 3 | 0 (0) | 30,000 / 30,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| spades | l02 | 3 | 0 (0) | 30,000 / 30,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| spades | l03 | 3 | 0 (0) | 30,000 / 30,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| spades | l04 | 3 | 0 (0) | 30,000 / 30,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |

### SNP distances

| Run | Method | Samples | Best convention | Exact pairs | Max error | Wrong pairs |
|---|---|---|---|---|---|---|
| samtools | minimap2 | 5/5 | events | 10/10 | 0 |  |
| samtools | ska | 5/5 | events | 10/10 | 0 |  |
| skesa | minimap2 | 5/5 | events | 10/10 | 0 |  |
| skesa | parsnp | 5/5 | events | 10/10 | 0 |  |
| skesa | parsnp_default | 5/5 | events | 10/10 | 0 |  |
| skesa | ska | 5/5 | events | 10/10 | 0 |  |
| skesa | ska_0.5 | 5/5 | events | 10/10 | 0 |  |
| spades | minimap2 | 5/5 | events | 10/10 | 0 |  |
| spades | parsnp | 5/5 | events | 10/10 | 0 |  |
| spades | parsnp_default | 5/5 | events | 10/10 | 0 |  |
| spades | ska | 5/5 | events | 10/10 | 0 |  |
| spades | ska_0.5 | 5/5 | events | 10/10 | 0 |  |

| Run | Wall time of the full run (s) |
|---|---|
| samtools | 2.35 |
| skesa | 6.5 |
| spades | 4.55 |
