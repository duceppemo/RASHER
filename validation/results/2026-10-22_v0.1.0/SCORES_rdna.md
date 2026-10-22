## Scenario `rdna`

### Assemblies

| Run | Sample | Contigs | Circular (truth) | Length / truth | Covered | Duplicated bp | Mismatches | Ins bp | Del bp | Errors /100 kb |
|---|---|---|---|---|---|---|---|---|---|---|
| samtools | t01 | 1 | NA (0) | 22,000 / 22,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| samtools | t02 | 1 | NA (0) | 22,000 / 22,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| samtools | t03 | 1 | NA (0) | 22,000 / 22,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| samtools | t04 | 1 | NA (0) | 22,000 / 22,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| skesa | t01 | 1 | 0 (0) | 6,958 / 22,000 | 0.4527 | 0 | 0 | 0 | 0 | 0.0 |
| skesa | t02 | 1 | 0 (0) | 6,958 / 22,000 | 0.4527 | 0 | 0 | 0 | 0 | 0.0 |
| skesa | t03 | 1 | 0 (0) | 6,958 / 22,000 | 0.4527 | 0 | 0 | 0 | 0 | 0.0 |
| skesa | t04 | 1 | 0 (0) | 6,958 / 22,000 | 0.4527 | 0 | 0 | 0 | 0 | 0.0 |
| spades | t01 | 2 | 0 (0) | 7,152 / 22,000 | 0.3251 | 0 | 0 | 0 | 0 | 0.0 |
| spades | t02 | 2 | 0 (0) | 7,152 / 22,000 | 0.3251 | 0 | 0 | 0 | 0 | 0.0 |
| spades | t03 | 2 | 0 (0) | 7,152 / 22,000 | 0.3251 | 0 | 0 | 0 | 0 | 0.0 |
| spades | t04 | 2 | 0 (0) | 7,152 / 22,000 | 0.3251 | 0 | 0 | 0 | 0 | 0.0 |

### SNP distances

| Run | Method | Samples | Best convention | Exact pairs | Max error | Wrong pairs |
|---|---|---|---|---|---|---|
| samtools | minimap2 | 5/5 | sites | 10/10 | 0 |  |
| samtools | ska | 5/5 | events | 10/10 | 0 |  |
| skesa | minimap2 | 5/5 | events | 2/10 | 2 | Reference-t02: 5 (truth 4); Reference-t03: 7 (truth 6); Reference-t04: 7 (truth 6); t01-t02: 5 (truth 4); t01-t03: 7 (truth 6); t01-t04: 7 (truth 6) … |
| skesa | parsnp | 5/5 | no_repeat | 10/10 | 0 |  |
| skesa | ska | 5/5 | events | 10/10 | 0 |  |
| skesa | ska_0.5 | 5/5 | events | 10/10 | 0 |  |
| spades | minimap2 | 5/5 | events | 10/10 | 0 |  |
| spades | parsnp | 5/5 | no_repeat | 10/10 | 0 |  |
| spades | ska | 5/5 | events | 10/10 | 0 |  |
| spades | ska_0.5 | 5/5 | events | 10/10 | 0 |  |

| Run | Wall time of the full run (s) |
|---|---|
| samtools | 2.49 |
| skesa | 6.59 |
| spades | 5.4 |
