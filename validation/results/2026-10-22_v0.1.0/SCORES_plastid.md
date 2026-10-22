## Scenario `plastid`

### Assemblies

| Run | Sample | Contigs | Circular (truth) | Length / truth | Covered | Duplicated bp | Mismatches | Ins bp | Del bp | Errors /100 kb |
|---|---|---|---|---|---|---|---|---|---|---|
| samtools | s01 | 1 | NA (1) | 70,000 / 70,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| samtools | s02 | 1 | NA (1) | 70,000 / 70,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| samtools | s03 | 1 | NA (1) | 70,000 / 70,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| samtools | s04 | 1 | NA (1) | 70,000 / 70,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| samtools | s05 | 1 | NA (1) | 70,000 / 70,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| samtools | s06 | 1 | NA (1) | 70,000 / 70,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| samtools | s07 | 1 | NA (1) | 70,005 / 70,005 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| samtools | s08 | 1 | NA (1) | 70,000 / 68,500 | 1.0256 | 0 | 0 | 1500 | 0 | 2142.9 |
| samtools | s09 | 1 | NA (1) | 70,000 / 70,000 | 0.9495 | 0 | 2053 | 0 | 0 | 3479.5 |
| skesa | s01 | 3 | 0 (1) | 59,954 / 70,000 | 0.9992 | 0 | 0 | 0 | 0 | 0.0 |
| skesa | s02 | 1 | 0 (1) | 59,998 / 70,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| skesa | s03 | 2 | 0 (1) | 59,976 / 70,000 | 0.9996 | 0 | 0 | 0 | 0 | 0.0 |
| skesa | s04 | 3 | 0 (1) | 59,954 / 70,000 | 0.9992 | 0 | 0 | 0 | 0 | 0.0 |
| skesa | s05 | 3 | 0 (1) | 59,954 / 70,000 | 0.9992 | 0 | 0 | 0 | 0 | 0.0 |
| skesa | s06 | 3 | 0 (1) | 59,954 / 70,000 | 0.9992 | 0 | 0 | 0 | 0 | 0.0 |
| skesa | s07 | 3 | 0 (1) | 59,959 / 70,005 | 0.9992 | 0 | 0 | 0 | 0 | 0.0 |
| skesa | s08 | 2 | 0 (1) | 58,476 / 68,500 | 0.9996 | 0 | 0 | 0 | 0 | 0.0 |
| skesa | s09 | 9 | 0 (1) | 59,143 / 70,000 | 0.9857 | 0 | 0 | 0 | 0 | 0.0 |
| spades | s01 | 3 | 0 (1) | 60,000 / 70,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| spades | s02 | 3 | 0 (1) | 60,000 / 70,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| spades | s03 | 3 | 0 (1) | 60,000 / 70,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| spades | s04 | 3 | 0 (1) | 60,000 / 70,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| spades | s05 | 3 | 0 (1) | 60,000 / 70,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| spades | s06 | 3 | 0 (1) | 60,000 / 70,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| spades | s07 | 3 | 0 (1) | 60,005 / 70,005 | 0.9999 | 0 | 0 | 0 | 0 | 0.0 |
| spades | s08 | 3 | 0 (1) | 58,500 / 68,500 | 0.9744 | 0 | 0 | 0 | 0 | 0.0 |
| spades | s09 | 5 | 0 (1) | 60,184 / 70,000 | 0.9979 | 0 | 0 | 0 | 0 | 0.0 |

### SNP distances

| Run | Method | Samples | Best convention | Exact pairs | Max error | Wrong pairs |
|---|---|---|---|---|---|---|
| samtools | minimap2 | 10/10 | no_repeat | 36/45 | 3 | Reference-s08: 37 (truth 40); s01-s08: 37 (truth 40); s02-s08: 37 (truth 40); s03-s08: 38 (truth 41); s04-s08: 39 (truth 42); s05-s08: 52 (truth 55) … |
| samtools | ska | 10/10 | no_repeat | 28/45 | 2 | Reference-s07: 21 (truth 20); Reference-s08: 38 (truth 40); s01-s07: 21 (truth 20); s01-s08: 38 (truth 40); s02-s07: 21 (truth 20); s02-s08: 38 (truth 40) … |
| skesa | minimap2 | 10/10 | no_repeat | 36/45 | 3 | Reference-s08: 37 (truth 40); s01-s08: 37 (truth 40); s02-s08: 37 (truth 40); s03-s08: 38 (truth 41); s04-s08: 39 (truth 42); s05-s08: 52 (truth 55) … |
| skesa | parsnp | 10/10 | events | 36/45 | 3 | Reference-s08: 37 (truth 40); s01-s08: 37 (truth 40); s02-s08: 37 (truth 40); s03-s08: 38 (truth 41); s04-s08: 39 (truth 42); s05-s08: 52 (truth 55) … |
| skesa | parsnp_default | 10/10 | events | 36/45 | 3 | Reference-s08: 37 (truth 40); s01-s08: 37 (truth 40); s02-s08: 37 (truth 40); s03-s08: 38 (truth 41); s04-s08: 39 (truth 42); s05-s08: 52 (truth 55) … |
| skesa | ska | 10/10 | events | 36/45 | 3 | Reference-s08: 37 (truth 40); s01-s08: 37 (truth 40); s02-s08: 37 (truth 40); s03-s08: 38 (truth 41); s04-s08: 39 (truth 42); s05-s08: 52 (truth 55) … |
| skesa | ska_0.5 | 10/10 | events | 44/45 | 3 | s08-s09: 58 (truth 61) |
| spades | minimap2 | 10/10 | sites | 45/45 | 0 |  |
| spades | parsnp | 10/10 | events | 36/45 | 2 | Reference-s08: 38 (truth 40); s01-s08: 38 (truth 40); s02-s08: 38 (truth 40); s03-s08: 39 (truth 41); s04-s08: 40 (truth 42); s05-s08: 53 (truth 55) … |
| spades | parsnp_default | 10/10 | events | 36/45 | 2 | Reference-s08: 38 (truth 40); s01-s08: 38 (truth 40); s02-s08: 38 (truth 40); s03-s08: 39 (truth 41); s04-s08: 40 (truth 42); s05-s08: 53 (truth 55) … |
| spades | ska | 10/10 | events | 36/45 | 1 | Reference-s08: 39 (truth 40); s01-s08: 39 (truth 40); s02-s08: 39 (truth 40); s03-s08: 40 (truth 41); s04-s08: 41 (truth 42); s05-s08: 54 (truth 55) … |
| spades | ska_0.5 | 10/10 | events | 44/45 | 1 | s08-s09: 60 (truth 61) |

| Run | Wall time of the full run (s) |
|---|---|
| samtools | 7.84 |
| skesa | 21.8 |
| spades | 20.55 |
