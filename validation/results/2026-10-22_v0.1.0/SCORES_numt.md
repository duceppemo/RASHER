## Scenario `numt`

### Assemblies

| Run | Sample | Contigs | Circular (truth) | Length / truth | Covered | Duplicated bp | Mismatches | Ins bp | Del bp | Errors /100 kb |
|---|---|---|---|---|---|---|---|---|---|---|
| samtools | n01 | 1 | NA (1) | 120,000 / 120,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| samtools | n02 | 1 | NA (1) | 120,000 / 120,000 | 1.0000 | 0 | 8 | 0 | 0 | 6.7 |
| samtools | n03 | 1 | NA (1) | 120,000 / 120,000 | 1.0000 | 0 | 507 | 0 | 0 | 422.5 |
| samtools | n04 | 1 | NA (1) | 120,000 / 120,000 | 1.0000 | 0 | 973 | 0 | 0 | 810.8 |
| samtools | n05 | 1 | NA (1) | 120,000 / 120,000 | 1.0000 | 0 | 40 | 0 | 0 | 33.3 |
| samtools | n06 | 1 | NA (1) | 120,000 / 120,000 | 1.0000 | 0 | 49 | 0 | 0 | 40.8 |
| skesa | n01 | 1 | 0 (1) | 107,998 / 120,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| skesa | n02 | 2 | 0 (1) | 107,976 / 120,000 | 0.9998 | 0 | 0 | 0 | 0 | 0.0 |
| skesa | n03 | 8 | 0 (1) | 111,229 / 120,000 | 0.9996 | 0 | 58 | 3 | 5 | 59.8 |
| skesa | n04 | 33 | 0 (1) | 136,708 / 120,000 | 0.9815 | 3,221 | 257 | 47 | 38 | 282.1 |
| skesa | n05 | 2 | 0 (1) | 107,976 / 120,000 | 0.9998 | 0 | 0 | 0 | 0 | 0.0 |
| skesa | n06 | 22 | 0 (1) | 106,568 / 120,000 | 0.9867 | 0 | 0 | 0 | 0 | 0.0 |
| spades | n01 | 8 | 0 (1) | 109,560 / 120,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| spades | n02 | 20 | 0 (1) | 114,728 / 120,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| spades | n03 | 13 | 0 (1) | 123,655 / 120,000 | 1.0000 | 0 | 58 | 5 | 5 | 61.9 |
| spades | n04 | 83 | 0 (1) | 137,161 / 120,000 | 1.0000 | 0 | 135 | 20 | 17 | 149.0 |
| spades | n05 | 6 | 0 (1) | 108,768 / 120,000 | 1.0000 | 0 | 0 | 0 | 0 | 0.0 |
| spades | n06 | 29 | 0 (1) | 111,342 / 120,000 | 1.0000 | 0 | 2 | 0 | 0 | 1.8 |

### SNP distances

| Run | Method | Samples | Best convention | Exact pairs | Max error | Wrong pairs |
|---|---|---|---|---|---|---|
| samtools | minimap2 | 7/7 | events | 3/21 | 5 | Reference-n02: 2 (truth 3); Reference-n03: 3 (truth 5); Reference-n04: 2 (truth 5); Reference-n05: 2 (truth 3); Reference-n06: 4 (truth 7); n01-n02: 2 (truth 3) … |
| samtools | ska | 7/7 | events | 3/21 | 5 | Reference-n02: 2 (truth 3); Reference-n03: 3 (truth 5); Reference-n04: 2 (truth 5); Reference-n05: 2 (truth 3); Reference-n06: 4 (truth 7); n01-n02: 2 (truth 3) … |
| skesa | minimap2 | 7/7 | events | 7/21 | 4 | Reference-n03: 4 (truth 5); Reference-n04: 2 (truth 5); Reference-n06: 4 (truth 7); n01-n03: 4 (truth 5); n01-n04: 2 (truth 5); n01-n06: 4 (truth 7) … |
| skesa | ska | 7/7 | events | 7/21 | 2 | Reference-n03: 4 (truth 5); Reference-n04: 4 (truth 5); Reference-n06: 6 (truth 7); n01-n03: 4 (truth 5); n01-n04: 4 (truth 5); n01-n06: 6 (truth 7) … |
| spades | minimap2 | 7/7 | events | 10/21 | 13 | Reference-n04: 13 (truth 5); Reference-n06: 6 (truth 7); n01-n04: 13 (truth 5); n01-n06: 6 (truth 7); n02-n04: 16 (truth 8); n02-n06: 9 (truth 10) … |
| spades | ska | 7/7 | events | 21/21 | 0 |  |

### Reads by origin

Reads (both mates of a pair counted) simulated / baited (1_extracted) / kept by the filter (2_filtered). NUPT reads: nuclear pairs whose fragment covers at least 100 bp of an insert; nuclear: the other nuclear reads. Share: of the filtered reads (bases).

| Run | Sample | Nuclear depth (genome) | Organelle reads | NUPT reads | Nuclear reads | NUPT share of filtered reads (bases) |
|---|---|---|---|---|---|---|
| samtools | n01 | 0.5x (A) | 40000/40000/40000 | 176/154/154 | 2376/2/2 | 0.4% (0.4%) |
| samtools | n02 | 2x (A) | 40000/40000/40000 | 628/552/552 | 9574/6/6 | 1.4% (1.4%) |
| samtools | n03 | 10x (A) | 40000/40000/40000 | 3124/2680/2680 | 47880/56/56 | 6.3% (6.3%) |
| samtools | n04 | 25x (A) | 40000/40000/40000 | 7752/6724/6724 | 119756/118/118 | 14.4% (14.4%) |
| samtools | n05 | 10x (B) | 40000/40000/40000 | 2940/2940/2940 | 48060/74/74 | 6.8% (6.8%) |
| samtools | n06 | 25x (B) | 40000/40000/40000 | 7712/7712/7712 | 119788/192/192 | 16.1% (16.1%) |
| skesa | n01 | 0.5x (A) | 40000/40000/40000 | 176/154/154 | 2376/2/2 | 0.4% (0.4%) |
| skesa | n02 | 2x (A) | 40000/40000/40000 | 628/552/552 | 9574/6/6 | 1.4% (1.4%) |
| skesa | n03 | 10x (A) | 40000/40000/40000 | 3124/2680/2680 | 47880/56/56 | 6.3% (6.3%) |
| skesa | n04 | 25x (A) | 40000/40000/40000 | 7752/6724/6724 | 119756/118/118 | 14.4% (14.4%) |
| skesa | n05 | 10x (B) | 40000/40000/40000 | 2940/2940/2940 | 48060/74/74 | 6.8% (6.8%) |
| skesa | n06 | 25x (B) | 40000/40000/40000 | 7712/7712/7712 | 119788/192/192 | 16.1% (16.1%) |
| spades | n01 | 0.5x (A) | 40000/40000/40000 | 176/154/154 | 2376/2/2 | 0.4% (0.4%) |
| spades | n02 | 2x (A) | 40000/40000/40000 | 628/552/552 | 9574/6/6 | 1.4% (1.4%) |
| spades | n03 | 10x (A) | 40000/40000/40000 | 3124/2680/2680 | 47880/56/56 | 6.3% (6.3%) |
| spades | n04 | 25x (A) | 40000/40000/40000 | 7752/6724/6724 | 119756/118/118 | 14.4% (14.4%) |
| spades | n05 | 10x (B) | 40000/40000/40000 | 2940/2940/2940 | 48060/74/74 | 6.8% (6.8%) |
| spades | n06 | 25x (B) | 40000/40000/40000 | 7712/7712/7712 | 119788/192/192 | 16.1% (16.1%) |

### NUPT reads baited, by insert

Baited / simulated reads of each insert (length, identity to the plastome in genome A; in genome B every insert is 99.9% identical).

| Run | Sample | #1 0.5 kb 99.9% | #2 0.5 kb 90% | #3 2 kb 99% | #4 2 kb 97% | #5 5 kb 97% | #6 5 kb 99.9% | #7 15 kb 90% | #8 15 kb 99% |
|---|---|---|---|---|---|---|---|---|---|
| samtools | n01 | 2/2 | 0/0 | 8/8 | 14/14 | 22/22 | 24/24 | 38/60 | 46/46 |
| samtools | n02 | 14/14 | 0/4 | 30/30 | 30/30 | 86/86 | 56/56 | 128/200 | 208/208 |
| samtools | n03 | 48/48 | 0/34 | 128/128 | 144/144 | 344/344 | 354/354 | 622/1032 | 1040/1040 |
| samtools | n04 | 102/102 | 0/100 | 340/340 | 386/388 | 876/876 | 888/888 | 1606/2532 | 2526/2526 |
| samtools | n05 | 36/36 | 26/26 | 154/154 | 138/138 | 298/298 | 338/338 | 1008/1008 | 942/942 |
| samtools | n06 | 118/118 | 108/108 | 358/358 | 344/344 | 838/838 | 808/808 | 2558/2558 | 2580/2580 |
| skesa | n01 | 2/2 | 0/0 | 8/8 | 14/14 | 22/22 | 24/24 | 38/60 | 46/46 |
| skesa | n02 | 14/14 | 0/4 | 30/30 | 30/30 | 86/86 | 56/56 | 128/200 | 208/208 |
| skesa | n03 | 48/48 | 0/34 | 128/128 | 144/144 | 344/344 | 354/354 | 622/1032 | 1040/1040 |
| skesa | n04 | 102/102 | 0/100 | 340/340 | 386/388 | 876/876 | 888/888 | 1606/2532 | 2526/2526 |
| skesa | n05 | 36/36 | 26/26 | 154/154 | 138/138 | 298/298 | 338/338 | 1008/1008 | 942/942 |
| skesa | n06 | 118/118 | 108/108 | 358/358 | 344/344 | 838/838 | 808/808 | 2558/2558 | 2580/2580 |
| spades | n01 | 2/2 | 0/0 | 8/8 | 14/14 | 22/22 | 24/24 | 38/60 | 46/46 |
| spades | n02 | 14/14 | 0/4 | 30/30 | 30/30 | 86/86 | 56/56 | 128/200 | 208/208 |
| spades | n03 | 48/48 | 0/34 | 128/128 | 144/144 | 344/344 | 354/354 | 622/1032 | 1040/1040 |
| spades | n04 | 102/102 | 0/100 | 340/340 | 386/388 | 876/876 | 888/888 | 1606/2532 | 2526/2526 |
| spades | n05 | 36/36 | 26/26 | 154/154 | 138/138 | 298/298 | 338/338 | 1008/1008 | 942/942 |
| spades | n06 | 118/118 | 108/108 | 358/358 | 344/344 | 838/838 | 808/808 | 2558/2558 | 2580/2580 |

### Consensus

Consensus errors (mismatched, inserted and deleted bp) inside the plastome segments copied into the nucleus and elsewhere, per aligned length; at the substitutions that distinguish the NUPTs from the plastome, what the assembly has (plastome base / NUPT base / N / other / split: different bases in different contigs / missing) and the number of sites where the NUPT base is in any contig; the sample's own SNPs inside the copied segments (where the NUPT reads carry the reference base): called / reference base / N / split or other; distance: SNPs to the reference reported by each method (truth).

| Run | Sample | Errors in copied segments (aligned bp) | Errors elsewhere (aligned bp) | NUPT sites: plastome / NUPT / N / other / split / missing | NUPT base in a contig | Own SNPs in copied segments | Distance to reference (truth) |
|---|---|---|---|---|---|---|---|
| samtools | n01 | 0 (45,000) | 0 (75,000) | 1666 / 0 / 0 / 0 / 0 / 0 | 0 | 0 / 0 / 0 / 0 | minimap2 0; ska 0 (0) |
| samtools | n02 | 8 (45,000) | 0 (75,000) | 1658 / 0 / 8 / 0 / 0 / 0 | 0 | 1 / 0 / 0 / 0 | minimap2 2; ska 2 (3) |
| samtools | n03 | 507 (45,000) | 0 (74,998) | 1163 / 0 / 503 / 0 / 0 / 0 | 0 | 0 / 0 / 2 / 0 | minimap2 3; ska 3 (5) |
| samtools | n04 | 973 (45,000) | 0 (75,000) | 704 / 0 / 962 / 0 / 0 / 0 | 0 | 0 / 0 / 3 / 0 | minimap2 2; ska 2 (5) |
| samtools | n05 | 40 (45,000) | 0 (75,000) | 7 / 0 / 39 / 0 / 0 / 0 | 0 | 0 / 0 / 1 / 0 | minimap2 2; ska 2 (3) |
| samtools | n06 | 49 (45,000) | 0 (75,000) | 0 / 0 / 46 / 0 / 0 / 0 | 0 | 0 / 0 / 3 / 0 | minimap2 4; ska 4 (7) |
| skesa | n01 | 0 (45,000) | 0 (62,998) | 1666 / 0 / 0 / 0 / 0 / 0 | 0 | 0 / 0 / 0 / 0 | minimap2 0; ska 0 (0) |
| skesa | n02 | 0 (45,000) | 0 (62,976) | 1666 / 0 / 0 / 0 / 0 / 0 | 0 | 1 / 0 / 0 / 0 | minimap2 3; ska 3 (3) |
| skesa | n03 | 66 (47,402) | 0 (62,976) | 1608 / 0 / 0 / 0 / 57 / 1 | 57 | 1 / 0 / 0 / 1 | minimap2 4; ska 4 (5) |
| skesa | n04 | 342 (58,246) | 0 (62,978) | 1387 / 0 / 0 / 0 / 254 / 25 | 253 | 1 / 0 / 0 / 2 | minimap2 2; ska 4 (5) |
| skesa | n05 | 0 (45,000) | 0 (62,976) | 46 / 0 / 0 / 0 / 0 / 0 | 0 | 1 / 0 / 0 / 0 | minimap2 3; ska 3 (3) |
| skesa | n06 | 0 (43,618) | 0 (62,950) | 26 / 0 / 0 / 0 / 0 / 20 | 0 | 2 / 0 / 0 / 1 | minimap2 4; ska 6 (7) |
| spades | n01 | 0 (45,000) | 0 (63,000) | 1666 / 0 / 0 / 0 / 0 / 0 | 0 | 0 / 0 / 0 / 0 | minimap2 0; ska 0 (0) |
| spades | n02 | 0 (45,000) | 0 (63,000) | 1666 / 0 / 0 / 0 / 0 / 0 | 0 | 1 / 0 / 0 / 0 | minimap2 3; ska 3 (3) |
| spades | n03 | 68 (46,783) | 0 (63,154) | 1608 / 0 / 0 / 0 / 58 / 0 | 58 | 2 / 0 / 0 / 0 | minimap2 5; ska 5 (5) |
| spades | n04 | 172 (52,200) | 0 (63,231) | 1530 / 11 / 0 / 0 / 122 / 3 | 132 | 1 / 0 / 0 / 2 | minimap2 13; ska 5 (5) |
| spades | n05 | 0 (45,000) | 0 (63,000) | 46 / 0 / 0 / 0 / 0 / 0 | 0 | 1 / 0 / 0 / 0 | minimap2 3; ska 3 (3) |
| spades | n06 | 2 (45,613) | 0 (63,231) | 41 / 2 / 0 / 0 / 0 / 3 | 2 | 3 / 0 / 0 / 0 | minimap2 6; ska 7 (7) |

| Run | Wall time of the full run (s) |
|---|---|
| samtools | 12.61 |
| skesa | 30.33 |
| spades | 26.13 |
