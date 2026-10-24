"""snippy core.aln (SNP-only fasta) -> snp_distances.tsv in RASHER's matrix format."""
import sys
from pathlib import Path

aln, out = sys.argv[1], sys.argv[2]
seqs, name = {}, None
for line in Path(aln).read_text().splitlines():
    if line.startswith(">"):
        name = line[1:].strip()
        seqs[name] = []
    else:
        seqs[name].append(line.strip())
seqs = {k: "".join(v) for k, v in seqs.items()}
names = sorted(seqs)
with open(out, "w") as fh:
    fh.write("snp-dists\t" + "\t".join(names) + "\n")
    for a in names:
        row = [str(sum(1 for x, y in zip(seqs[a], seqs[b]) if x != y)) for b in names]
        fh.write(a + "\t" + "\t".join(row) + "\n")
print(out, len(names), "samples,", len(next(iter(seqs.values()))), "SNP columns")
