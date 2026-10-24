"""Score competitor assemblies: pick each sample's assembly fasta, map to the reference, report
reference coverage %, contig count, assembly length. Usage: comp_score.py TOOL DATASET REF"""
import statistics
import subprocess
import sys
from pathlib import Path

tool, ds, ref = sys.argv[1:4]
base = Path("/media/8tb/analyses/bacon_benchmark/rasher") / tool / ds

def pick(d: Path):
    if tool == "getorganelle":
        pats = ["*complete.graph1.1.path_sequence.fasta", "*graph1.1.path_sequence.fasta"]
    else:
        pats = ["Circularized_assembly_1_*.fasta", "Option_1_*.fasta", "Contigs_1_*.fasta"]
    for p in pats:
        hits = sorted(d.glob(p))
        if hits:
            return hits[0]
    return None

rows = []
for d in sorted(p for p in base.iterdir() if p.is_dir()):
    fa = pick(d)
    if fa is None:
        rows.append((d.name, 0.0, 0, 0, "NO_ASSEMBLY"))
        continue
    text = fa.read_text()
    ncont = text.count(">")
    length = sum(len(line) for line in text.splitlines() if not line.startswith(">"))
    paf = subprocess.run(["minimap2", "-x", "asm5", "--secondary=no", ref, str(fa)],
                         capture_output=True, text=True).stdout
    cov, reflen = set(), 1
    for line in paf.splitlines():
        f = line.split("\t")
        reflen = int(f[6])
        cov.update(range(int(f[7]), int(f[8])))
    rows.append((d.name, 100 * len(cov) / reflen, ncont, length, fa.name.split("_")[0]))

for r in rows:
    print(f"{r[0]}\t{r[1]:.2f}\t{r[2]}\t{r[3]}\t{r[4]}")
covs = [r[1] for r in rows]
conts = [r[2] for r in rows]
lens = [r[3] for r in rows]
full = sum(1 for c in covs if c >= 99.9)
print(f"# {tool}/{ds}: n={len(rows)}, ref-cov median {statistics.median(covs):.2f}% "
      f"(min {min(covs):.2f}), >=99.9%: {full}/{len(rows)}, contigs median {statistics.median(conts):.0f} "
      f"({min(conts)}-{max(conts)}), length median {statistics.median(lens)} ({min(lens)}-{max(lens)})",
      file=sys.stderr)
