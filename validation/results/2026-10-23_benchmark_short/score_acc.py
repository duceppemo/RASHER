#!/usr/bin/env python3
"""Cross-compare SNP distance matrices from runs of the same real dataset under different tools.

Usage: python score_acc.py LABEL=path/to/snp_distances.tsv [LABEL=path ...]

Sample names are normalized by stripping a trailing _ERR<digits> run id (the Illumina runs carry it, the
ONT runs do not); the Reference row/column is dropped. Prints each run's distance summary, then for every
pair of runs: shared sample pairs, identical pairwise distances, mean/max absolute difference, and the
Pearson correlation.
"""

from __future__ import annotations

import re
import statistics
import sys
from itertools import combinations
from pathlib import Path


def norm(name: str) -> str:
    return re.sub(r"_ERR\d+$", "", name)


def matrix(path: str) -> dict[tuple[str, str], int]:
    lines = Path(path).read_text().splitlines()
    header = [norm(h) for h in lines[0].split("\t")]
    out: dict[tuple[str, str], int] = {}
    for line in lines[1:]:
        cells = line.split("\t")
        row = norm(cells[0])
        for col, val in zip(header[1:], cells[1:]):
            if row != col and "Reference" not in (row, col):
                out[tuple(sorted((row, col)))] = int(val)
    return out


def pearson(x: list[int], y: list[int]) -> float:
    if len(x) < 3 or statistics.pstdev(x) == 0 or statistics.pstdev(y) == 0:
        return float("nan")
    mx, my = statistics.fmean(x), statistics.fmean(y)
    num = sum((a - mx) * (b - my) for a, b in zip(x, y))
    den = (sum((a - mx) ** 2 for a in x) * sum((b - my) ** 2 for b in y)) ** 0.5
    return num / den


def main() -> None:
    runs = dict(arg.split("=", 1) for arg in sys.argv[1:])
    mats = {}
    for label, path in runs.items():
        mats[label] = matrix(path)
        vals = list(mats[label].values())
        samples = {s for pair in mats[label] for s in pair}
        print(
            f"{label}: {len(samples)} samples, median pairwise distance {statistics.median(vals):.0f}, "
            f"range {min(vals)}-{max(vals)}"
        )
    print()
    for a, b in combinations(mats, 2):
        common = sorted(set(mats[a]) & set(mats[b]))
        if not common:
            print(f"{a}  vs  {b}: no samples in common")
            continue
        xa = [mats[a][k] for k in common]
        xb = [mats[b][k] for k in common]
        same = sum(1 for u, v in zip(xa, xb) if u == v)
        diffs = [abs(u - v) for u, v in zip(xa, xb)]
        print(f"{a}  vs  {b}:")
        print(
            f"  {len(common)} shared pairs; identical {same}/{len(common)}; "
            f"mean |diff| {statistics.fmean(diffs):.1f}; max |diff| {max(diffs)}; r={pearson(xa, xb):.4f}"
        )


if __name__ == "__main__":
    main()
