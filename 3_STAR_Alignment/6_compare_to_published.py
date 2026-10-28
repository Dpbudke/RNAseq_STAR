#!/usr/bin/env python3
"""Compare two count matrices.

Reports whether the gene and sample sets match, per-sample totals and their
difference, how many gene x sample cells are exactly equal, Pearson and Spearman
correlation per sample, and the genes with the largest disagreement.

Useful for checking a run against a previously published matrix for the same
study.

  ./run.sh compare -- /path/to/other_counts.txt
"""

import argparse
import sys

import numpy as np
import pandas as pd


def load(path):
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = df.index.astype(str)
    df.columns = [str(c) for c in df.columns]
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mine", required=True)
    ap.add_argument("--published", required=True)
    ap.add_argument("--out")
    ap.add_argument("--top", type=int, default=15)
    a = ap.parse_args()

    out = []

    def say(line=""):
        print(line)
        out.append(line)

    mine = load(a.mine)
    pub = load(a.published)

    say("=== shapes ===")
    say(f"  mine      {mine.shape[0]} genes x {mine.shape[1]} samples   {a.mine}")
    say(f"  published {pub.shape[0]} genes x {pub.shape[1]} samples   {a.published}")
    say()

    say("=== gene IDs ===")
    gm, gp = set(mine.index), set(pub.index)
    if gm == gp:
        say(f"  identical sets, {len(gm)} genes")
    else:
        say(f"  DIFFER: {len(gm - gp)} only in mine, {len(gp - gm)} only in published")
        for g in sorted(gm - gp)[:5]:
            say(f"    only mine:      {g}")
        for g in sorted(gp - gm)[:5]:
            say(f"    only published: {g}")
    say()

    say("=== samples ===")
    sm, sp = set(mine.columns), set(pub.columns)
    if sm == sp:
        say(f"  identical sets, {len(sm)} samples")
    else:
        say(f"  DIFFER: {len(sm - sp)} only in mine, {len(sp - sm)} only in published")
        if sm - sp:
            say(f"    only mine:      {', '.join(sorted(sm - sp))}")
        if sp - sm:
            say(f"    only published: {', '.join(sorted(sp - sm))}")
    say()

    genes = sorted(gm & gp)
    samples = sorted(sm & sp)
    if not genes or not samples:
        say("nothing in common, stopping")
        return 1

    A = mine.loc[genes, samples].to_numpy(dtype=np.int64)
    B = pub.loc[genes, samples].to_numpy(dtype=np.int64)

    say(f"=== comparing {len(genes)} genes x {len(samples)} samples ===")
    total_cells = A.size
    equal_cells = int((A == B).sum())
    say(f"  cells exactly equal: {equal_cells:,} / {total_cells:,} "
        f"({100.0 * equal_cells / total_cells:.4f}%)")
    diff = A.astype(np.int64) - B.astype(np.int64)
    say(f"  max absolute difference in any cell: {int(np.abs(diff).max()):,}")
    say(f"  total counts mine:      {int(A.sum()):,}")
    say(f"  total counts published: {int(B.sum()):,}")
    if B.sum():
        say(f"  overall difference:     {100.0 * (A.sum() - B.sum()) / B.sum():+.4f}%")
    say()

    say("=== per sample ===")
    say(f"  {'sample':<10} {'mine':>14} {'published':>14} {'diff':>12} {'diff %':>9} "
        f"{'cells equal':>12} {'pearson':>9} {'spearman':>9}")
    rows = []
    for j, s in enumerate(samples):
        x, y = A[:, j], B[:, j]
        eq = int((x == y).sum())
        d = int(x.sum() - y.sum())
        dpct = 100.0 * d / y.sum() if y.sum() else float("nan")
        # correlations on the raw counts, and on ranks for spearman
        pear = float(np.corrcoef(x, y)[0, 1]) if x.std() and y.std() else float("nan")
        xr = pd.Series(x).rank().to_numpy()
        yr = pd.Series(y).rank().to_numpy()
        spear = float(np.corrcoef(xr, yr)[0, 1]) if xr.std() and yr.std() else float("nan")
        rows.append((s, int(x.sum()), int(y.sum()), d, dpct, eq, pear, spear))
        say(f"  {s:<10} {x.sum():>14,} {y.sum():>14,} {d:>12,} {dpct:>8.3f}% "
            f"{eq:>12,} {pear:>9.6f} {spear:>9.6f}")
    say()

    pears = [r[6] for r in rows if not np.isnan(r[6])]
    spears = [r[7] for r in rows if not np.isnan(r[7])]
    if pears:
        say(f"  pearson  min {min(pears):.6f}  median {float(np.median(pears)):.6f}  max {max(pears):.6f}")
    if spears:
        say(f"  spearman min {min(spears):.6f}  median {float(np.median(spears)):.6f}  max {max(spears):.6f}")
    say()

    if equal_cells != total_cells:
        say(f"=== {a.top} genes with the largest total absolute disagreement ===")
        per_gene = np.abs(diff).sum(axis=1)
        order = np.argsort(per_gene)[::-1][: a.top]
        say(f"  {'gene':<24} {'sum |diff|':>12} {'mine total':>14} {'pub total':>14}")
        for i in order:
            if per_gene[i] == 0:
                break
            say(f"  {genes[i]:<24} {int(per_gene[i]):>12,} "
                f"{int(A[i].sum()):>14,} {int(B[i].sum()):>14,}")
        say()

    if a.out:
        with open(a.out, "w") as fh:
            fh.write("\n".join(out) + "\n")
        print(f"wrote {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
