#!/usr/bin/env python3
"""Summarise the per-sample HTStream JSON files into one table.

Each JSON is a list of ten objects, one per pipeline stage, in order. The
values pulled out per sample:

  initial fragments        hts_Stats  'initial stats'      Fragment.in
  phiX read pairs          hts_SeqScreener 'screen phix'   Paired_end.hits
  rRNA read pairs          hts_SeqScreener 'count ... rRNA' Paired_end.hits
  PCR duplicates removed   hts_SuperDeduper                Fragment.duplicate
  final fragments          hts_Stats 'final stats'         Fragment.out
  basepairs kept           hts_Stats 'final stats'         Fragment.basepairs_out

Note the rRNA stage runs with -r, so its hits are counted and every read is
passed through. rRNA is measured here, never removed.

usage: hts_summary.py [--hts-dir DIR] [--samples FILE] [--out FILE]
Defaults are read from config.sh via the environment when run through run.sh,
otherwise from the command line.
"""

import argparse
import json
import os
import sys


STAGE_NOTES = {
    "initial": "initial stats",
    "phix": "screen phix",
    "rrna": "count the number of rRNA reads",
    "dedup": "remove PCR duplicates",
    "adapter": "trim adapters",
    "final": "final stats",
}


def stage(entries, note):
    """Return the stage whose notes match, or None."""
    for e in entries:
        if e.get("Program_details", {}).get("options", {}).get("notes") == note:
            return e
    return None


def hits(stage):
    """Screened hits: read pairs for paired runs, reads for single-end runs (the other count is zero)."""
    return stage.get("Paired_end", {}).get("hits", 0) + stage.get("Single_end", {}).get("hits", 0)


def summarise(path):
    with open(path) as fh:
        entries = json.load(fh)
    if not isinstance(entries, list):
        raise ValueError(f"{path}: expected a list of stages, got {type(entries).__name__}")

    found, missing = {}, []
    for key, note in STAGE_NOTES.items():
        st = stage(entries, note)
        if st is None and key != "dedup":  # absent when HTS_DEDUP=no
            missing.append(key)
        else:
            found[key] = st
    if missing:
        raise ValueError(f"{path}: stages not found: {', '.join(missing)}")

    init, phix, rrna = found["initial"], found["phix"], found["rrna"]
    dedup, adapt, final = found.get("dedup"), found["adapter"], found["final"]
    dups = dedup["Fragment"].get("duplicate", 0) if dedup else 0

    n_in = init["Fragment"]["in"]
    bp_in = init["Fragment"]["basepairs_in"]
    n_out = final["Fragment"]["out"]
    bp_out = final["Fragment"]["basepairs_out"]

    return {
        "reads_in": n_in,
        "phix": hits(phix),
        "rrna": hits(rrna),
        "rrna_pct": 100.0 * hits(rrna) / n_in if n_in else 0.0,
        "dups": dups,
        "dups_pct": 100.0 * dups / n_in if n_in else 0.0,
        "adapter_bp_trimmed": adapt["Fragment"]["basepairs_in"] - adapt["Fragment"]["basepairs_out"],
        "reads_out": n_out,
        "reads_kept_pct": 100.0 * n_out / n_in if n_in else 0.0,
        "bp_out": bp_out,
        "bp_kept_pct": 100.0 * bp_out / bp_in if bp_in else 0.0,
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--hts-dir", default=os.environ.get("HTS_OUT", "01-HTS_Preproc"))
    p.add_argument("--samples", default=os.environ.get("SAMPLES_TXT", "samples.txt"))
    p.add_argument("--out", default=None, help="write TSV here as well as stdout")
    a = p.parse_args()

    if not os.path.isfile(a.samples):
        sys.exit(f"samples file not found: {a.samples}")
    with open(a.samples) as fh:
        samples = [ln.strip() for ln in fh if ln.strip()]
    if not samples:
        sys.exit(f"{a.samples} is empty")

    cols = ["sample", "reads_in", "phix", "rrna", "rrna_pct", "dups", "dups_pct",
            "adapter_bp_trimmed", "reads_out", "reads_kept_pct", "bp_out", "bp_kept_pct"]
    rows, problems = [], []

    for s in samples:
        path = os.path.join(a.hts_dir, s, f"{s}.json")
        if not os.path.isfile(path):
            problems.append(f"{s}: no {path}")
            continue
        try:
            r = summarise(path)
        except Exception as exc:                      # noqa: BLE001 - report and continue
            problems.append(f"{s}: {exc}")
            continue
        r["sample"] = s
        rows.append(r)

    def fmt(v):
        if isinstance(v, float):
            return f"{v:.2f}"
        return str(v)

    widths = {c: max(len(c), *(len(fmt(r[c])) for r in rows)) for c in cols} if rows else {}
    if not rows:
        sys.exit("No complete HTStream reports: " + "; ".join(problems))
    print("  ".join(c.rjust(widths[c]) for c in cols))
    for r in rows:
        print("  ".join(fmt(r[c]).rjust(widths[c]) for c in cols))

    if rows:
        tot_in = sum(r["reads_in"] for r in rows)
        tot_out = sum(r["reads_out"] for r in rows)
        print()
        print(f"{len(rows)} samples   {tot_in:,} read pairs in   {tot_out:,} out "
              f"({100.0 * tot_out / tot_in:.1f}% kept)")
        print(f"rRNA   {min(r['rrna_pct'] for r in rows):.2f}% - {max(r['rrna_pct'] for r in rows):.2f}%")
        print(f"dups   {min(r['dups_pct'] for r in rows):.2f}% - {max(r['dups_pct'] for r in rows):.2f}%")
        print(f"phiX   {sum(r['phix'] for r in rows)} read pairs across all samples")

    if a.out:
        with open(a.out, "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in rows:
                fh.write("\t".join(fmt(r[c]) for c in cols) + "\n")
        print(f"\nwrote {a.out}")

    if problems:
        print("\nproblems:", file=sys.stderr)
        for pr in problems:
            print(f"  {pr}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
