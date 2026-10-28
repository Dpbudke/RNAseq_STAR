#!/usr/bin/env python3
"""Select exact lane files from the validated input manifest."""

import csv
import sys
from pathlib import Path


def select(manifest, sample, mate, lanes):
    if mate not in {"R1", "R2"}:
        raise ValueError("Mate must be R1 or R2")
    lane_order = [str(Path(lane).resolve()) for lane in lanes]
    if len(set(lane_order)) != len(lane_order):
        raise ValueError("Lane directories are duplicated")
    field = mate.lower()
    with Path(manifest).open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if reader.fieldnames != ["sample", "lane", "r1", "r2", "r1_bytes", "r2_bytes"]:
            raise ValueError("Input manifest headers differ")
        rows = [row for row in reader if row["sample"] == sample]
    if len(rows) != len(lane_order) or {row["lane"] for row in rows} != set(lane_order):
        raise ValueError("Sample does not have exactly one manifest row per lane")
    by_lane = {row["lane"]: row for row in rows}
    paths = []
    for lane in lane_order:
        row = by_lane[lane]
        if not row[field]:
            raise ValueError(f"The input manifest has no {mate} files (single-end run)")
        path = Path(row[field])
        if (
            not path.is_absolute()
            or str(path.parent.resolve()) != lane
            or not path.is_file()
            or path.stat().st_size != int(row[field + "_bytes"])
            or "\n" in str(path)
            or "\r" in str(path)
        ):
            raise ValueError("Manifest input is missing, changed or outside its lane")
        paths.append(str(path))
    return paths


if __name__ == "__main__":
    print("\n".join(select(*sys.argv[1:4], sys.argv[4:])))
