#!/usr/bin/env python3
"""Join the count tables of several finished runs into one table.

Use this when a study's libraries were sequenced as separate pools on different
lanes. Each pool runs as its own pipeline run (lanes inside one run must hold the
same samples); this joins the validated results afterwards.

Every run must have finished finalization (pipeline_state.txt says validated),
its counts.txt must still match the SHA-256 in count_validation.json, all runs
must list the same genes in the same order and use the same STAR count column,
and no sample may appear in more than one run. The output directory must not
exist. After writing, every column is read back and compared with its source.

    python3 scripts/combine_runs.py OUT_DIR RUN_DIR [RUN_DIR ...]

Writes OUT_DIR/counts.txt, OUT_DIR/samples.tsv (sample, run directory, reads
counted) and OUT_DIR/SHA256SUMS.txt.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def read_matrix(path):
    with path.open(newline="") as stream:
        rows = list(csv.reader(stream, delimiter="\t"))
    if not rows or rows[0][:1] != ["gene_id"]:
        raise ValueError(f"{path}: first column is not gene_id")
    header, body = rows[0], rows[1:]
    if any(len(row) != len(header) for row in body):
        raise ValueError(f"{path}: rows have different lengths")
    genes = [row[0] for row in body]
    columns = {name: [row[i] for row in body] for i, name in enumerate(header) if i}
    return header[1:], genes, columns


def load_run(run):
    run = Path(run).resolve()
    state = run / "pipeline_state.txt"
    if not state.is_file() or state.read_text().strip() != "validated":
        raise ValueError(f"{run}: finalization has not passed (pipeline_state.txt)")
    report = json.loads((run / "count_validation.json").read_text())
    counts = run / "03-Counts" / "counts.txt"
    if report.get("status") != "passed" or report.get("matrix_sha256") != digest(counts):
        raise ValueError(f"{run}: count validation did not pass or counts.txt changed since")
    samples, genes, columns = read_matrix(counts)
    return run, report["count_column"], samples, genes, columns


def combine(out_dir, run_dirs):
    out_dir = Path(out_dir)
    if out_dir.exists():
        raise ValueError(f"{out_dir} already exists")
    if len(run_dirs) < 2:
        raise ValueError("Give at least two run directories")
    runs = [load_run(r) for r in run_dirs]
    if len({r[0] for r in runs}) != len(runs):
        raise ValueError("A run directory is listed twice")
    genes = runs[0][3]
    if len({r[1] for r in runs}) != 1:
        raise ValueError("Runs used different STAR count columns; check the library preparation")
    order, columns, sources = [], {}, {}
    for run, _, samples, run_genes, run_columns in runs:
        if run_genes != genes:
            raise ValueError(f"{run}: gene IDs or their order differ; use the same reference for every run")
        for sample in samples:
            if sample in columns:
                raise ValueError(f"{sample} is in more than one run")
            order.append(sample)
            columns[sample] = run_columns[sample]
            sources[sample] = run
    out_dir.mkdir(parents=True)
    with (out_dir / "counts.txt").open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene_id", *order])
        for i, gene in enumerate(genes):
            writer.writerow([gene, *(columns[s][i] for s in order)])
    back_samples, back_genes, back = read_matrix(out_dir / "counts.txt")
    if back_samples != order or back_genes != genes or any(back[s] != columns[s] for s in order):
        raise ValueError("Read-back of the joined table differs from its sources")
    with (out_dir / "samples.tsv").open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample", "run_dir", "reads_counted"])
        for sample in order:
            writer.writerow([sample, sources[sample], sum(int(x) for x in columns[sample])])
    sums = "".join(f"{digest(out_dir / n)}  {n}\n" for n in ("counts.txt", "samples.tsv"))
    (out_dir / "SHA256SUMS.txt").write_text(sums)
    return {"samples": len(order), "genes": len(genes), "runs": len(runs), "count_column": runs[0][1]}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("out_dir")
    parser.add_argument("run_dirs", nargs="+")
    args = parser.parse_args()
    print(json.dumps(combine(args.out_dir, args.run_dirs)))


if __name__ == "__main__":
    main()
