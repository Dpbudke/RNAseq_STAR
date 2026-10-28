#!/usr/bin/env python3
"""Compare every matrix cell with its STAR table using the selected count column."""
import argparse
from array import array
import csv
import hashlib
import json
from pathlib import Path


def compare_counts(star_dir, counts_file, samples_file, column):
    if column not in (2, 3, 4):
        raise ValueError('Count column must be 2, 3, or 4')
    samples = Path(samples_file).read_text().splitlines()
    if not samples or any(not s for s in samples) or len(samples) != len(set(samples)):
        raise ValueError('Sample list is empty or contains duplicates')
    dirs = {p.name for p in Path(star_dir).iterdir() if p.is_dir()}
    if dirs != set(samples):
        raise ValueError(f'STAR sample coverage differs: missing={set(samples)-dirs}, extra={dirs-set(samples)}')
    genes, seen = [], set()
    columns = [array('q') for _ in samples]
    with Path(counts_file).open() as f:
        rows = csv.reader(f, delimiter='\t')
        if next(rows, None) != ['gene_id', *samples]:
            raise ValueError('Matrix sample names or order differ from samples.txt')
        for row in rows:
            if len(row) != len(samples)+1 or not row[0] or row[0] in seen:
                raise ValueError('Invalid row length or duplicate gene in matrix')
            genes.append(row[0]); seen.add(row[0])
            for dest, count in zip(columns, row[1:]):
                if not count.isdecimal():
                    raise ValueError('Matrix counts must be nonnegative integers')
                dest.append(int(count))
    if not genes:
        raise ValueError('Matrix has no genes')
    sums = {}
    for sample, counts in zip(samples, columns):
        path = Path(star_dir)/sample/(sample+'_ReadsPerGene.out.tab')
        with path.open() as f:
            rows = csv.reader(f, delimiter='\t')
            for label in ['N_unmapped', 'N_multimapping', 'N_noFeature', 'N_ambiguous']:
                row = next(rows, [])
                if len(row) != 4 or row[0] != label:
                    raise ValueError(f'Invalid STAR summary rows: {sample}')
            n = 0
            for n, row in enumerate(rows, 1):
                if n > len(genes) or len(row) != 4 or row[0] != genes[n-1]:
                    raise ValueError(f'STAR gene order or coverage differs: {sample}')
                if not row[column-1].isdecimal() or int(row[column-1]) != counts[n-1]:
                    raise ValueError(f'Count differs: {sample}, {row[0]}')
            if n != len(genes):
                raise ValueError(f'STAR table is incomplete: {sample}')
        sums[sample] = sum(counts)
    with Path(counts_file).open('rb') as f:
        digest = hashlib.file_digest(f, 'sha256').hexdigest()
    return {'status':'passed', 'samples':len(samples), 'genes':len(genes),
            'cells_compared':len(samples)*len(genes), 'count_column':column,
            'sample_totals':sums, 'matrix_sha256':digest}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--star-dir', type=Path, required=True)
    p.add_argument('--counts', type=Path, required=True)
    p.add_argument('--samples', type=Path, required=True)
    p.add_argument('--column', type=int, choices=(2,3,4), required=True)
    p.add_argument('--report', type=Path)
    a = p.parse_args()
    report = compare_counts(a.star_dir, a.counts, a.samples, a.column)
    if a.report:
        a.report.write_text(json.dumps(report, indent=2)+'\n')
    print(f"Validated {report['cells_compared']} cells: {report['genes']} genes x {report['samples']} samples, STAR column {a.column}")


if __name__ == '__main__':
    main()
