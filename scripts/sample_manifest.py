#!/usr/bin/env python3
"""Require one R1 (and, for paired runs, one R2) per sample in every lane and record the exact inputs."""
import argparse
import csv
import gzip
import re
import subprocess
from pathlib import Path


def read_length(path):
    with gzip.open(path, 'rt', encoding='ascii') as stream:
        header, sequence, plus, quality = [stream.readline().rstrip('\r\n') for _ in range(4)]
    if not header.startswith('@') or not plus.startswith('+') or not sequence or len(sequence) != len(quality):
        raise ValueError(f'Invalid first FASTQ record: {path}')
    return len(sequence)


def collect(lanes, regex, suffixes, excluded):
    if len({str(Path(lane).resolve()) for lane in lanes}) != len(lanes):
        raise ValueError('Lane directories are duplicated')
    records, sample_sets = [], []
    for lane in lanes:
        pairs = {}
        for mate, suffix in enumerate(suffixes, 1):
            paths = sorted(Path(lane).glob('*' + suffix))
            bases = [p.name[:-len(suffix)] for p in paths]
            names = subprocess.run(['sed', '-E', regex], input='\n'.join(bases)+'\n',
                                   capture_output=True, text=True, check=True).stdout.splitlines()
            for path, sample in zip(paths, names):
                if sample in excluded:
                    continue
                if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', sample):
                    raise ValueError(f'Unsupported sample name: {sample!r}')
                if not path.is_file() or not path.stat().st_size:
                    raise ValueError(f'Missing or empty FASTQ: {path}')
                key = (sample, mate)
                if key in pairs:
                    raise ValueError(f'Multiple files for {sample} R{mate} in {lane}')
                pairs[key] = path
        samples = {s for s, _ in pairs}
        if not samples:
            raise ValueError(f'No samples in {lane}')
        for sample in sorted(samples):
            if any((sample, mate) not in pairs for mate in range(1, len(suffixes) + 1)):
                raise ValueError(f'Missing mate for {sample} in {lane}')
            r2 = pairs.get((sample, 2))
            records.append([sample, str(Path(lane).resolve()),
                            str(pairs[sample, 1].absolute()), str(r2.absolute()) if r2 else '',
                            pairs[sample, 1].stat().st_size, r2.stat().st_size if r2 else 0])
        sample_sets.append(samples)
    if not sample_sets or any(s != sample_sets[0] for s in sample_sets[1:]):
        raise ValueError('Lane sample sets differ; no samples were dropped')
    def natural(s):
        return [(1,int(x)) if x.isdigit() else (0,x) for x in re.split(r'(\d+)', s)]
    return sorted(sample_sets[0], key=natural), records


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--lane', action='append', required=True)
    p.add_argument('--regex', required=True)
    p.add_argument('--r1', default='_R1_001.fastq.gz')
    p.add_argument('--r2', default='_R2_001.fastq.gz')
    p.add_argument('--layout', choices=['paired', 'single'], default='paired')
    p.add_argument('--exclude', action='append', default=[])
    p.add_argument('--samples', type=Path, required=True)
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--expected', type=Path)
    a = p.parse_args()
    suffixes = [a.r1, a.r2] if a.layout == 'paired' else [a.r1]
    samples, records = collect(a.lane, a.regex, suffixes, set(a.exclude))
    if a.expected:
        expected = a.expected.read_text().splitlines()
        if len(expected) != len(set(expected)) or set(expected) != set(samples):
            raise ValueError('Observed sample names differ from the approved sample list')
    a.samples.parent.mkdir(parents=True, exist_ok=True)
    a.samples.write_text(''.join(s+'\n' for s in samples))
    with a.manifest.open('w') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(['sample', 'lane', 'r1', 'r2', 'r1_bytes', 'r2_bytes'])
        w.writerows(records)
    what = 'complete pairs' if a.layout == 'paired' else 'single-end files'
    print(f'{len(samples)} samples, {len(a.lane)} lanes; {what} and identical sample sets')


if __name__ == '__main__':
    main()
