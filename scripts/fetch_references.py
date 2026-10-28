#!/usr/bin/env python3
"""Fetch a pinned genome and annotation into a reference folder, check them, and record SHA-256.

Sources (--source):
  gencode_mouse  GENCODE mouse release (--release M35). File names, manifest and rRNA file are those of the
                 lab's mouse runs, so existing reference folders and STAR indexes stay valid.
  gencode_human  GENCODE human release (--release 47): GRCh38 primary assembly and its annotation.
  ensembl        Ensembl release for any species (--species rattus_norvegicus --release 113): the primary
                 assembly (or toplevel when there is none) and the release's GTF.
  custom         Your own FASTA and GTF (--genome-fa, --gtf; plain or gzip).

Downloads are checked against the provider's published checksums (GENCODE MD5SUMS, Ensembl CHECKSUMS) and every
file's SHA-256 is recorded in reference_manifest.json. A folder already holding another reference, or a file whose
SHA-256 no longer matches its record, is refused. The rRNA sequences for the HTStream rRNA count are the repository's
mouse file for gencode_mouse, --rrna when given, and otherwise the annotation's rRNA genes (an estimate for QC only;
it does not change gene counts). Finally the genome and annotation must share chromosome names.
"""
import argparse
import gzip
import hashlib
import json
import re
import shutil
import subprocess
import urllib.request
from pathlib import Path

GENCODE = 'https://ftp.ebi.ac.uk/pub/databases/gencode'
ENSEMBL = 'https://ftp.ensembl.org/pub'
RRNA_TYPES = {'rRNA', 'Mt_rRNA'}


def digest(path, algorithm='sha256'):
    with path.open('rb') as f:
        return hashlib.file_digest(f, algorithm).hexdigest()


def fetch_text(url):
    return urllib.request.urlopen(url, timeout=60).read().decode('ascii')


def gencode_mouse(a):
    """The lab's mouse references, unchanged: same files, same manifest bytes, same rRNA file."""
    base = f'{GENCODE}/Gencode_mouse/release_{a.release}'
    manifest = fetch_text(base + '/MD5SUMS')
    md5 = {line.split()[1].lstrip('*'): line.split()[0] for line in manifest.splitlines() if len(line.split()) == 2}
    report = a.directory / 'reference_manifest.json'
    old = json.loads(report.read_text()) if report.exists() else {}
    recorded = old.get('files', {}) if old.get('release') == a.release else {}
    outputs = {}
    for name in ['GRCm39.primary_assembly.genome.fa', f'gencode.v{a.release}.annotation.gtf']:
        dest = a.directory / name
        expected = md5[name + '.gz']
        previous = recorded.get(name, {})
        if dest.exists():
            if previous.get('source_md5') != expected or previous.get('sha256') != digest(dest):
                raise ValueError(f'Existing reference lacks matching provenance: {dest}')
        else:
            part = dest.with_name(name + '.download.part')
            unpacked = dest.with_name(name + '.unpacked.part')
            with urllib.request.urlopen(base + '/' + name + '.gz', timeout=120) as src, part.open('wb') as dst:
                shutil.copyfileobj(src, dst, 8 * 1024 * 1024)
            if digest(part, 'md5') != expected:
                raise ValueError(f'Published MD5 mismatch: {name}')
            with gzip.open(part) as src, unpacked.open('wb') as dst:
                shutil.copyfileobj(src, dst, 8 * 1024 * 1024)
            unpacked.replace(dest)
            part.unlink()
        outputs[name] = {'source_url': base + '/' + name + '.gz', 'source_md5': expected, 'sha256': digest(dest),
                         'bytes': dest.stat().st_size}
        report.write_text(json.dumps({'assembly': 'GRCm39', 'release': a.release, 'files': {**recorded, **outputs}},
                                     indent=2) + '\n')
    rrna = a.directory / 'mouse_rrna.fasta'
    if rrna.exists() and digest(rrna) != digest(a.rrna):
        raise ValueError('Existing rRNA reference differs from the repository copy')
    if not rrna.exists():
        shutil.copyfile(a.rrna, rrna)
    outputs[rrna.name] = {'source': 'repository References/mouse_rrna.fasta', 'sha256': digest(rrna),
                          'bytes': rrna.stat().st_size}
    report.write_text(json.dumps({'assembly': 'GRCm39', 'release': a.release, 'files': outputs}, indent=2) + '\n')
    (a.directory / 'GENCODE_MD5SUMS.txt').write_text(manifest)
    return a.directory / 'GRCm39.primary_assembly.genome.fa', a.directory / f'gencode.v{a.release}.annotation.gtf'


def bsd_sum(path):
    """Ensembl's CHECKSUMS values: the BSD sum (GNU `sum` default) and 1 KiB block count."""
    out = subprocess.run(['sum', str(path)], check=True, capture_output=True, text=True).stdout.split()
    return int(out[0]), int(out[1])


def download(url, dest, check):
    """Download url, run check(part) on the compressed file, then unpack gzip into dest."""
    part = dest.with_name(dest.name + '.download.part')
    unpacked = dest.with_name(dest.name + '.unpacked.part')
    with urllib.request.urlopen(url, timeout=120) as src, part.open('wb') as dst:
        shutil.copyfileobj(src, dst, 8 * 1024 * 1024)
    check(part)
    with gzip.open(part) as src, unpacked.open('wb') as dst:
        shutil.copyfileobj(src, dst, 8 * 1024 * 1024)
    unpacked.replace(dest)
    part.unlink()


def copy_local(source, dest):
    """Copy a user's file into the reference folder, unpacking gzip."""
    part = dest.with_name(dest.name + '.copy.part')
    with source.open('rb') as f:
        compressed = f.read(2) == b'\x1f\x8b'
    with (gzip.open(source) if compressed else source.open('rb')) as src, part.open('wb') as dst:
        shutil.copyfileobj(src, dst, 8 * 1024 * 1024)
    part.replace(dest)


def plan(a):
    """What to fetch: identity of the reference, and for each output file its source and checker."""
    if a.source == 'gencode_human':
        base = f'{GENCODE}/Gencode_human/release_{a.release}'
        listing = fetch_text(base + '/MD5SUMS')
        md5 = {line.split()[1].lstrip('*'): line.split()[0] for line in listing.splitlines() if len(line.split()) == 2}
        names = {'genome.fa': 'GRCh38.primary_assembly.genome.fa.gz',
                 'annotation.gtf': f'gencode.v{a.release}.annotation.gtf.gz'}
        files = {}
        for out, name in names.items():
            if name not in md5:
                raise ValueError(f'{name} is not in GENCODE human release {a.release}')
            files[out] = {'source_url': f'{base}/{name}', 'source_md5': md5[name]}
        ident = {'source': 'gencode_human', 'release': a.release, 'assembly': 'GRCh38'}
        return ident, files, {'GENCODE_MD5SUMS.txt': listing}
    if a.source == 'ensembl':
        if not re.fullmatch(r'[a-z0-9_]+', a.species or ''):
            raise ValueError('--species takes the Ensembl folder name, for example rattus_norvegicus')
        dna = f'{ENSEMBL}/release-{a.release}/fasta/{a.species}/dna'
        gtf = f'{ENSEMBL}/release-{a.release}/gtf/{a.species}'
        dna_sums, gtf_sums = fetch_text(dna + '/CHECKSUMS'), fetch_text(gtf + '/CHECKSUMS')

        def table(text):
            return {f[2]: (int(f[0]), int(f[1])) for f in (line.split() for line in text.splitlines()) if len(f) == 3}
        dna_t, gtf_t = table(dna_sums), table(gtf_sums)
        genome = ([n for n in dna_t if n.endswith('.dna.primary_assembly.fa.gz')]
                  or [n for n in dna_t if n.endswith('.dna.toplevel.fa.gz')])
        annotation = [n for n in gtf_t
                      if n.endswith(f'.{a.release}.gtf.gz') and '.chr' not in n and 'abinitio' not in n]
        if len(genome) != 1 or len(annotation) != 1:
            raise ValueError(f'Cannot pick one genome and one GTF for {a.species} release {a.release}: '
                             f'{genome} {annotation}')
        assembly = genome[0].split('.')[1]
        files = {'genome.fa': {'source_url': f'{dna}/{genome[0]}', 'source_bsd_sum': list(dna_t[genome[0]])},
                 'annotation.gtf': {'source_url': f'{gtf}/{annotation[0]}',
                                    'source_bsd_sum': list(gtf_t[annotation[0]])}}
        ident = {'source': 'ensembl', 'species': a.species, 'release': a.release, 'assembly': assembly}
        return ident, files, {'ENSEMBL_CHECKSUMS_dna.txt': dna_sums, 'ENSEMBL_CHECKSUMS_gtf.txt': gtf_sums}
    if a.source == 'custom':
        if not (a.genome_fa and a.gtf and a.name):
            raise ValueError('custom needs --name, --genome-fa and --gtf')
        files = {}
        for out, src in (('genome.fa', a.genome_fa), ('annotation.gtf', a.gtf)):
            src = src.resolve()
            files[out] = {'source_path': str(src), 'source_sha256': digest(src)}
        ident = {'source': 'custom', 'name': a.name}
        return ident, files, {}
    raise ValueError(f'Unknown --source {a.source}')


def rrna_from_annotation(genome, gtf, dest):
    """Write the exons of the annotation's rRNA genes (rRNA, Mt_rRNA) as FASTA; returns the record count."""
    wanted = {}
    with gtf.open() as f:
        for line in f:
            if line.startswith('#'):
                continue
            col = line.rstrip('\n').split('\t')
            if len(col) < 9 or col[2] != 'exon':
                continue
            kind = re.search(r'gene_(?:type|biotype) "([^"]+)"', col[8])
            if not kind or kind.group(1) not in RRNA_TYPES:
                continue
            gene = re.search(r'gene_id "([^"]+)"', col[8])
            gene = gene.group(1) if gene else 'rRNA'
            wanted.setdefault(col[0], set()).add((int(col[3]), int(col[4]), col[6], gene))
    complement = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
    count = 0
    part = dest.with_name(dest.name + '.part')
    with genome.open() as f, part.open('w') as out:
        def flush(name, chunks):
            nonlocal count
            if name not in wanted:
                return
            seq = ''.join(chunks)
            for start, end, strand, gene in sorted(wanted[name]):
                s = seq[start - 1:end]
                if strand == '-':
                    s = s.translate(complement)[::-1]
                out.write(f'>{gene}|{name}:{start}-{end}({strand})\n')
                out.write('\n'.join(s[i:i + 60] for i in range(0, len(s), 60)) + '\n')
                count += 1
        name, chunks = None, []
        for line in f:
            if line.startswith('>'):
                flush(name, chunks)
                name, chunks = line[1:].split()[0], []
            elif name in wanted:
                chunks.append(line.strip())
        flush(name, chunks)
    if count == 0:
        part.unlink()
        raise ValueError('The annotation has no rRNA genes; give an rRNA FASTA with --rrna')
    part.replace(dest)
    return count


def check_contigs(genome, gtf):
    """The genome and the annotation must name chromosomes the same way, or STAR counts nothing."""
    with genome.open() as f:
        fasta = {line[1:].split()[0] for line in f if line.startswith('>')}
    annotated = set()
    with gtf.open() as f:
        for line in f:
            if not line.startswith('#'):
                annotated.add(line.split('\t', 1)[0])
    missing = annotated - fasta
    if not annotated & fasta:
        raise ValueError(f'No annotated chromosome is in the genome FASTA (for example {sorted(annotated)[:3]})')
    if missing:
        print(f'note: {len(missing)} annotated sequences are not in the genome FASTA, '
              f'for example {sorted(missing)[:3]}')


def other_source(a):
    ident, files, side = plan(a)
    report = a.directory / 'reference_manifest.json'
    old = json.loads(report.read_text()) if report.exists() else {}
    if old and {k: old.get(k) for k in ident} != ident:
        held = {k: old.get(k) for k in ident}
        raise ValueError(f'{a.directory} already holds another reference ({held}); use a new folder')
    recorded = old.get('files', {})
    outputs = {}
    for out, info in files.items():
        dest = a.directory / out
        previous = recorded.get(out, {})
        if dest.exists():
            same_source = all(previous.get(k) == v for k, v in info.items())
            if not same_source or previous.get('sha256') != digest(dest):
                raise ValueError(f'Existing reference lacks matching provenance: {dest}')
        elif 'source_path' in info:
            copy_local(Path(info['source_path']), dest)
        else:
            def check(part, info=info, out=out):
                if 'source_md5' in info and digest(part, 'md5') != info['source_md5']:
                    raise ValueError(f'Published MD5 mismatch: {out}')
                if 'source_bsd_sum' in info and list(bsd_sum(part)) != info['source_bsd_sum']:
                    raise ValueError(f'Published Ensembl checksum mismatch: {out}')
            download(info['source_url'], dest, check)
        outputs[out] = {**info, 'sha256': digest(dest), 'bytes': dest.stat().st_size}
        report.write_text(json.dumps({**ident, 'files': {**recorded, **outputs}}, indent=2) + '\n')
    genome, gtf = a.directory / 'genome.fa', a.directory / 'annotation.gtf'
    rrna = a.directory / 'rrna.fasta'
    if a.rrna:
        info = {'source': f'user file {a.rrna.resolve()}', 'source_sha256': digest(a.rrna)}
        if rrna.exists() and digest(rrna) != digest(a.rrna):
            raise ValueError('Existing rrna.fasta differs from --rrna')
        if not rrna.exists():
            shutil.copyfile(a.rrna, rrna)
    else:
        info = {'source': 'exons of the annotation genes with gene type rRNA or Mt_rRNA'}
        previous = recorded.get(rrna.name, {})
        if rrna.exists() and previous.get('sha256') != digest(rrna):
            raise ValueError(f'Existing reference lacks matching provenance: {rrna}')
        if not rrna.exists():
            info['records'] = rrna_from_annotation(genome, gtf, rrna)
        else:
            info['records'] = previous.get('records')
    outputs[rrna.name] = {**info, 'sha256': digest(rrna), 'bytes': rrna.stat().st_size}
    report.write_text(json.dumps({**ident, 'files': outputs}, indent=2) + '\n')
    for name, text in side.items():
        (a.directory / name).write_text(text)
    return genome, gtf


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--directory', type=Path, required=True)
    p.add_argument('--source', default='gencode_mouse', choices=['gencode_mouse', 'gencode_human', 'ensembl', 'custom'])
    p.add_argument('--release', help='GENCODE release (M35, 47) or Ensembl release (113)')
    p.add_argument('--species', help='ensembl: the species folder name, for example rattus_norvegicus')
    p.add_argument('--name', help='custom: a short name for this reference')
    p.add_argument('--genome-fa', type=Path, help='custom: genome FASTA (plain or gzip)')
    p.add_argument('--gtf', type=Path, help='custom: annotation GTF (plain or gzip)')
    p.add_argument('--rrna', type=Path, help='rRNA FASTA (gencode_mouse: required, the repository file)')
    a = p.parse_args()
    a.directory.mkdir(parents=True, exist_ok=True)
    if a.source in ('gencode_mouse', 'gencode_human', 'ensembl') and not a.release:
        p.error(f'{a.source} needs --release')
    if a.source == 'gencode_mouse':
        if not a.rrna:
            p.error('gencode_mouse needs --rrna (References/mouse_rrna.fasta)')
        genome, gtf = gencode_mouse(a)
    else:
        genome, gtf = other_source(a)
    check_contigs(genome, gtf)
    print('References verified:', a.directory / 'reference_manifest.json')


if __name__ == '__main__':
    main()
