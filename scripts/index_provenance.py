#!/usr/bin/env python3
"""Bind a complete STAR index to reference hashes, version, and overhang."""
import hashlib,json,sys
from pathlib import Path

def main():
    index,manifest=map(Path,sys.argv[1:3]);overhang,version,mode=sys.argv[3:]
    required=['SA','SAindex','Genome','chrName.txt','chrLength.txt','chrStart.txt','exonInfo.tab','geneInfo.tab','sjdbList.out.tab','genomeParameters.txt']
    for name in required:
        if not (index/name).is_file() or not (index/name).stat().st_size:raise ValueError('Incomplete STAR index: '+name)
    expected={'reference_manifest_sha256':hashlib.sha256(manifest.read_bytes()).hexdigest(),'overhang':int(overhang),'star_version':version}
    report=index/'index_manifest.json'
    if mode=='record':
        expected['files']={}
        for name in required:
            with (index/name).open('rb') as f:h=hashlib.file_digest(f,'sha256').hexdigest()
            expected['files'][name]={'sha256':h,'bytes':(index/name).stat().st_size}
        report.write_text(json.dumps(expected,indent=2)+'\n')
    else:
        old=json.loads(report.read_text())
        if any(old.get(k)!=v for k,v in expected.items()):raise ValueError('STAR index provenance differs; use a new index directory')
        for name in required:
            with (index/name).open('rb') as f:h=hashlib.file_digest(f,'sha256').hexdigest()
            if old['files'][name]['sha256']!=h:raise ValueError('STAR index content changed: '+name)
    print('STAR index verified:',report)

if __name__=='__main__':main()
