#!/usr/bin/env python3
"""Collect STAR mapping metrics and check every BAM header and end-of-file block."""
import csv,json,subprocess,sys
from pathlib import Path
work=Path(sys.argv[1]);samples=(work/'samples.txt').read_text().splitlines();rows=[];files={}
fields=['Number of input reads','Uniquely mapped reads number','Uniquely mapped reads %','Number of reads mapped to multiple loci','% of reads mapped to multiple loci','% of reads unmapped: too short']
for s in samples:
    d=work/'02-STAR_alignment'/s;p=d/(s+'_Log.final.out')
    values={}
    for line in p.read_text().splitlines():
        if '|' in line:
            k,v=line.split('|',1);values[k.strip()]=v.strip()
    if any(k not in values for k in fields) or int(values['Number of input reads'])<=0:
        raise ValueError('Missing mapping metrics or no input reads: '+s)
    bam=d/(s+'_Aligned.sortedByCoord.out.bam')
    subprocess.run(['samtools','quickcheck','-v',str(bam)],check=True)
    subprocess.run(['samtools','index','-@','2',str(bam)],check=True)
    for path in [bam,Path(str(bam)+'.bai')]:
        stat=path.stat()
        if stat.st_size<=0:raise ValueError('Empty alignment or index: '+str(path))
        files[str(path.relative_to(work))]=[stat.st_size,stat.st_mtime_ns,stat.st_ctime_ns,stat.st_ino]
    rows.append([s,*[values[k] for k in fields]])
with (work/'mapping_qc.tsv').open('w') as f:
    w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['sample',*fields]);w.writerows(rows)
(work/'bam_validation.json').write_text(json.dumps({'status':'passed','samples':len(samples),'checks':['samtools quickcheck','coordinate BAM index'],'fingerprint_fields':['bytes','mtime_ns','ctime_ns','inode'],'files':files},indent=2)+'\n')
print('Mapping metrics and BAM validation:',len(samples),'samples')
