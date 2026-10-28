#!/usr/bin/env python3
"""Record per-sample strandedness evidence and reject mixed automatic selections."""
import csv,json,sys
from pathlib import Path
star,samples,out=map(Path,sys.argv[1:4]);requested=sys.argv[4]
rows=[]
for s in samples.read_text().splitlines():
    sums=[0,0,0]
    with (star/s/(s+'_ReadsPerGene.out.tab')).open() as f:
        for i,line in enumerate(f):
            if i<4:continue
            fields=line.rstrip().split('\t')
            for j in range(3):sums[j]+=int(fields[j+1])
    u,fwd,rev=sums
    if max(sums)<=0:raise ValueError('No assigned gene counts: '+s)
    inferred=4 if rev>3*fwd else 3 if fwd>3*rev else 2
    rows.append([s,*sums,inferred])
with out.open('w') as f:
    w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['sample','unstranded','forward','reverse','inferred_column']);w.writerows(rows)
choices={r[-1] for r in rows}
if requested=='auto':
    if len(choices)!=1:raise ValueError('Samples have mixed strandedness evidence; review strandedness.tsv and set COUNT_COLUMN explicitly')
    column=choices.pop()
else:
    column=int(requested)
    if column not in (2,3,4):raise ValueError('Invalid COUNT_COLUMN')
print(column)
