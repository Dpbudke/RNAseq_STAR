# Check a completed run

Use the same pipeline copy and run configuration that created the results. With
`EXECUTOR=slurm`, `confirm` and the other checks below submit a small Slurm job (a
login node is only for submitting and checking jobs); read its output in
`WORKDIR/slurmout/` once `sacct` shows it finished.

```bash
./run.sh status
./run.sh confirm
```

`confirm` compares every gene/sample cell against the chosen STAR count column.
Its report is `WORKDIR/count_validation.json`. All sample names, gene IDs and
integer counts must agree. Matching totals alone does not establish equality.

Check every recorded Slurm job with `sacct -X -j JOB_IDS` and require `COMPLETED`
with exit `0:0` for every expected array task. In local mode, inspect `local_jobs.json`:
all requested tasks' latest attempts must have completed with exit 0. Older failed
attempts stay in that journal. A saved status report does not replace output checks.

Finalization writes:

- `count_validation.json`: dimensions, selected column and matrix SHA-256.
- `bam_validation.json`: quickcheck/index results and file identities for every sample.
- `mapping_qc.tsv`, `summary_htstream.txt`, `strandedness.tsv`: per-sample metrics.
- `pipeline_state.txt`: `validated` after the final checks pass.

To repeat BAM checks and indexing from the pipeline directory:

```bash
source config.sh
load_software "$MODULE_PYTHON" "$MODULE_SAMTOOLS"
python3 scripts/final_qc.py "$WORKDIR"
```

This reads BAMs and may rebuild indexes. Save the existing validation reports before
repeating a check that writes them. Review unusual QC measurements with the study's
library protocol; there is no universal mapping or rRNA percentage that proves a
biologically sound experiment.

`run_records/code_sha256.txt` binds the original source paths and selected profile.
Every later stage checks those hashes. `run_records/code/` preserves the source
snapshot; `settings_snapshot.sh` records selected effective settings, and
`software.json` identifies executables. The original source paths remain historical
if a completed result directory is moved. Configure a new WORKDIR for a rerun.

Keep source MD5 reports, samples/input_manifest, count tables, references and index
hashes with the results. A copy to permanent storage needs an explicit file manifest
and checksum comparison at the destination. Gzip checks prove stream integrity;
they do not prove another copy contains the same bytes. Slurm creates .err files
even for successful jobs; read nonempty messages together with exit codes and checks.


## Completed checks, September 10, 2026

The earlier recorded Slurm workflow completed a 264-sample, two-lane Ceres run:
57,132 genes, every count cell checked, and all BAMs checked and indexed.

The later reusable version completed local and Slurm runs on the same two samples,
two lanes and 20,000 read pairs per sample/lane, using custom sample names. Each
mode completed 12 tasks. Count matrices, all STAR count columns and QC matched;
BAM checks and the source-change guard passed. Both modes ran on Ceres Linux
compute nodes with the same GRCm39/M35 index and pinned STAR executable. These
checks do not establish compatibility with every Linux distribution or tool build.

The 30 automated tests also passed. Repeat them from the repository directory:

```bash
python3 -B -m unittest discover -s tests -v
```
