# 1_Pre_Preprocessing

Lists the samples and combines each sample's sequencing lanes into one FASTQ per read before QC. `run.sh` runs both
steps (`./run.sh samples`, `./run.sh concat`, or `./run.sh all`).

## Scripts

### `1_make_samples_txt.sh`
Finds every sample in every lane directory by its full file name and records the exact input files.

- **Input:** `LANE_DIRS` (one directory per lane) with `<sample>_S<n>_L<lane>_R1_001.fastq.gz` and, for paired-end
  runs, the matching `R2`; `SAMPLE_NAME_STYLE` sets how a file name becomes a sample name (`./run.sh scan` shows it)
- **Output:** `samples.txt` (one sample per line) and `input_manifest.tsv` (every sample's file in every lane, with
  its size)
- **Checks:** every lane holds the same samples, each with all its reads; names that begin other names (`MOUSE_1`,
  `MOUSE_10`) keep their own files; `EXPECTED_SAMPLES_FILE`, when set, must list exactly the samples found

### `2_concatenate_rnaseq.slurm`
Slurm array job, one task per sample: joins that sample's lane files listed in `input_manifest.tsv`.

- **Output:** `concatenated_reads/<sample>/<sample>_R1_concatenated.fastq.gz` (and `_R2_` for paired-end runs)
- **Checks:** each joined file decompresses cleanly and holds whole FASTQ records; R1 and R2 hold the same number of
  reads

## Expected output structure

```
concatenated_reads/
├── Sample1/
│   ├── Sample1_R1_concatenated.fastq.gz
│   └── Sample1_R2_concatenated.fastq.gz
├── Sample2/
│   └── ...
```
