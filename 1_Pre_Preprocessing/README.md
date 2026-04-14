# 1_Pre_Preprocessing

Merges multi-lane FASTQs into per-sample files and organizes them into individual directories before QC.

## Scripts

### `1_concatenate_rnaseq.sh`
Concatenates R1 and R2 reads from two sequencing lanes into a single file per read per sample.

- **Input:** Two lane directories containing `*_L008_R{1,2}_001.fastq.gz` and `*_L004_R{1,2}_001.fastq.gz`
- **Output:** `concatenated_reads/<sample>_R{1,2}_concatenated.fastq.gz`
- **Usage:** Edit `LANE1_DIR`, `LANE2_DIR`, and `OUTPUT_DIR` at the top of the script, then run:
  ```bash
  bash 1_concatenate_rnaseq.sh
  ```

### `2_organize_samples.sh`
Moves concatenated FASTQ pairs into per-sample subdirectories.

- **Input:** `*_R1_concatenated.fastq.gz` / `*_R2_concatenated.fastq.gz` in current directory
- **Output:** `<sample>/` directories each containing the R1 and R2 files
- **Usage:**
  ```bash
  cd concatenated_reads
  bash ../2_organize_samples.sh
  ```

> **After running:** move `2_organize_samples.sh` out of `concatenated_reads/` back to the base directory (e.g. `CCPups/`) before generating `samples.txt`, so the script filename doesn't appear as a sample name.

### Generating `samples.txt`

`samples.txt` (one sample name per line) is required by the HTStream and STAR SLURM array jobs. Generate it from inside `concatenated_reads/` after organizing:

```bash
cd concatenated_reads
ls > ../samples.txt
cat ../samples.txt   # verify contents look correct
```

## Expected Output Structure

```
concatenated_reads/
├── Sample1/
│   ├── Sample1_R1_concatenated.fastq.gz
│   └── Sample1_R2_concatenated.fastq.gz
├── Sample2/
│   └── ...
```
