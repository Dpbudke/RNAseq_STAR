# 3_STAR_Alignment

STAR-based read alignment to a reference genome, followed by read count matrix construction for downstream DE analysis.

## Scripts

### `1_star_index.slurm`
Builds the STAR genome index (run once per reference genome).

- **Input:** GRCm39 genome FASTA + GENCODE M35 GTF annotation
- **Output:** `References/star.overlap100.gencode.M35/` index directory
- **Usage:**
  ```bash
  mkdir -p slurmout
  sbatch 1_star_index.slurm
  ```
- **Runtime:** ~1–2 hours, 8 cores, 128 GB RAM

---

### `2_star.slurm`
SLURM array job — aligns each sample's trimmed reads to the STAR index.

- **Input:** `01-HTS_Preproc/<sample>/` trimmed FASTQs; `samples.txt`
- **Output:** `02-STAR_alignment/<sample>/` — BAM files + `*_ReadsPerGene.out.tab` count tables
- **Usage:**
  ```bash
  # Edit --array=1-N to match number of samples
  sbatch 2_star.slurm
  ```

---

### `3_star_stats.sh`
Summarizes alignment statistics across all samples into a single tab-separated table.

- **Input:** `02-STAR_alignment/*/` (all STAR output directories)
- **Output:** `summary_star_alignments.txt`
- **Usage:**
  ```bash
  bash 3_star_stats.sh
  ```

---

### `4_reads2counts.sh`
Aggregates per-sample `ReadsPerGene.out.tab` files into a single gene × sample count matrix.

- **Input:** `02-STAR_alignment/` (STAR output)
- **Output:** `03-Counts/` directory with count matrix
- **Usage:**
  ```bash
  bash 4_reads2counts.sh
  ```

---

### `5_confirm_reads2counts.py`
Validates the count matrix by cross-checking totals between the STAR output and the merged count file.

- **Usage:**
  ```bash
  python3 5_confirm_reads2counts.py
  ```

## Expected Output Structure

```
02-STAR_alignment/
├── Sample1/
│   ├── Sample1_Aligned.sortedByCoord.out.bam
│   └── Sample1_ReadsPerGene.out.tab
03-Counts/
└── <prefix>_counts.txt    # gene × sample count matrix
summary_star_alignments.txt
```
