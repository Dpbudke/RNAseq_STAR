# 3_STAR_Alignment

STAR alignment to the reference genome, the gene count matrix, and the checks on it. `./run.sh all` runs every step;
each step also has its own `run.sh` name (below).

## Scripts

### `1_star_index.slurm`
Builds the STAR genome index once per reference and read length (`./run.sh index`), or verifies an existing one.

- **Input:** the genome FASTA and GTF in `REF_DIR` (`./run.sh references`; GENCODE mouse M35 by default)
- **Output:** `REF_DIR/star.overlap<N>.<reference>/` with `index_manifest.json` (reference hashes, STAR version,
  overhang)
- **Runtime:** about 1 hour for mouse or human, 16 cores, 128 GB RAM

---

### `2_star.slurm`
Slurm array job: aligns each sample's trimmed reads (`./run.sh star`).

- **Input:** `01-HTS_Preproc/<sample>/` trimmed reads; `samples.txt`
- **Output:** `02-STAR_alignment/<sample>/`: sorted BAM and `*_ReadsPerGene.out.tab` count tables

---

### `3_star_stats.sh`
Summarizes alignment statistics across all samples (`./run.sh stats`).

- **Output:** `summary_star_alignments.txt` and `strandedness.tsv` (which count column fits the library)

---

### `4_reads2counts.sh`
Builds the gene × sample count matrix from the `ReadsPerGene.out.tab` files (`./run.sh counts`).

- **Output:** `03-Counts/counts.txt`

### `5_confirm_reads2counts.py`, `6_compare_to_published.py`, `7_finalize.slurm`
Check every count cell against its STAR table (`confirm`), compare with an earlier count table (`compare`), and run
the final statistics, counts, checks, `mapping_qc.tsv` and BAM indexes as one job (`finalize`, the last step of
`all`).
