# 2_HTStream_Preprocessing

Quality control and trimming with [HTStream](https://ibest.github.io/HTStream/), followed by a k-mer check that
preprocessing kept each sample's identity. `./run.sh hts` (or `all`) runs the preprocessing; `htsqc` and `kmerqc`
summarize it.

## Scripts

### `1_hts_preproc.slurm`
Slurm array job that runs the HTStream steps on each sample, paired-end (`-1`/`-2`) or single-end (`-U`).

**HTStream steps applied (in order):**
1. `hts_Stats` — initial read statistics
2. `hts_SeqScreener` — screen and remove PhiX reads
3. `hts_SeqScreener` — count rRNA reads (counted only; nothing is removed)
4. `hts_SuperDeduper` — remove PCR duplicates (`HTS_DEDUP=yes`, the default)
5. `hts_AdapterTrimmer` — trim Illumina adapters
6. `hts_PolyATTrim` — remove poly-A/T tails
7. `hts_NTrimmer` — remove ambiguous N bases
8. `hts_QWindowTrim` — quality window trimming
9. `hts_LengthFilter` — discard reads shorter than `HTS_MIN_LENGTH` (default 50 bp)
10. `hts_Stats` — final read statistics

- **Requirements:** `samples.txt`; the rRNA FASTA prepared by the reference stage (`mouse_rrna.fasta` for GENCODE
  mouse, `rrna.fasta` otherwise)
- **Output:** `01-HTS_Preproc/<sample>/<sample>.json` (QC stats) and the trimmed reads (`_R1`, `_R2` and the
  orphaned `_SE` for paired-end runs; `_SE` for single-end runs)

---

### `2_HTStream_Kmer_QC.py`
Post-QC validation: 7-mer frequency profiles of raw and processed reads, with Spearman correlations to confirm each
sample still matches itself best. Run with `./run.sh kmerqc`.

- **Logic:** each sample's raw reads should correlate most with its own processed reads; cross-sample correlations
  should be lower
- **Output:** `sample_correlations_R1.csv` (and `_R2.csv` for paired-end runs)

---

### `3_hts_summary.py`
One table of the HTStream statistics of every sample (`./run.sh htsqc`).

- **Output:** `summary_htstream.txt`: reads in, PhiX, rRNA, duplicates, adapter bases trimmed, reads and bases kept
