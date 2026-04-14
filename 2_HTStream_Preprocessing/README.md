# 2_HTStream_Preprocessing

Quality control and adapter trimming using [HTStream](https://ibest.github.io/HTStream/), followed by k-mer correlation QC to verify preprocessing did not distort sample identity.

## Scripts

### `1_hts_preproc.slurm`
SLURM array job that runs the full HTStream QC pipeline on each sample.

**HTStream steps applied (in order):**
1. `hts_Stats` — initial read statistics
2. `hts_SeqScreener` — screen and remove PhiX reads
3. `hts_SeqScreener` — count/remove rRNA reads (mouse rRNA reference required)
4. `hts_SuperDeduper` — remove PCR duplicates
5. `hts_AdapterTrimmer` — trim Illumina adapters
6. `hts_PolyATTrim` — remove poly-A/T tails
7. `hts_NTrimmer` — remove ambiguous N bases
8. `hts_QWindowTrim` — quality window trimming
9. `hts_LengthFilter` — discard reads < 50 bp
10. `hts_Stats` — final read statistics

**Requirements:**
- `samples.txt` — one sample name per line (matches subdirectory names in `concatenated_reads/`)
- `References/mouse_rrna.fasta` — mouse rRNA reference FASTA
- `slurmout/` directory must exist

**Usage:**
```bash
mkdir -p slurmout
# Edit --array=1-N to match number of samples
sbatch 1_hts_preproc.slurm
```

**Output:** `01-HTS_Preproc/<sample>/<sample>.json` (QC stats) and trimmed FASTQ files

---

### `2_HTStream_Kmer_QC.py`
Post-QC validation: computes 7-mer frequency profiles from raw vs. processed reads and calculates Spearman correlations to confirm preprocessing preserved sample identity.

**Logic:** Each sample's raw R1/R2 should correlate highly with its own processed reads. Cross-sample correlations should be lower.

**Requirements:** `biopython`, `scipy`

**Usage:**
```bash
# Edit the samples list in main() before running
python3 2_HTStream_Kmer_QC.py
```

**Output:** `sample_correlations_R1.csv`, `sample_correlations_R2.csv`
