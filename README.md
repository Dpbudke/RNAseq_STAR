# RNAseq_STAR

RNA-seq preprocessing and differential expression analysis pipeline for Collaborative Cross (CC) mouse studies. This repository covers the full workflow from raw FASTQ concatenation through STAR alignment and downstream DE analysis.

## Pipeline Overview

```
Raw FASTQs (multi-lane)
    │
    ▼
1_Pre_Preprocessing     — Concatenate lanes, organize by sample
    │
    ▼
2_HTStream_Preprocessing — Quality trimming, adapter removal, deduplication
    │
    ▼
3_STAR_Alignment        — Genome index build, read alignment, read-to-count matrix
    │
    ▼
4_DE_Analysis           — Differential expression (R/RMarkdown)
```

## Dependencies

| Stage | Tools |
|-------|-------|
| Pre-processing | bash |
| QC/Trimming | [HTStream v1.3.3](https://ibest.github.io/HTStream/) |
| Alignment | [STAR](https://github.com/alexdobin/STAR) |
| DE Analysis | R, DESeq2 / edgeR (see individual Rmd files) |

All HPC scripts are written for SLURM schedulers and assume module-based software loading (`module load htstream`, `module load star`).

## Repository Structure

```
RNAseq_STAR/
├── 1_Pre_Preprocessing/       # Lane concatenation and sample organization
├── 2_HTStream_Preprocessing/  # HTStream QC and k-mer validation
├── 3_STAR_Alignment/          # STAR index, alignment, and count matrix generation
└── 4_DE_Analysis/             # Differential expression R Markdown scripts
```

## Quick Start

1. **Concatenate lanes** — edit paths in `1_Pre_Preprocessing/1_concatenate_rnaseq.sh` and run
2. **Organize samples** — run `1_Pre_Preprocessing/2_organize_samples.sh` in the concatenated reads directory
3. **HTStream QC** — create `samples.txt` (one sample per line), then `sbatch 2_HTStream_Preprocessing/1_hts_preproc.slurm`
4. **STAR index** — `sbatch 3_STAR_Alignment/1_star_index.slurm` (one-time, ~1-2h)
5. **STAR alignment** — `sbatch 3_STAR_Alignment/2_star.slurm`
6. **Build count matrix** — `bash 3_STAR_Alignment/4_reads2counts.sh`
7. **DE analysis** — open and knit the relevant Rmd in `4_DE_Analysis/`

## Study Context

Scripts were developed for RNA-seq analysis of CC pup samples (~246 samples) studying diet effects on gene expression. Reference genome: combined CC multi-strain FASTA + annotation GTF.

## Related Repository

GBRS-based allele-specific expression analysis has been separated into [RNAseq_GBRS](https://github.com/Dpbudke/RNAseq_GBRS).

## Author

Dawson Budke — dpbudke@ucdavis.edu
