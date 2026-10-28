# RNAseq_STAR

By Dawson Budke ([@Dpbudke](https://github.com/Dpbudke)), Bennett Lab.

RNA-seq from paired-end or single-end reads, for mouse, human, any Ensembl species
or your own genome: combine sequencing lanes, preprocess with HTStream, align with
STAR, and produce checked gene counts and BAM files.

## Pipeline overview

```
Raw FASTQs (one directory per sequencing lane)
    │
    ▼
0_Setup                  — references, Python environment, provider checksums, STAR build
    │
    ▼
1_Pre_Preprocessing      — list samples by exact file name, combine each sample's lanes
    │
    ▼
2_HTStream_Preprocessing — PhiX screen, rRNA count, duplicates, adapters, quality trimming
    │
    ▼
3_STAR_Alignment         — STAR index, alignment, count matrix and checks
    │
    ▼
4_DE_Analysis            — differential expression notebooks (R Markdown, run separately)
```

Each numbered folder has a README for its scripts. `run.sh` runs the stages in order.

## Setup and run

On **Ceres (SCINet)**, use **Clusters > Ceres Shell Access** in
[Open OnDemand](https://scinet.usda.gov/guides/use/open-ondemand), or connect by SSH.
Work in a directory you can write to under your project allocation. Raw data
belongs in persistent project storage; processing outputs can go under `/90daydata`.

Clone the repository:

```bash
git clone https://github.com/BJBennettlab/RNAseq_STAR.git
cd RNAseq_STAR
```

Keep raw FASTQs in separate lane directories. Each lane must contain the same
samples and the provider's checksum file (`@md5Sum.md5` by default; set `SOURCE_MD5`
to its name, or to `none` when the provider gave none). If a study was sequenced as
separate pools on different lanes, run each pool as its own run and join the
validated count tables afterwards with `scripts/combine_runs.py` (see
[SETUP_INSTRUCTIONS.txt](SETUP_INSTRUCTIONS.txt)). Continue on Ceres:

```bash
cp examples/ceres.sh run_settings.sh
nano run_settings.sh
python3 0_Setup/4_fetch_star.py
./run.sh scan
./run.sh check
./run.sh all --dry-run
./run.sh all
./run.sh status
```

In `run_settings.sh`, set the output directory, the raw lane paths, your SCINet
project as `SLURM_ACCOUNT`, the read layout (`LIBRARY_LAYOUT`, paired or single) and
the reference (`REFERENCE_SOURCE`, below). `check` refuses an account you do not
have, and work or reference folders outside that project's `/project` or `/90daydata`.
`scan` shows filenames and sample naming; `check` checks settings and software.
`all` submits dependent Slurm jobs. The numbered scripts still perform each stage.
The STAR setup downloads a fixed, checksum-verified Linux build into this copy.

`SOFTWARE_MODE=container` runs STAR, SAMtools and HTStream from the lab's pinned
image (`ghcr.io/bjbennettlab/rnaseq_star`) instead of Ceres modules; it gives the same
results as the modules ([runtime setup](docs/runtime.md#container-image-software_modecontainer)).

On **another Linux system**, use `examples/linux.sh`. Install or activate Python
3.11+, HTStream, STAR, SAMtools and GNU command-line tools first. Local mode runs
one sample at a time without Slurm; keep the process running until it finishes.
A workstation needs enough RAM for the mouse index and BAM sorting, typically at
least 48 GiB. Another Slurm cluster can use either site modules or programs on PATH.

[SETUP_INSTRUCTIONS.txt](SETUP_INSTRUCTIONS.txt) covers each setting, stage-by-stage
commands, checks and recovery. [Runtime setup](docs/runtime.md) describes software
installation and recorded builds. A new study needs a new output directory.

## Results

| Run directory | Contents |
| --- | --- |
| `03-Counts/counts.txt` | Integer gene-by-sample count matrix |
| `02-STAR_alignment/` | Sorted BAMs, BAM indexes, STAR count tables and logs |
| `summary_htstream.txt`, `mapping_qc.tsv`, `strandedness.tsv` | Per-sample QC |
| `count_validation.json`, `bam_validation.json` | Count comparison and BAM checks |
| `samples.txt`, `input_manifest.tsv` | Sample order and exact inputs from every lane |
| `run_records/` | Source snapshot, configuration and executable identities |
| `jobs.tsv` or `local_jobs.json`, `slurmout/` | Task records and output/error logs |

Each FASTQ is assigned to a sample by its full name (`<sample>_S<n>_L<lane>_R1_001.fastq.gz`
by default), so a name that begins another name, such as `MOUSE_1` and `MOUSE_10`,
always keeps its own files. `input_manifest.tsv` records the exact file used for
every sample and lane.

Processing does not modify raw files. Missing mates, different sample sets, name collisions,
failed source checks and failed pipeline commands stop processing. Finalization
compares every count cell with its STAR table, checks BAMs and builds indexes.
Require successful task completion and passed validation reports. An empty error
log alone is not a completion check. [Validation commands](docs/validation.md)
show how to repeat the checks.

## References and analysis choices

`REFERENCE_SOURCE` picks the genome and annotation, pinned by release:

| Settings | Reference |
| --- | --- |
| `gencode_mouse`, `GENCODE_RELEASE="M35"` (default) | GRCm39, GENCODE M35 |
| `gencode_human`, `GENCODE_RELEASE="47"` | GRCh38, GENCODE 47 |
| `ensembl`, `ENSEMBL_SPECIES="rattus_norvegicus"`, `ENSEMBL_RELEASE="113"` | any Ensembl species and release |
| `custom`, `CUSTOM_NAME`, `CUSTOM_GENOME_FA`, `CUSTOM_GTF` | your own FASTA and GTF |

The reference stage downloads the files, checks the provider's published checksums,
records SHA-256 in `REF_DIR/reference_manifest.json`, and checks that the genome and
annotation name chromosomes the same way. One `REF_DIR` holds one reference; point
several runs at it and it is downloaded and indexed once. Keep the reference fixed
within a study: providers keep old releases, so no update is needed, and a different
release needs a new `REF_DIR`, index and run. HTStream's rRNA count is a QC figure:
it uses `References/mouse_rrna.fasta` for GENCODE mouse, `RRNA_FASTA` when given, and
otherwise the annotation's rRNA genes. See [reference inputs](docs/references.md) for
details and the separate inputs needed by the original R analyses.

Sample names retain underscores by default. Any number of lanes is supported.
`COUNT_COLUMN=auto` checks strand estimates per sample and stops if they disagree;
confirm the library preparation when interpreting this choice. HTStream's original
duplicate removal and trimming are retained (`HTS_DEDUP=no` turns duplicate removal off; `check` warns about it for
single-end reads, where many reads starting at the same place are not PCR duplicates). rRNA is measured, not removed.
In paired-end runs, trimmed singletons are saved but not aligned. `HTS_MIN_LENGTH`
(default 50) is the shortest read kept after trimming; use about 30 for 50 bp reads.

`4_DE_Analysis/` contains the original study-specific R analyses. They are not run
by `all`; they need their own metadata, design and contrasts. GBRS and QTL mapping
are separate analyses. [CHANGELOG.md](CHANGELOG.md) records pipeline changes.

## Study context

The scripts were first written for RNA-seq of CC and CC-RIX pup samples (about 246 samples) studying diet effects on
gene expression, on GRCm39. The lab's Diversity Outbred studies were later processed with the same workflow, which
now also takes other organisms and single-end libraries.

## Related repository

GBRS-based allele-specific expression analysis is in [RNAseq_GBRS](https://github.com/BJBennettlab/RNAseq_GBRS).

## Author

Dawson Budke — dpbudke@ucdavis.edu. His original repository: [Dpbudke/RNAseq_STAR](https://github.com/Dpbudke/RNAseq_STAR).
