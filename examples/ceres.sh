# Copy to run_settings.sh and edit paths/account before running.
EXECUTOR="slurm"
SOFTWARE_MODE="modules"
WORKDIR="/90daydata/YOUR_PROJECT/rnaseq_run"
LANE_DIRS=(
    "/project/YOUR_PROJECT/study/raw/lane1"
    "/project/YOUR_PROJECT/study/raw/lane2"
)
SLURM_ACCOUNT="YOUR_PROJECT"
SLURM_PARTITION="ceres"
SLURM_MAIL=""
SAMPLE_NAME_STYLE="illumina"
LIBRARY_LAYOUT="paired"            # or single
SOURCE_MD5="@md5Sum.md5"           # the provider's checksum file in each lane, or none
REFERENCE_SOURCE="gencode_mouse"   # gencode_human, ensembl or custom (config.sh lists their settings)
GENCODE_RELEASE="M35"              # gencode_human: for example 47
# ENSEMBL_SPECIES="rattus_norvegicus"; ENSEMBL_RELEASE="113"   # for REFERENCE_SOURCE="ensembl"
STAR_SJDB_OVERHANG="auto"
COUNT_COLUMN="auto"
# Optional: run the tools from the lab's image instead of modules (docs/runtime.md).
# SOFTWARE_MODE="container"
# STAR_SIF="/project/YOUR_PROJECT/containers/rnaseq_star_1.0.0.sif"
# Optional: put references in a shared, writable reference directory.
# REF_DIR="/90daydata/YOUR_PROJECT/references/gencode_m35"
# STAR_CPUS=8
# STAR_MEMORY="48G"
# INDEX_CPUS=16
# INDEX_MEMORY="128G"
