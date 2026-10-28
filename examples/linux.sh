# Copy to run_settings.sh. Activate the required tools before running.
EXECUTOR="local"
SOFTWARE_MODE="path"
LOCAL_THREADS=8
WORKDIR="/data/YOUR_STUDY/rnaseq_run"
LANE_DIRS=(
    "/data/YOUR_STUDY/raw/lane1"
    "/data/YOUR_STUDY/raw/lane2"
)
SAMPLE_NAME_STYLE="illumina"
GENCODE_RELEASE="M35"
STAR_SJDB_OVERHANG="auto"
COUNT_COLUMN="auto"
# Add/remove lane entries as needed. Each lane must have the same samples.
# For another Slurm cluster, set EXECUTOR="slurm", then set these:
# SLURM_ACCOUNT="your_account"   # empty uses the site's default
# SLURM_PARTITION="your_partition"  # empty uses the site's default
# SLURM_MAIL=""
# Programs may be supplied by the site's modules or an activated environment.
