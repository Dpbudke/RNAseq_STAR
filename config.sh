# Defaults for the workflow. Copy an examples/ profile to run_settings.sh and
# edit that file for each run; its values override these defaults.
# Run ./run.sh scan to inspect naming, then ./run.sh check before submission.


# ======================= EDIT FROM HERE =======================

# Where output goes. Needs room for roughly 3x your raw FASTQ size.
# /90daydata is purged at 90 days by access time. /project is not.
WORKDIR="/90daydata/CHANGEME/rnaseq_run"

# Raw FASTQ directories, one per sequencing lane. Any number, not only two.
# The workflow reads these directories without writing into them.
LANE_DIRS=(
    "/path/to/lane1"
    "/path/to/lane2"
)

# Execution: slurm submits jobs; local runs one stage/sample at a time on Linux.
EXECUTOR="slurm"
SOFTWARE_MODE="modules"  # modules or path (programs already on PATH)
LOCAL_THREADS=8

# SLURM.
#   account:    sacctmgr -nP show assoc user=$USER format=Account | sort -u
#   partition:  scontrol show partition | grep PartitionName
#               Ceres has ceres and scavenger.
SLURM_ACCOUNT="CHANGEME"
SLURM_MAIL="CHANGEME@usda.gov"
SLURM_PARTITION="ceres"

# How a FASTQ filename becomes a sample name. ./run.sh scan shows both.
#   illumina      strip the _S<n>_L<n> block
#                   101_S13_L008_R1_001.fastq.gz          -> 101
#                   STUDY_120_S23_L006_R1_001.fastq.gz    -> STUDY_120
#                   245_ColumnCleaned_S167_L008_R1_001... -> 245_ColumnCleaned
#   first_field   everything before the first underscore
#                   101_S13_L008_R1_001.fastq.gz          -> 101
#                   STUDY_120_S23_L006_R1_001.fastq.gz    -> STUDY   <- collides
#                   245_ColumnCleaned_S167_L008_R1_001... -> 245
#   custom        set SAMPLE_REGEX below, applied with sed -E after the read
#                 suffix is stripped
SAMPLE_NAME_STYLE="illumina"
SAMPLE_REGEX=""

# Samples to leave out of samples.txt.
# Undetermined is the demultiplexer's bin of reads whose index matched no
# sample. It is not a sample, and it is usually one of the largest files.
EXCLUDE_SAMPLES=( "Undetermined" )

# Reference genome and annotation, pinned by release (scripts/fetch_references.py, docs/references.md):
#   gencode_mouse  GENCODE mouse, GENCODE_RELEASE like M35 (GRCm39). The lab's default.
#   gencode_human  GENCODE human, GENCODE_RELEASE like 47 (GRCh38).
#   ensembl        Any Ensembl species: ENSEMBL_SPECIES (the Ensembl folder name, for example
#                  rattus_norvegicus) and ENSEMBL_RELEASE (for example 113).
#   custom         Your own files: CUSTOM_NAME, CUSTOM_GENOME_FA and CUSTOM_GTF (plain or gzip).
# The HTStream rRNA count (a QC figure; it removes nothing) uses the repository's mouse rRNA file for
# gencode_mouse. For other sources it uses RRNA_FASTA when set, otherwise the annotation's rRNA genes.
REFERENCE_SOURCE="gencode_mouse"
GENCODE_RELEASE="M35"
ENSEMBL_SPECIES=""
ENSEMBL_RELEASE=""
CUSTOM_NAME=""
CUSTOM_GENOME_FA=""
CUSTOM_GTF=""
RRNA_FASTA=""

# Where the genome, annotation and STAR index live. Point several runs at the
# same directory and they are downloaded and built once, not per run.
# The index is about 26 GB and takes roughly 50 min to build.
REF_DIR="${WORKDIR}/References"

# Read length minus 1. "auto" reads it from the first FASTQ record.
# Set a number instead to reuse an index that was built with that value.
STAR_SJDB_OVERHANG="auto"

# ======================== EDIT TO HERE ========================
#
# Everything below is derived. Change it only if you have a reason.


default_ref_dir="${WORKDIR}/References"
# Optional run-specific settings, kept beside the template or selected explicitly.
settings="${PIPELINE_CONFIG:-$(dirname "${BASH_SOURCE[0]}")/run_settings.sh}"
if [[ -f "$settings" ]]; then
    export PIPELINE_CONFIG="$(cd "$(dirname "$settings")" && pwd)/$(basename "$settings")"
    source "$PIPELINE_CONFIG"
elif [[ -n "${PIPELINE_CONFIG:-}" ]]; then
    echo "config.sh: configuration file not found: $PIPELINE_CONFIG" >&2
    return 1
fi
case "$EXECUTOR" in slurm|local) ;; *) echo "EXECUTOR must be slurm or local" >&2; return 1;; esac
case "$SOFTWARE_MODE" in modules|path|container) ;; *) echo "SOFTWARE_MODE must be modules, path or container" >&2; return 1;; esac
# SOFTWARE_MODE=container: STAR, SAMtools, HTStream and the QC scripts' Python libraries run in the image STAR_SIF
# (container/Dockerfile, ghcr.io/bjbennettlab/rnaseq_star); the pipeline's own scripts use the host's Python 3.11+
# (MODULE_PYTHON where modules exist). Apptainer runs the image.
STAR_SIF="${STAR_SIF:-}"
MODULE_APPTAINER="${MODULE_APPTAINER:-apptainer}"
container_tools() {
    # One small wrapper per tool in WORKDIR/.container_bin, first on PATH, so every script runs them unchanged.
    local dir="${WORKDIR}/.container_bin" t inner
    [[ -r "$STAR_SIF" ]] || { echo "STAR_SIF is not a readable image: ${STAR_SIF:-unset}" >&2; return 1; }
    mkdir -p "$dir" || return $?
    for t in STAR samtools hts_Stats hts_SeqScreener hts_SuperDeduper hts_AdapterTrimmer hts_PolyATTrim \
             hts_NTrimmer hts_QWindowTrim hts_LengthFilter python; do
        inner=$t; [[ $t == python ]] && inner=python3
        printf '#!/bin/bash\nexec apptainer exec "%s" %s "$@"\n' "$STAR_SIF" "$inner" > "$dir/$t.$$" || return $?
        chmod 0755 "$dir/$t.$$" && mv -f "$dir/$t.$$" "$dir/$t" || return $?
    done
    case ":$PATH:" in *":$dir:"*) ;; *) PATH="$dir:$PATH" ;; esac
}
load_software() {
    if [[ "$SOFTWARE_MODE" == modules ]]; then
        type module >/dev/null 2>&1 || { echo "Environment modules are unavailable; use SOFTWARE_MODE=path for installed programs" >&2; return 1; }
        module load "$@" || return $?
        # Some tool modules include another python3. Keep the selected Python first.
        module unload "$MODULE_PYTHON" || return $?
        module load "$MODULE_PYTHON" || return $?
    elif [[ "$SOFTWARE_MODE" == container ]]; then
        if type module >/dev/null 2>&1; then
            module load "$MODULE_PYTHON" || return $?
            command -v apptainer >/dev/null || module load "$MODULE_APPTAINER" >/dev/null 2>&1 || true
        fi
        container_tools || return $?
    fi
}
# The optional setup installs a checksum-verified STAR build inside this copy.
# Explicit module overrides in the environment or run settings remain available.
runtime="$(dirname "${BASH_SOURCE[0]}")/tools/star_runtime.sh"
if [[ -z "${MODULE_STAR:-}" && -f "$runtime" ]]; then
    source "$runtime" || return $?
fi
if [[ "$REF_DIR" == "$default_ref_dir" ]]; then REF_DIR="${WORKDIR}/References"; fi

# --- FASTQ naming ------------------------------------------------------------
FASTQ_SUFFIX_R1="${FASTQ_SUFFIX_R1:-_R1_001.fastq.gz}"
FASTQ_SUFFIX_R2="${FASTQ_SUFFIX_R2:-_R2_001.fastq.gz}"
# Read layout: paired (an R1 and an R2 file per sample in every lane) or single (R1 files only).
LIBRARY_LAYOUT="${LIBRARY_LAYOUT:-paired}"
case "$LIBRARY_LAYOUT" in
    paired) MATES=(R1 R2) ;;
    single) MATES=(R1) ;;
    *) echo "config.sh: LIBRARY_LAYOUT must be paired or single" >&2; return 1 ;;
esac
# The provider's checksum file in each lane directory (UC Davis Genome Center: @md5Sum.md5). Set SOURCE_MD5=none
# only when the provider gave none; the md5 step then records each FASTQ's SHA-256 instead of checking it.
SOURCE_MD5="${SOURCE_MD5:-@md5Sum.md5}"
# HTStream drops reads shorter than this after trimming. 50 suits 100 to 150 bp reads; use about 30 for 50 bp reads.
HTS_MIN_LENGTH="${HTS_MIN_LENGTH:-50}"
# HTStream duplicate removal (hts_SuperDeduper): yes, the lab's method, or no. Without UMIs, single-end reads that
# start at the same place are often not PCR duplicates (highly expressed genes); check warns for that case.
HTS_DEDUP="${HTS_DEDUP:-yes}"
case "$HTS_DEDUP" in yes|no) ;; *) echo "config.sh: HTS_DEDUP must be yes or no" >&2; return 1 ;; esac
[[ "$HTS_MIN_LENGTH" =~ ^[0-9]+$ ]] || { echo "config.sh: HTS_MIN_LENGTH must be a number" >&2; return 1; }

case "${SAMPLE_NAME_STYLE}" in
    illumina)    SAMPLE_REGEX='s/_S[0-9]+_L[0-9]+$//' ;;
    first_field) SAMPLE_REGEX='s/_.*$//' ;;
    custom)      : ;;   # SAMPLE_REGEX set above
    *) echo "config.sh: SAMPLE_NAME_STYLE must be illumina, first_field or custom" >&2 ;;
esac

# --- stage outputs -----------------------------------------------------------
CONCAT_DIR="${WORKDIR}/concatenated_reads"
HTS_OUT="${WORKDIR}/01-HTS_Preproc"
STAR_OUT="${WORKDIR}/02-STAR_alignment"
COUNTS_DIR="${WORKDIR}/03-Counts"
COUNTS_FILE="${COUNTS_DIR}/counts.txt"
SAMPLES_TXT="${WORKDIR}/samples.txt"
SLURMOUT="${WORKDIR}/slurmout"

# --- references --------------------------------------------------------------
# One reference per REF_DIR. GENCODE mouse keeps the file and index names of the lab's earlier runs, so those
# reference folders and indexes are reused as they are; every other source uses the names below.
case "$REFERENCE_SOURCE" in
    gencode_mouse)
        REFERENCE_TAG="gencode.${GENCODE_RELEASE}"
        GENOME_FA="${REF_DIR}/GRCm39.primary_assembly.genome.fa"
        GTF="${REF_DIR}/gencode.v${GENCODE_RELEASE}.annotation.gtf"
        # mouse_rrna.fasta ships with this repo (UC Davis RNA-seq workshop); the fetch job copies it into REF_DIR.
        RRNA_FA="${REF_DIR}/mouse_rrna.fasta" ;;
    gencode_human) REFERENCE_TAG="gencode_human.${GENCODE_RELEASE}" ;;
    ensembl)       REFERENCE_TAG="ensembl.${ENSEMBL_SPECIES}.${ENSEMBL_RELEASE}" ;;
    custom)        REFERENCE_TAG="custom.${CUSTOM_NAME}" ;;
    *) echo "config.sh: REFERENCE_SOURCE must be gencode_mouse, gencode_human, ensembl or custom" >&2; return 1 ;;
esac
if [[ "$REFERENCE_SOURCE" != gencode_mouse ]]; then
    GENOME_FA="${REF_DIR}/genome.fa"
    GTF="${REF_DIR}/annotation.gtf"
    RRNA_FA="${REF_DIR}/rrna.fasta"
fi

STAR_INDEX="${REF_DIR}/star.overlap${STAR_SJDB_OVERHANG}.${REFERENCE_TAG}"

# --- counts ------------------------------------------------------------------
# Column of STAR's ReadsPerGene.out.tab.
#   2 unstranded, 3 stranded forward, 4 stranded reverse (dUTP kits)
# "auto" picks it from the alignments. 3_star_stats.sh prints all three totals.
COUNT_COLUMN="${COUNT_COLUMN:-auto}"

# --- SLURM -------------------------------------------------------------------
SLURM_MAIL_TYPE="${SLURM_MAIL_TYPE:-NONE}"
# How many array tasks run at once, per stage. Raise on an idle cluster.
THROTTLE_CONCAT="${THROTTLE_CONCAT:-25}"
THROTTLE_HTS="${THROTTLE_HTS:-30}"
THROTTLE_STAR="${THROTTLE_STAR:-20}"

# --- modules -----------------------------------------------------------------
MODULE_HTSTREAM="${MODULE_HTSTREAM:-htstream/1.3.3}"    # binaries inside report 1.3.2
MODULE_STAR="${MODULE_STAR:-star/2.7.11b}"
MODULE_PYTHON="${MODULE_PYTHON:-python_3/3.11.1}"
MODULE_SAMTOOLS="${MODULE_SAMTOOLS:-samtools/1.17}"

# Ceres has no module carrying scipy or pandas, so the two .py scripts need an
# environment. 0_Setup/2_python_env.slurm builds one here.
PY_VENV="${WORKDIR}/.venv"

# Group-writable output so collaborators in the project group can use the run.
umask 007

# Values worked out by ./run.sh check (read length, index path, count column).
# Written into WORKDIR so a run stays reproducible. Overrides the "auto"
# settings above.
# (an `if` here, not `[[ ]] && source`, so sourcing config.sh never returns 1
# when resolved.sh does not exist yet)
if [[ -f "${WORKDIR}/resolved.sh" ]]; then
    source "${WORKDIR}/resolved.sh"
fi
if [[ -f "$WORKDIR/run_records/code_sha256.txt" ]]; then
    sha256sum --strict --status -c "$WORKDIR/run_records/code_sha256.txt" || {
        echo "Pipeline source or run settings changed after the run was recorded" >&2
        return 1
    }
fi
true
