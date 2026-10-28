#!/bin/bash
# Build the STAR genome index. Once per genome, annotation and read length.
#
# The index directory is named for the overhang and the reference (REFERENCE_TAG), so
# indexes for different read lengths can sit side by side in the same REF_DIR.
# If the index is already there the job reports it and exits without rebuilding.
#
#   ./run.sh index
#
# Mouse GRCm39 needs roughly 32 GB for genomeGenerate. This asks for 128 GB
# because an out-of-memory failure here costs an hour. Takes about 50 minutes on
# 16 cores.
#SBATCH --job-name=star_index
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=06:00:00
#SBATCH --output=slurmout/star-index_%A.out
#SBATCH --error=slurmout/star-index_%A.err

set -euo pipefail

# Locate config.sh. SLURM copies the batch script into a spool directory before
# running it, so BASH_SOURCE points at /var/spool/slurmd/... inside a job and
# cannot be used to find anything. run.sh exports PIPELINE_ROOT; SLURM_SUBMIT_DIR
# covers a plain `sbatch` issued from the run directory; BASH_SOURCE covers
# running the script directly with bash.
find_pipeline_root() {
    local c
    for c in "${PIPELINE_ROOT:-}" "${SLURM_SUBMIT_DIR:-}" \
             "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." 2>/dev/null && pwd)" \
             "$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)"; do
        if [[ -n "$c" && -f "$c/config.sh" ]]; then echo "$c"; return 0; fi
    done
    echo "ERROR: cannot locate config.sh. Run via ./run.sh, or export PIPELINE_ROOT=/path/to/run" >&2
    return 1
}
HERE="$(find_pipeline_root)" || exit 1
source "${HERE}/config.sh"
load_software "${MODULE_PYTHON}"

start=$(date +%s)
echo "host: $HOSTNAME"

for f in "${GENOME_FA}" "${GTF}"; do
    [[ -s "$f" ]] || { echo "ERROR: missing $f - run './run.sh references' first" >&2; exit 1; }
done

load_software "${MODULE_STAR}"
if [[ -e "${STAR_INDEX}" ]]; then
    python3 "$HERE/scripts/index_provenance.py" "$STAR_INDEX" "$REF_DIR/reference_manifest.json" "$STAR_SJDB_OVERHANG" "$(STAR --version)" check
    exit 0
fi

mkdir -p "${STAR_INDEX}"

load_software "${MODULE_STAR}"
STAR --version

threads="${SLURM_CPUS_PER_TASK:-8}"
if [[ "${STAR_SJDB_OVERHANG}" == "auto" ]]; then
    echo "ERROR: STAR_SJDB_OVERHANG is still \"auto\". Run ./run.sh check first." >&2
    exit 1
fi
# Small genomes need a smaller suffix-array index (STAR manual: min(14, log2(genome length)/2 - 1)). Mouse and
# human give 14, STAR's default, so their command is unchanged; the option is added only below 14.
genome_bp=$(awk '!/^>/ { n += length($0) } END { print n }' "${GENOME_FA}")
sa_nbases=$(awk -v n="$genome_bp" 'BEGIN { v = int(log(n) / log(2) / 2 - 1); print (v < 14 ? v : 14) }')
extra=()
(( sa_nbases < 14 )) && extra=(--genomeSAindexNbases "$sa_nbases")
echo "threads: ${threads}   sjdbOverhang: ${STAR_SJDB_OVERHANG}   genome: ${genome_bp} bp   SAindexNbases: ${sa_nbases}"
echo "fasta:   ${GENOME_FA}"
echo "gtf:     ${GTF}"
echo "out:     ${STAR_INDEX}"
echo

# --outTmpDir must not exist and must be on a fast local filesystem; STAR
# defaults to ./_STARtmp which lands on the shared filesystem. Left as default
# here because Ceres nodes do not expose a guaranteed local scratch path.
STAR \
    --runThreadN "${threads}" \
    --runMode genomeGenerate \
    --genomeDir "${STAR_INDEX}" \
    --genomeFastaFiles "${GENOME_FA}" \
    --sjdbGTFfile "${GTF}" \
    --sjdbOverhang "${STAR_SJDB_OVERHANG}" \
    "${extra[@]+"${extra[@]}"}"

echo
echo "=== index contents ==="
ls -la "${STAR_INDEX}"
for f in SA SAindex Genome chrName.txt exonInfo.tab sjdbList.out.tab; do
    [[ -s "${STAR_INDEX}/${f}" ]] || { echo "ERROR: index incomplete, missing ${f}" >&2; exit 1; }
done
echo "index looks complete"
echo "annotated junctions: $(wc -l < "${STAR_INDEX}/sjdbList.out.tab")"

echo "runtime: $(( $(date +%s) - start ))s"

python3 "$HERE/scripts/index_provenance.py" "$STAR_INDEX" "$REF_DIR/reference_manifest.json" "$STAR_SJDB_OVERHANG" "$(STAR --version)" record
