#!/bin/bash
# HTStream preprocessing, one SLURM array task per sample.
#
#   hts_Stats          initial read statistics
#   hts_SeqScreener    remove PhiX (built-in reference, no file needed)
#   hts_SeqScreener    count rRNA. The -r flag records the hit and passes every
#                      read through, so nothing is removed at this step.
#   hts_SuperDeduper   remove PCR duplicates (HTS_DEDUP=yes, the default)
#   hts_AdapterTrimmer trim Illumina adapter read-through
#   hts_PolyATTrim     trim poly-A/poly-T tails
#   hts_NTrimmer       trim remaining N bases
#   hts_QWindowTrim    quality trim read ends
#   hts_LengthFilter   discard reads shorter than HTS_MIN_LENGTH (default 50 bp)
#   hts_Stats          final statistics, and writes the trimmed FASTQs
#
# Output per sample: <sample>.json and, for paired runs, <sample>_R1.fastq.gz,
# <sample>_R2.fastq.gz and <sample>_SE.fastq.gz. The alignment step reads R1 and
# R2 only, so the orphaned single-end reads the trimmers rescue are not carried
# forward. For single-end runs the output is <sample>_SE.fastq.gz.
#
#   ./run.sh hts
#SBATCH --job-name=htstream
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=9
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --output=slurmout/htstream_%A_%a.out
#SBATCH --error=slurmout/htstream_%A_%a.err

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

start=$(date +%s)
echo "host: $HOSTNAME"
echo "array task: ${SLURM_ARRAY_TASK_ID:-<none>}"

[[ -s "${SAMPLES_TXT}" ]] || { echo "ERROR: ${SAMPLES_TXT} missing or empty" >&2; exit 1; }
sample=$(sed "${SLURM_ARRAY_TASK_ID}q;d" "${SAMPLES_TXT}")
if [[ -z "$sample" ]]; then
    echo "ERROR: task ${SLURM_ARRAY_TASK_ID} maps to no line in ${SAMPLES_TXT}" >&2
    echo "       ${SAMPLES_TXT} has $(grep -c . "${SAMPLES_TXT}") lines." >&2
    exit 1
fi
echo "sample: ${sample}"

R1="${CONCAT_DIR}/${sample}/${sample}_R1_concatenated.fastq.gz"
R2="${CONCAT_DIR}/${sample}/${sample}_R2_concatenated.fastq.gz"
dedup() {  # HTS_DEDUP=no passes the stream through unchanged
    if [[ "$HTS_DEDUP" == yes ]]; then hts_SuperDeduper -A "${json}" -N 'remove PCR duplicates'; else cat; fi
}
if [[ "$LIBRARY_LAYOUT" == paired ]]; then
    inputs=("$R1" "$R2"); reads=(-1 "$R1" -2 "$R2"); outputs=(R1 R2)
else
    inputs=("$R1"); reads=(-U "$R1"); outputs=(SE)
fi
for f in "${inputs[@]}" "${RRNA_FA}"; do
    [[ -s "$f" ]] || { echo "ERROR: required input missing or empty: $f" >&2; exit 1; }
done

outdir="${HTS_OUT}/${sample}"
mkdir -p "${outdir}"
json="${outdir}/${sample}.json"

load_software "${MODULE_HTSTREAM}"
echo "htstream: $(hts_Stats --version 2>&1 | grep -i version | head -1)"

hts_Stats          -L "${json}" -N 'initial stats' "${reads[@]}" | \
hts_SeqScreener    -A "${json}" -N 'screen phix' | \
hts_SeqScreener    -A "${json}" -N 'count the number of rRNA reads' -r -s "${RRNA_FA}" | \
dedup               | \
hts_AdapterTrimmer -A "${json}" -N 'trim adapters' | \
hts_PolyATTrim     -A "${json}" -N 'remove polyAT' | \
hts_NTrimmer       -A "${json}" -N 'remove any remaining N characters' | \
hts_QWindowTrim    -A "${json}" -N 'quality trim the ends of reads' | \
hts_LengthFilter   -A "${json}" -N "remove reads < ${HTS_MIN_LENGTH}bp" -n -m "${HTS_MIN_LENGTH}" | \
hts_Stats          -A "${json}" -N 'final stats' -f "${outdir}/${sample}"

# The pipeline is one shell pipeline; without pipefail a failure in any stage
# but the last would be invisible. set -o pipefail above catches it.

echo
echo "=== outputs ==="
ls -la "${outdir}"
for o in "${outputs[@]}"; do
    f="${outdir}/${sample}_${o}.fastq.gz"
    [[ -s "$f" ]] || { echo "ERROR: expected output missing: $f" >&2; exit 1; }
    gzip -t "$f" || { echo "ERROR: corrupt output: $f" >&2; exit 1; }
done
echo "${outputs[*]} present and pass gzip -t"

echo "runtime: $(( $(date +%s) - start ))s"
