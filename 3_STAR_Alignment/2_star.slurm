#!/bin/bash
# STAR alignment, one SLURM array task per sample.
#
# Reads the trimmed FASTQs from HTS_OUT and aligns them against STAR_INDEX.
# Inputs and the index are checked before STAR starts, and the exit status and
# expected outputs are checked afterwards.
#
# Output per sample: <sample>_Aligned.sortedByCoord.out.bam,
# <sample>_ReadsPerGene.out.tab (--quantMode GeneCounts), <sample>_Log.final.out.
#
#   ./run.sh star
#SBATCH --job-name=star
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=slurmout/star_%A_%a.out
#SBATCH --error=slurmout/star_%A_%a.err

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
[[ -n "$sample" ]] || { echo "ERROR: task ${SLURM_ARRAY_TASK_ID} maps to no line in ${SAMPLES_TXT}" >&2; exit 1; }
echo "sample: ${sample}"

if [[ "$LIBRARY_LAYOUT" == paired ]]; then
    reads=("${HTS_OUT}/${sample}/${sample}_R1.fastq.gz" "${HTS_OUT}/${sample}/${sample}_R2.fastq.gz")
else
    reads=("${HTS_OUT}/${sample}/${sample}_SE.fastq.gz")
fi
for f in "${reads[@]}"; do
    [[ -s "$f" ]] || { echo "ERROR: missing HTStream output: $f" >&2; exit 1; }
done
[[ -s "${STAR_INDEX}/SA" ]] || { echo "ERROR: no STAR index at ${STAR_INDEX} - run './run.sh index'" >&2; exit 1; }

outdir="${STAR_OUT}/${sample}"
mkdir -p "${outdir}"

load_software "${MODULE_STAR}"
STAR --version

threads="${SLURM_CPUS_PER_TASK:-8}"
echo "threads: ${threads}"
echo "index:   ${STAR_INDEX}"

STAR \
    --runThreadN "${threads}" \
    --genomeDir "${STAR_INDEX}" \
    --outSAMtype BAM SortedByCoordinate \
    --readFilesCommand zcat \
    --readFilesIn "${reads[@]}" \
    --quantMode GeneCounts \
    --outFileNamePrefix "${outdir}/${sample}_" \
    > "${outdir}/${sample}-STAR.stdout" 2> "${outdir}/${sample}-STAR.stderr"

for f in "${outdir}/${sample}_Aligned.sortedByCoord.out.bam" \
         "${outdir}/${sample}_ReadsPerGene.out.tab"; do
    [[ -s "$f" ]] || {
        echo "ERROR: STAR produced no ${f##*/}" >&2
        tail -20 "${outdir}/${sample}-STAR.stderr" >&2 || true
        exit 1
    }
done

echo
echo "=== alignment summary ==="
grep -E "Number of input reads|Uniquely mapped reads %|% of reads mapped to multiple loci|% of reads unmapped: too short" \
    "${outdir}/${sample}_Log.final.out" || true
echo
echo "=== ReadsPerGene header (STAR's 4 summary rows) ==="
head -4 "${outdir}/${sample}_ReadsPerGene.out.tab"

echo "runtime: $(( $(date +%s) - start ))s"
