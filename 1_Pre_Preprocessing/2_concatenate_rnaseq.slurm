#!/bin/bash
# Join one sample's reads across every lane. One SLURM array task per sample.
#
# Writes straight into CONCAT_DIR/<sample>/, so no separate organize step is
# needed. A sample missing from any lane fails that task.
#
# Verifies its own output: the concatenated file must decompress cleanly, the
# line count must be a multiple of 4, and R1 and R2 must hold the same number of
# records.
#
# cat on .gz files is correct here. The gzip format allows members to be joined
# end to end, and zcat reads them back as one stream.
#
#   ./run.sh concat
#
# Read only on LANE_DIRS.
#SBATCH --job-name=concat
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=06:00:00
#SBATCH --mem=4G
#SBATCH --output=slurmout/concat_%A_%a.out
#SBATCH --error=slurmout/concat_%A_%a.err

set -euo pipefail

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
load_software "$MODULE_PYTHON"

start=$(date +%s)
echo "host: $HOSTNAME"
echo "array task: ${SLURM_ARRAY_TASK_ID:-<none>}"

[[ -s "${SAMPLES_TXT}" ]] || { echo "ERROR: ${SAMPLES_TXT} missing or empty - run './run.sh samples' first" >&2; exit 1; }
sample=$(sed "${SLURM_ARRAY_TASK_ID}q;d" "${SAMPLES_TXT}")
[[ -n "$sample" ]] || { echo "ERROR: task ${SLURM_ARRAY_TASK_ID} maps to no line in ${SAMPLES_TXT}" >&2; exit 1; }
echo "sample: ${sample}"

outdir="${CONCAT_DIR}/${sample}"
mkdir -p "${outdir}"

for rd in "${MATES[@]}"; do
    selected=$(python3 "$HERE/scripts/manifest_inputs.py" "$WORKDIR/input_manifest.tsv" "$sample" "$rd" "${LANE_DIRS[@]}")
    mapfile -t files <<< "$selected"

    out="${outdir}/${sample}_${rd}_concatenated.fastq.gz"
    echo "${rd} <- ${#files[@]} lane(s):"
    printf '    %s\n' "${files[@]}"
    cat "${files[@]}" > "${out}"

    # One decompression pass does both jobs: zcat fails on a corrupt or
    # truncated member (set -o pipefail propagates that), and the line count
    # gives the record count. A truncated source file therefore fails here,
    # which is why a separate gzip -t is not also run - it would be a second
    # full read of the same data for no extra coverage.
    lines=$(zcat "${out}" | wc -l)
    if (( lines % 4 != 0 )); then
        echo "ERROR: ${out} has ${lines} lines, not a multiple of 4" >&2
        exit 1
    fi
    got=$(( lines / 4 ))
    echo "    -> ${out}  ${got} records, decompresses clean"
    declare "count_${rd}=${got}"
done

if [[ "$LIBRARY_LAYOUT" == paired ]]; then
    if [[ "${count_R1}" -ne "${count_R2}" ]]; then
        echo "ERROR: ${sample} R1=${count_R1} R2=${count_R2} - reads are not paired" >&2
        exit 1
    fi
    echo
    echo "${sample}: ${count_R1} read pairs across ${#LANE_DIRS[@]} lanes"
else
    echo
    echo "${sample}: ${count_R1} single-end reads across ${#LANE_DIRS[@]} lanes"
fi
ls -la "${outdir}"
echo "runtime: $(( $(date +%s) - start ))s"
