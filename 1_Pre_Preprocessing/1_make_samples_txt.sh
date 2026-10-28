#!/bin/bash
# Build samples.txt from the lane directories.
#
# Runs before concatenation so the sample list exists for the array jobs to
# index. Requires every sample to have an R1 and an R2 in every lane, drops
# names listed in EXCLUDE_SAMPLES, and refuses to write a list it cannot back
# with files.
#
#   ./run.sh samples
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
load_software "${MODULE_PYTHON}"

args=(--regex "${SAMPLE_REGEX}" --r1 "${FASTQ_SUFFIX_R1}" --r2 "${FASTQ_SUFFIX_R2}" --layout "${LIBRARY_LAYOUT}"
      --samples "${SAMPLES_TXT}" --manifest "${WORKDIR}/input_manifest.tsv")
for lane in "${LANE_DIRS[@]}"; do args+=(--lane "$lane"); done
for sample in "${EXCLUDE_SAMPLES[@]}"; do args+=(--exclude "$sample"); done
[[ -z "${EXPECTED_SAMPLES_FILE:-}" ]] || args+=(--expected "$EXPECTED_SAMPLES_FILE")
python3 "${HERE}/scripts/sample_manifest.py" "${args[@]}"
