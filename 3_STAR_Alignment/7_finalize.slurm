#!/bin/bash
# Build and validate counts, summaries, and BAM indexes so the pipeline can be
# chained with --dependency and finish without anyone watching it.
#
# Both are cheap enough to run by hand instead:
#   ./run.sh stats
#   ./run.sh counts

#SBATCH --job-name=finalize
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=08:00:00
#SBATCH --mem=8G
#SBATCH --output=slurmout/finalize_%A.out
#SBATCH --error=slurmout/finalize_%A.err

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

echo "############ 3_star_stats.sh ############"
bash "${HERE}/3_STAR_Alignment/3_star_stats.sh"

echo
echo "############ 4_reads2counts.sh ############"
bash "${HERE}/3_STAR_Alignment/4_reads2counts.sh"

source "${HERE}/config.sh"
python3 "$HERE/3_STAR_Alignment/5_confirm_reads2counts.py" \
    --star-dir "$STAR_OUT" --counts "$COUNTS_FILE" --samples "$SAMPLES_TXT" \
    --column "$COUNT_COLUMN" --report "$WORKDIR/count_validation.json"
python3 "$HERE/2_HTStream_Preprocessing/3_hts_summary.py" \
    --hts-dir "$HTS_OUT" --samples "$SAMPLES_TXT" --out "$WORKDIR/summary_htstream.txt"
load_software "${MODULE_SAMTOOLS:-samtools/1.17}"
python3 "$HERE/scripts/final_qc.py" "$WORKDIR"
printf 'validated\n' > "$WORKDIR/pipeline_state.txt"
