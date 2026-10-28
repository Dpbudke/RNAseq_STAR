#!/bin/bash
# Download the genome FASTA, the annotation GTF and the rRNA FASTA into REF_DIR.
# Anything already present is left alone, so this is safe to re-run.
#
# Checks each download against the provider's published checksum, records SHA-256,
# then checks that the FASTA and GTF use the same contig naming. A mismatch there
# builds an index with no annotation and a count matrix of zeros, so the job stops.
#
#   ./run.sh references
#SBATCH --job-name=fetch_refs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=03:00:00
#SBATCH --mem=8G
#SBATCH --output=slurmout/fetch_refs_%A.out
#SBATCH --error=slurmout/fetch_refs_%A.err

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

args=(--directory "$REF_DIR" --source "$REFERENCE_SOURCE")
case "$REFERENCE_SOURCE" in
    gencode_mouse) args+=(--release "$GENCODE_RELEASE" --rrna "$HERE/References/mouse_rrna.fasta") ;;
    gencode_human) args+=(--release "$GENCODE_RELEASE") ;;
    ensembl)       args+=(--species "$ENSEMBL_SPECIES" --release "$ENSEMBL_RELEASE") ;;
    custom)        args+=(--name "$CUSTOM_NAME" --genome-fa "$CUSTOM_GENOME_FA" --gtf "$CUSTOM_GTF") ;;
esac
if [[ -n "$RRNA_FASTA" && "$REFERENCE_SOURCE" != gencode_mouse ]]; then args+=(--rrna "$RRNA_FASTA"); fi
python3 "${HERE}/scripts/fetch_references.py" "${args[@]}"
