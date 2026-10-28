#!/bin/bash
# Summarise STAR alignment statistics into one table, and work out which
# ReadsPerGene column matches the library.
#
# STAR --quantMode GeneCounts writes four columns:
#   col1 gene id
#   col2 counts, library treated as unstranded
#   col3 counts, read 1 on the same strand as the RNA  (stranded forward)
#   col4 counts, read 1 on the opposite strand         (stranded reverse)
#
# With COUNT_COLUMN set to auto, this sums all three across every sample, picks
# the one that fits, and writes it into resolved.sh so the count matrix uses the
# same value. The totals are printed so the choice can be checked.
#
#   ./run.sh stats
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

[[ -s "${SAMPLES_TXT}" ]] || { echo "ERROR: ${SAMPLES_TXT} missing or empty" >&2; exit 1; }
[[ -d "${STAR_OUT}" ]]    || { echo "ERROR: ${STAR_OUT} not found" >&2; exit 1; }

tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
out="${WORKDIR}/summary_star_alignments.txt"

missing=0
while read -r s; do
    [[ -s "${STAR_OUT}/${s}/${s}_ReadsPerGene.out.tab" ]] || {
        echo "MISSING: ${STAR_OUT}/${s}/${s}_ReadsPerGene.out.tab" >&2; missing=1; }
done < "${SAMPLES_TXT}"
[[ $missing -eq 0 ]] || { echo "aborting" >&2; exit 1; }

# Resolve COUNT_COLUMN. "auto" measures it from the alignments and writes the
# answer into resolved.sh so 4_reads2counts.sh uses the same value.
col=$(python3 "$HERE/scripts/strand_report.py" "$STAR_OUT" "$SAMPLES_TXT" "$WORKDIR/strandedness.tsv" "$COUNT_COLUMN")
echo "COUNT_COLUMN=${col}" >> "${WORKDIR}/resolved.sh"
echo "Selected STAR column ${col}; per-sample evidence: ${WORKDIR}/strandedness.tsv"

# Row labels: STAR's four summary rows plus the in-feature total.
first=$(head -1 "${SAMPLES_TXT}")
{ head -4 "${STAR_OUT}/${first}/${first}_ReadsPerGene.out.tab" | cut -f1
  echo "total_in_feature"; } > "${tmp}/labels"

printf 'sample_name' > "${tmp}/header"
while read -r s; do
    printf '\t%s' "$s" >> "${tmp}/header"
    f="${STAR_OUT}/${s}/${s}_ReadsPerGene.out.tab"
    { head -4 "$f" | cut -f"${col}"
      awk -v c="${col}" 'NR>4 {t+=$c} END {print t+0}' "$f"; } > "${tmp}/col_${s}"
done < "${SAMPLES_TXT}"
printf '\n' >> "${tmp}/header"

paste_args=("${tmp}/labels")
while read -r s; do paste_args+=("${tmp}/col_${s}"); done < "${SAMPLES_TXT}"
{ cat "${tmp}/header"; paste "${paste_args[@]}"; } > "${out}"

echo "wrote ${out}  (using ReadsPerGene column ${col})"
echo
column -t -s $'\t' "${out}" | cut -c1-200

echo
echo "=== strandedness check, summed over all samples ==="
u=0; fwd=0; rev=0
while read -r s; do
    f="${STAR_OUT}/${s}/${s}_ReadsPerGene.out.tab"
    read -r a b c < <(awk 'NR>4 {u+=$2; f+=$3; r+=$4} END {print u+0, f+0, r+0}' "$f")
    u=$((u+a)); fwd=$((fwd+b)); rev=$((rev+c))
done < "${SAMPLES_TXT}"
printf '  col2 unstranded        : %d\n' "$u"
printf '  col3 stranded forward  : %d\n' "$fwd"
printf '  col4 stranded reverse  : %d\n' "$rev"
if   [[ $rev -gt $((fwd * 3)) ]]; then verdict="col4 (reverse-stranded, e.g. TruSeq Stranded mRNA) - matches the default"
elif [[ $fwd -gt $((rev * 3)) ]]; then verdict="col3 (forward-stranded) - set COUNT_COLUMN=3 in config.sh"
else                                   verdict="col2 (unstranded) - set COUNT_COLUMN=2 in config.sh"
fi
echo "  suggested column: ${verdict}"
