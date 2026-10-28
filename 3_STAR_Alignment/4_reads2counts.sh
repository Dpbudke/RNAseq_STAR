#!/bin/bash
# Build the gene x sample count matrix from the per-sample ReadsPerGene tables,
# then verify it against them.
#
# Sample order comes from samples.txt. Every sample's gene list is compared
# against the first sample's before anything is pasted, and paste runs in chunks
# so the open file limit is not a constraint.
#
# Verification: column sums in the matrix against the source tables for every
# sample, plus a full column comparison for the first and last sample.
#
#   ./run.sh counts
#
# Needs COUNT_COLUMN resolved, which 3_star_stats.sh does.
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
col="${COUNT_COLUMN:-4}"
if [[ "$col" == "auto" ]]; then
    echo "ERROR: COUNT_COLUMN is still \"auto\". Run 3_star_stats.sh first, it" >&2
    echo "       measures the strandedness and writes the answer to resolved.sh." >&2
    exit 1
fi
mkdir -p "${COUNTS_DIR}"
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT

echo "count column: ${col}"
echo

# --- gene ids from the first sample, then require every sample to match -------
first=$(head -1 "${SAMPLES_TXT}")
tail -n +5 "${STAR_OUT}/${first}/${first}_ReadsPerGene.out.tab" | cut -f1 > "${tmp}/geneids"
ngenes=$(wc -l < "${tmp}/geneids")
echo "genes (from ${first}): ${ngenes}"

fail=0
while read -r s; do
    f="${STAR_OUT}/${s}/${s}_ReadsPerGene.out.tab"
    [[ -s "$f" ]] || { echo "MISSING: $f" >&2; fail=1; continue; }
    tail -n +5 "$f" | cut -f1 > "${tmp}/g_${s}"
    if ! cmp -s "${tmp}/geneids" "${tmp}/g_${s}"; then
        echo "GENE LIST MISMATCH: ${s} does not match ${first}" >&2
        echo "  $(wc -l < "${tmp}/g_${s}") genes vs ${ngenes}" >&2
        fail=1
    fi
    tail -n +5 "$f" | cut -f"${col}" > "${tmp}/c_${s}"
done < "${SAMPLES_TXT}"
[[ $fail -eq 0 ]] || { echo "aborting, no matrix written" >&2; exit 1; }
echo "all samples share the same gene list"

# --- paste in chunks ---------------------------------------------------------
cp "${tmp}/geneids" "${tmp}/acc"
chunk=()
flush() {
    [[ ${#chunk[@]} -eq 0 ]] && return
    paste "${tmp}/acc" "${chunk[@]}" > "${tmp}/acc.new"
    mv "${tmp}/acc.new" "${tmp}/acc"
    chunk=()
}
while read -r s; do
    chunk+=("${tmp}/c_${s}")
    [[ ${#chunk[@]} -ge 200 ]] && flush
done < "${SAMPLES_TXT}"
flush

{
    printf 'gene_id'
    while read -r s; do printf '\t%s' "$s"; done < "${SAMPLES_TXT}"
    printf '\n'
    cat "${tmp}/acc"
} > "${COUNTS_FILE}"

echo
echo "wrote ${COUNTS_FILE}"
echo "  $(wc -l < "${COUNTS_FILE}") lines (1 header + ${ngenes} genes)"
echo "  $(head -1 "${COUNTS_FILE}" | awk -F'\t' '{print NF-1}') sample columns"

# --- verify against the source tables ----------------------------------------
echo
echo "=== verification: column sums in the matrix vs the STAR tables ==="
vfail=0
i=1
while read -r s; do
    i=$((i+1))
    a=$(awk -F'\t' -v c="$i" 'NR>1 {t+=$c} END {print t+0}' "${COUNTS_FILE}")
    b=$(awk -F'\t' -v c="$col" 'NR>4 {t+=$c} END {print t+0}' "${STAR_OUT}/${s}/${s}_ReadsPerGene.out.tab")
    if [[ "$a" == "$b" ]]; then
        printf '  %-20s %12s  ok\n' "$s" "$a"
    else
        printf '  %-20s matrix=%s star=%s  MISMATCH\n' "$s" "$a" "$b"
        vfail=1
    fi
done < "${SAMPLES_TXT}"

# A sum check alone would miss two genes swapping values, so also compare the
# full column for the first and last sample byte for byte.
echo
echo "=== verification: full column comparison, first and last sample ==="
last=$(tail -1 "${SAMPLES_TXT}")
for s in "$first" "$last"; do
    idx=$(awk -F'\t' -v n="$s" 'NR==1 {for(j=2;j<=NF;j++) if($j==n) {print j; exit}}' "${COUNTS_FILE}")
    awk -F'\t' -v c="$idx" 'NR>1 {print $c}' "${COUNTS_FILE}" > "${tmp}/m_${s}"
    if cmp -s "${tmp}/m_${s}" "${tmp}/c_${s}"; then
        echo "  ${s}: identical"
    else
        echo "  ${s}: DIFFERS"; vfail=1
    fi
done

echo
[[ $vfail -eq 0 ]] && echo "count matrix verified" || { echo "VERIFICATION FAILED" >&2; exit 1; }
