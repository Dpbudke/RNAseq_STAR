#!/bin/bash
# Check the raw lane FASTQs against the provider's checksum file (SOURCE_MD5,
# default @md5Sum.md5), one array task per lane. The file must cover every FASTQ.
# With SOURCE_MD5=none the job records each FASTQ's SHA-256 instead.
#
#   ./run.sh md5
#
# Read only. Nothing in this job writes into LANE_DIRS.
# Allow about half an hour per 500 GB.
#SBATCH --job-name=md5check
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=08:00:00
#SBATCH --mem=4G
#SBATCH --output=slurmout/md5_%A_%a.out
#SBATCH --error=slurmout/md5_%A_%a.err

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
idx="${SLURM_ARRAY_TASK_ID:-1}"
lane="${LANE_DIRS[$((idx - 1))]}"
echo "host: $HOSTNAME"
echo "lane: ${lane}"

outdir="${WORKDIR}/md5_check"
mkdir -p "${outdir}"
report="${outdir}/lane$(basename "${lane}").txt"

if [[ "$SOURCE_MD5" == none ]]; then
    # No provider checksums to check against: record each FASTQ's SHA-256 so later copies can be compared.
    echo "SOURCE_MD5=none: no provider checksums; recording SHA-256 of every FASTQ"
    (cd "${lane}" && find . -maxdepth 1 -name '*.fastq.gz' -printf '%f\n' | sort | xargs -r sha256sum) > "${report%.txt}.sha256"
    echo "files: $(grep -c . "${report%.txt}.sha256")  report: ${report%.txt}.sha256"
    exit 0
fi
manifest="${lane}/${SOURCE_MD5}"
[[ -s "$manifest" ]] || { echo "no checksum file at ${manifest} (set SOURCE_MD5 to its name, or none)" >&2; exit 1; }

echo "manifest entries: $(grep -c . "$manifest")"
echo "fastq.gz present: $(find "${lane}" -maxdepth 1 -name '*.fastq.gz' | wc -l)"
echo

# md5sum -c wants to run from the directory the manifest's relative paths are
# rooted at. cd into the lane, but write the report outside it.
cd "${lane}"
python3 - "$manifest" <<'PYMD5'
import pathlib,re,sys
names=[]
for line in pathlib.Path(sys.argv[1]).read_text().splitlines():
    m=re.fullmatch(r'[0-9a-fA-F]{32} [ *](?:\./)?([^/]+)',line)
    if not m:raise SystemExit('Invalid source checksum entry')
    names.append(m[1])
actual={p.name for p in pathlib.Path('.').glob('*.fastq.gz')}
if len(names)!=len(set(names)) or set(names)!=actual:raise SystemExit('Source MD5 does not cover every FASTQ exactly once')
PYMD5
status=0
md5sum --strict -c "$manifest" > "${report}" 2>&1 || status=$?

ok=$(grep -c ': OK$'     "${report}" || true)
bad=$(grep -c ': FAILED$' "${report}" || true)
missing=$(grep -c 'No such file' "${report}" || true)

echo "OK:      ${ok}"
echo "FAILED:  ${bad}"
echo "missing: ${missing}"
echo "report:  ${report}"

if [[ "${bad}" -gt 0 || "${missing}" -gt 0 ]]; then
    echo
    echo "=== problems ==="
    grep -E ': FAILED$|No such file' "${report}" | head -40
fi

echo
echo "runtime: $(( $(date +%s) - start ))s"
exit "${status}"
