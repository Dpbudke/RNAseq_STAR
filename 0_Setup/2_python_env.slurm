#!/bin/bash
# Build the Python environment for the two .py scripts.
#
#   2_HTStream_Preprocessing/2_HTStream_Kmer_QC.py  needs biopython, scipy
#   3_STAR_Alignment/5_confirm_reads2counts.py      needs pandas, numpy
#
# No module on Ceres carries scipy or pandas, so the venv is built here rather
# than loaded. It goes under WORKDIR, not $HOME, which has a 30 GB quota.
#
#   ./run.sh pyenv

#SBATCH --job-name=py_env
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --output=slurmout/py_env_%A.out
#SBATCH --error=slurmout/py_env_%A.err

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
echo "python: $(python3 --version)"

if [[ "$SOFTWARE_MODE" == container ]]; then
    # The image holds requirements.txt's packages; PY_VENV/bin/python runs its Python.
    mkdir -p "${PY_VENV}/bin"
    ln -sfn "${WORKDIR}/.container_bin/python" "${PY_VENV}/bin/python"
    "${PY_VENV}/bin/python" - > "${WORKDIR}/python_packages.txt" <<'PYLIST'
import importlib.metadata as m
print("\n".join(sorted(f"{d.metadata['Name']}=={d.version}" for d in m.distributions())))
PYLIST
else
    if [[ -x "${PY_VENV}/bin/python" ]]; then
        echo "venv already at ${PY_VENV}"
    else
        python3 -m venv "${PY_VENV}"
    fi

    # pip caches into $HOME by default; keep it inside the run directory.
    export PIP_CACHE_DIR="${PY_VENV}/.pip-cache"

    "${PY_VENV}/bin/python" -m pip install -r "${HERE}/requirements.txt"
    "${PY_VENV}/bin/python" -m pip freeze > "${WORKDIR}/python_packages.txt"
fi

echo
echo "=== installed ==="
"${PY_VENV}/bin/python" - <<'PY'
import Bio, scipy, pandas, numpy, sys
print("python    ", sys.version.split()[0])
print("biopython ", Bio.__version__)
print("scipy     ", scipy.__version__)
print("pandas    ", pandas.__version__)
print("numpy     ", numpy.__version__)
PY

echo
echo "use it with:  ${PY_VENV}/bin/python <script>"
