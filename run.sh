#!/bin/bash
# Single entry point. Edit config.sh, then:
#
#   ./run.sh scan      look at the lane directories, print what is there
#   ./run.sh check     validate config.sh, work out read length and index path
#   ./run.sh all       submit every stage, chained so each waits for the last
#   ./run.sh status    squeue
#
# Individual stages, if you want to drive it by hand:
#   references pyenv md5 samples concat hts htsqc index star stats counts
#   finalize confirm kmerqc compare
# With EXECUTOR=slurm every stage that reads data runs as a Slurm job, including
# the short ones (htsqc stats counts confirm kmerqc compare): a cluster login node
# is only for submitting and checking jobs.
#
# --dry-run prints the sbatch line without submitting.
# Anything after -- is passed through to sbatch:
#   ./run.sh hts -- --array=1-3

set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${HERE}/config.sh"

JOURNAL=""
DRY=0
STEP="${1:-}"
[[ -n "$STEP" ]] || { sed -n '2,18p' "${BASH_SOURCE[0]}" | sed 's/^# \?//'; exit 2; }
shift || true
EXTRA=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run) DRY=1; shift ;;
        --) shift; EXTRA+=("$@"); break ;;
        *) EXTRA+=("$1"); shift ;;
    esac
done

if [[ "$EXECUTOR" == local ]]; then
    case "$STEP" in
        all|references|pyenv|md5|concat|hts|index|star|finalize|status)
            [[ ${#EXTRA[@]} -eq 0 ]] || { echo "Local execution does not accept sbatch options" >&2; exit 2; }
            load_software "$MODULE_PYTHON"
            args=("$STEP" --root "$HERE" --workdir "$WORKDIR" --lanes "${#LANE_DIRS[@]}" --threads "$LOCAL_THREADS")
            [[ "$DRY" -eq 0 ]] || args+=(--dry-run)
            exec python3 "$HERE/scripts/local_run.py" "${args[@]}"
            ;;
    esac
fi

red()  { printf '  FAIL  %s\n' "$*"; }
ok()   { printf '  ok    %s\n' "$*"; }
warn() { printf '  note  %s\n' "$*"; }

sample_name() {                     # sample_name <fastq basename> <regex>
    local b="$1"
    b="${b%"${FASTQ_SUFFIX_R1}"}"
    b="${b%"${FASTQ_SUFFIX_R2}"}"
    sed -E "$2" <<< "$b"
}

read_length() {                     # read_length <fastq.gz>
    python3 - "$HERE" "$1" <<'PYREAD'
import sys
sys.path.insert(0, sys.argv[1] + '/scripts')
from sample_manifest import read_length
print(read_length(sys.argv[2]))
PYREAD
}

# ---------------------------------------------------------------- scan
do_scan() {
    load_software "${MODULE_PYTHON}"
    echo "=== lane directories ==="
    local bad=0
    for d in "${LANE_DIRS[@]}"; do
        if [[ -d "$d" ]]; then
            local n1 n2 tags
            n1=$(find "$d" -maxdepth 1 -name "*${FASTQ_SUFFIX_R1}" | wc -l)
            n2=$(find "$d" -maxdepth 1 -name "*${FASTQ_SUFFIX_R2}" | wc -l)
            tags=$(find "$d" -maxdepth 1 -name "*.fastq.gz" -printf '%f\n' 2>/dev/null \
                   | sed -nE 's/.*_(L[0-9]+)_R[12]_001\.fastq\.gz/\1/p' | sort -u | tr '\n' ' ')
            printf '  %s\n    R1=%s R2=%s  lane tags: %s\n' "$d" "$n1" "$n2" "${tags:-none}"
            [[ "$LIBRARY_LAYOUT" == single || "$n1" == "$n2" ]] || { red "R1 and R2 counts differ in $d"; bad=1; }
            [[ "$LIBRARY_LAYOUT" == paired || "$n2" == 0 ]] || warn "R2 files in $d are ignored (LIBRARY_LAYOUT=single)"
            [[ "$n1" -gt 0 ]] || { red "no ${FASTQ_SUFFIX_R1} files in $d"; bad=1; }
        else
            red "not a directory: $d"; bad=1
        fi
    done
    [[ $bad -eq 0 ]] || return 1

    echo
    echo "=== sample names, both styles ==="
    local first="${LANE_DIRS[0]}" f
    printf '  %-46s %-22s %s\n' "filename" "illumina" "first_field"
    local i=0
    while IFS= read -r f; do
        printf '  %-46s %-22s %s\n' "$(basename "$f")" \
            "$(sample_name "$(basename "$f")" 's/_S[0-9]+_L[0-9]+$//')" \
            "$(sample_name "$(basename "$f")" 's/_.*$//')"
        i=$((i+1)); [[ $i -ge 5 ]] && break
    done < <(find "$first" -maxdepth 1 -name "*${FASTQ_SUFFIX_R1}" | sort)

    echo
    local ni nf tot
    tot=$(find "$first" -maxdepth 1 -name "*${FASTQ_SUFFIX_R1}" | wc -l)
    ni=$(find "$first" -maxdepth 1 -name "*${FASTQ_SUFFIX_R1}" -printf '%f\n' \
         | while read -r f; do sample_name "$f" 's/_S[0-9]+_L[0-9]+$//'; done | sort -u | wc -l)
    nf=$(find "$first" -maxdepth 1 -name "*${FASTQ_SUFFIX_R1}" -printf '%f\n' \
         | while read -r f; do sample_name "$f" 's/_.*$//'; done | sort -u | wc -l)
    printf '  %s files -> illumina gives %s distinct names, first_field gives %s\n' "$tot" "$ni" "$nf"
    [[ "$ni" -lt "$tot" ]] && warn "illumina style collapses names, check the table above"
    [[ "$nf" -lt "$tot" ]] && warn "first_field style collapses names, check the table above"
    echo "  config.sh has SAMPLE_NAME_STYLE=\"${SAMPLE_NAME_STYLE}\""

    echo
    echo "=== excluded by EXCLUDE_SAMPLES ==="
    local found=0
    for e in "${EXCLUDE_SAMPLES[@]:-}"; do
        [[ -z "$e" ]] && continue
        local hits
        hits=$(find "$first" -maxdepth 1 -name "${e}*" | wc -l)
        [[ "$hits" -gt 0 ]] && { printf '  %s  (%s files in lane 1)\n' "$e" "$hits"; found=1; }
    done
    [[ $found -eq 0 ]] && echo "  none of the excluded names are present"

    echo
    echo "=== read length ==="
    local one rl
    one=$(find "$first" -maxdepth 1 -name "*${FASTQ_SUFFIX_R1}" | sort | awk 'NR==1{print}')
    rl=$(read_length "$one")
    printf '  %s bp  (from %s)\n' "$rl" "$(basename "$one")"
    printf '  --sjdbOverhang will be %s\n' "$((rl - 1))"
}

# ---------------------------------------------------------------- check
do_check() {
    local bad=0
    echo "=== config ==="
    [[ "${WORKDIR}"      == *CHANGEME* ]] && { red "WORKDIR still says CHANGEME"; bad=1; } || ok "WORKDIR ${WORKDIR}"
    if [[ "$EXECUTOR" == slurm ]]; then
        command -v sbatch >/dev/null || { red "sbatch is unavailable"; bad=1; }
        [[ "$SLURM_ACCOUNT" != *CHANGEME* ]] || { red "SLURM_ACCOUNT still says CHANGEME"; bad=1; }
        [[ "$SLURM_MAIL" != *CHANGEME* ]] || { red "Set SLURM_MAIL to an address or an empty string"; bad=1; }
        if [[ -n "$SLURM_ACCOUNT" ]] && command -v sacctmgr >/dev/null; then
            if sacctmgr -nP show assoc user="$USER" format=Account 2>/dev/null | sort -u | grep -qx "$SLURM_ACCOUNT"; then
                ok "account $SLURM_ACCOUNT"
            else
                red "account $SLURM_ACCOUNT is not available to $USER"; bad=1
            fi
        fi
        # On SCINet a project is the Slurm account and its storage, so the run writes only under that project.
        if [[ -n "$SLURM_ACCOUNT" && -d /90daydata ]]; then
            for d in "$WORKDIR" "$REF_DIR"; do
                case "$(realpath -m "$d")/" in
                    /project/"$SLURM_ACCOUNT"/*|/90daydata/"$SLURM_ACCOUNT"/*) ok "$d is in project $SLURM_ACCOUNT" ;;
                    *) red "$d is not under /project/$SLURM_ACCOUNT or /90daydata/$SLURM_ACCOUNT"; bad=1 ;;
                esac
            done
        fi
        if [[ -n "$SLURM_PARTITION" ]] && command -v scontrol >/dev/null; then
            if scontrol show partition 2>/dev/null | grep -q "PartitionName=${SLURM_PARTITION}\b"; then
                ok "partition $SLURM_PARTITION"
            else
                red "partition $SLURM_PARTITION not found"; bad=1
            fi
        fi
    fi

    echo "=== lanes ($LIBRARY_LAYOUT-end) ==="
    for d in "${LANE_DIRS[@]}"; do
        [[ -d "$d" ]] && ok "$d" || { red "missing $d"; bad=1; }
        if [[ "$SOURCE_MD5" == none ]]; then
            warn "SOURCE_MD5=none: $d is not checked against provider checksums (SHA-256 is recorded)"
        elif [[ -d "$d" && ! -s "$d/$SOURCE_MD5" ]]; then
            red "no $SOURCE_MD5 in $d (set SOURCE_MD5 to the provider's checksum file, or none)"; bad=1
        fi
    done
    [[ ${#LANE_DIRS[@]} -ge 1 ]] || { red "LANE_DIRS is empty"; bad=1; }

    echo "=== software ==="
    if [[ "$SOFTWARE_MODE" == modules ]]; then
        for m in "$MODULE_PYTHON" "$MODULE_HTSTREAM" "$MODULE_STAR" "$MODULE_SAMTOOLS"; do
            if (load_software "$m") >/dev/null 2>&1; then ok "$m"; else red "cannot load module $m"; bad=1; fi
        done
    elif [[ "$SOFTWARE_MODE" == container ]]; then
        # Apptainer may not run on a login node; the jobs run the image. Here: the image is there, apptainer exists.
        [[ -r "$STAR_SIF" ]] && ok "image $STAR_SIF" || { red "STAR_SIF is not a readable image: ${STAR_SIF:-unset}"; bad=1; }
        if command -v apptainer >/dev/null || { type module >/dev/null 2>&1 && module -t avail "$MODULE_APPTAINER" 2>&1 | grep -q .; }; then
            ok "apptainer"
        else
            red "apptainer is unavailable (MODULE_APPTAINER=$MODULE_APPTAINER)"; bad=1
        fi
    else
        for tool in python3 STAR samtools hts_Stats hts_SeqScreener hts_SuperDeduper hts_AdapterTrimmer hts_PolyATTrim hts_NTrimmer hts_QWindowTrim hts_LengthFilter md5sum zcat; do
            command -v "$tool" >/dev/null || { red "$tool is not on PATH"; bad=1; }
        done
    fi
    if [[ "$(uname -s)" != Linux ]]; then
        red "Processing requires Linux and GNU command-line tools"; bad=1
    fi

    if load_software "$MODULE_PYTHON"; then
        python3 - "$WORKDIR" "${LANE_DIRS[@]}" <<'PYCHECK' || bad=1
import sys
from pathlib import Path
if sys.version_info < (3, 11):
    raise SystemExit('FAIL: Python 3.11 or newer is required')
work = Path(sys.argv[1]).resolve()
for value in sys.argv[2:]:
    lane = Path(value).resolve()
    if work == lane or work in lane.parents or lane in work.parents:
        raise SystemExit('FAIL: WORKDIR and raw lane directories must not overlap')
PYCHECK
    else
        bad=1
    fi
    [[ $bad -eq 0 ]] || { echo "check failed; fix the FAIL lines above"; return 1; }
    echo "=== workdir ==="
    mkdir -p "${WORKDIR}" 2>/dev/null || true
    if [[ -d "${WORKDIR}" && -w "${WORKDIR}" ]]; then ok "writable"; else red "${WORKDIR} not writable"; bad=1; fi
    local avail
    avail=$(df -BG --output=avail "${WORKDIR}" 2>/dev/null | tail -1 | tr -dc '0-9')
    [[ -n "$avail" ]] && ok "${avail} GB free on that filesystem"

    [[ $bad -eq 0 ]] || { echo; echo "check failed, fix the FAIL lines above"; return 1; }

    load_software "${MODULE_PYTHON}"
    # resolve the auto settings
    echo "=== resolving ==="
    local one rl overhang
    one=$(find "${LANE_DIRS[0]}" -maxdepth 1 -name "*${FASTQ_SUFFIX_R1}" | sort | awk 'NR==1{print}')
    rl=$(read_length "$one")
    if [[ "${STAR_SJDB_OVERHANG}" == "auto" ]]; then overhang=$((rl - 1)); else overhang="${STAR_SJDB_OVERHANG}"; fi
    ok "read length ${rl} bp, --sjdbOverhang ${overhang}"
    if [[ "$LIBRARY_LAYOUT" == single && "$HTS_DEDUP" == yes ]]; then
        warn "HTS_DEDUP=yes on single-end reads: without UMIs, reads starting at the same place are often not PCR"
        warn "duplicates (highly expressed genes), so many may be removed; HTS_DEDUP=no keeps them"
    fi
    if (( rl < HTS_MIN_LENGTH )); then
        red "reads (${rl} bp) are shorter than HTS_MIN_LENGTH (${HTS_MIN_LENGTH}); every read would be dropped"; return 1
    elif (( rl < HTS_MIN_LENGTH + 20 )); then
        warn "reads (${rl} bp) are close to HTS_MIN_LENGTH (${HTS_MIN_LENGTH}); trimmed reads below it are dropped"
    fi

    mkdir -p "${WORKDIR}"
    {
        echo "# written by ./run.sh check on $(date '+%Y-%m-%d %H:%M')"
        echo "READ_LENGTH=${rl}"
        echo "STAR_SJDB_OVERHANG=${overhang}"
        echo "STAR_INDEX=\"${REF_DIR}/star.overlap${overhang}.${REFERENCE_TAG}\""
    } > "${WORKDIR}/resolved.sh"
    ok "wrote ${WORKDIR}/resolved.sh"
    ok "STAR index will be ${REF_DIR}/star.overlap${overhang}.${REFERENCE_TAG}"
    echo
    echo "config looks good. next: ./run.sh all"
}

record_run() (
    load_software "$MODULE_PYTHON" "$MODULE_HTSTREAM" "$MODULE_STAR" "$MODULE_SAMTOOLS"
    mkdir -p "$WORKDIR/run_records"
    declare -p WORKDIR LANE_DIRS REF_DIR REFERENCE_SOURCE REFERENCE_TAG GENCODE_RELEASE ENSEMBL_SPECIES \
        ENSEMBL_RELEASE CUSTOM_NAME CUSTOM_GENOME_FA CUSTOM_GTF RRNA_FASTA STAR_SJDB_OVERHANG COUNT_COLUMN \
        SAMPLE_NAME_STYLE SAMPLE_REGEX EXCLUDE_SAMPLES FASTQ_SUFFIX_R1 FASTQ_SUFFIX_R2 LIBRARY_LAYOUT SOURCE_MD5 HTS_MIN_LENGTH HTS_DEDUP \
        EXECUTOR SOFTWARE_MODE STAR_SIF LOCAL_THREADS SLURM_ACCOUNT SLURM_PARTITION SLURM_MAIL \
        MODULE_STAR MODULE_HTSTREAM MODULE_PYTHON MODULE_SAMTOOLS \
        > "$WORKDIR/run_records/settings_snapshot.sh"
    SOFTWARE_MODE="$SOFTWARE_MODE" STAR_SIF="$STAR_SIF" python3 "$HERE/scripts/record_run.py" "$HERE" "$WORKDIR" "${PIPELINE_CONFIG:-}"
)

# ---------------------------------------------------------------- submit
common_opts() {
    [[ -z "$SLURM_ACCOUNT" ]] || printf '%s\n' "--account=$SLURM_ACCOUNT"
    [[ -z "$SLURM_PARTITION" ]] || printf '%s\n' "--partition=$SLURM_PARTITION"
    [[ -z "$SLURM_MAIL" ]] || printf '%s\n' "--mail-user=$SLURM_MAIL" "--mail-type=$SLURM_MAIL_TYPE"
    printf '%s\n' "--chdir=$WORKDIR" "--kill-on-invalid-dep=yes" "--export=ALL,PIPELINE_ROOT=$HERE"
}

n_samples() {
    if [[ ! -s "${SAMPLES_TXT}" ]]; then
        # under --dry-run samples.txt has not been written yet, so show a
        # placeholder rather than stopping
        if [[ $DRY -eq 1 ]]; then echo "N"; return; fi
        echo "ERROR: ${SAMPLES_TXT} missing or empty. Run './run.sh samples' first." >&2
        exit 1
    fi
    grep -c . "${SAMPLES_TXT}"
}
sb() {                              # sb <script> [opts...]  -> prints job id
    local script="$1"; shift
    local -a cmd=(sbatch --parsable)
    while IFS= read -r o; do cmd+=("$o"); done < <(common_opts)
    case "$script" in
        3_STAR_Alignment/2_star.slurm) cmd+=("--cpus-per-task=${STAR_CPUS:-8}" "--mem=${STAR_MEMORY:-48G}");;
        3_STAR_Alignment/1_star_index.slurm) cmd+=("--cpus-per-task=${INDEX_CPUS:-16}" "--mem=${INDEX_MEMORY:-128G}");;
    esac
    cmd+=("$@"); cmd+=("${EXTRA[@]+"${EXTRA[@]}"}"); cmd+=("${script}")
    if [[ $DRY -eq 1 ]]; then echo "+ ${cmd[*]}" >&2; echo "DRYRUN"; return; fi
    echo "+ ${cmd[*]}" >&2
    cd "${HERE}"
    local job
    job=$("${cmd[@]}")
    [[ "$job" =~ ^[0-9]+([\;].*)?$ ]] || { echo "Invalid sbatch response" >&2; exit 1; }
    [[ -z "$JOURNAL" ]] || printf '%s\t%s\n' "$script" "$job" >> "$JOURNAL"
    echo "$job"
}
loc() { echo "+ bash $*" >&2; [[ $DRY -eq 1 ]] && return 0; cd "${HERE}"; bash "$@"; }

# Short summaries and checks: a small Slurm job with EXECUTOR=slurm, run here with EXECUTOR=local.
short_job() {                       # short_job NAME CPUS MEMORY TIME COMMAND...
    local name="$1" cpus="$2" mem="$3" time="$4"; shift 4
    if [[ "$EXECUTOR" == local ]]; then
        echo "+ $*" >&2; [[ $DRY -eq 1 ]] && return 0
        (cd "${HERE}" && "$@"); return
    fi
    local -a cmd=(sbatch --parsable "--job-name=${name}" "--cpus-per-task=${cpus}" "--mem=${mem}"
                  "--time=${time}" "--output=${SLURMOUT}/${name}_%j.out")
    while IFS= read -r o; do cmd+=("$o"); done < <(common_opts)
    cmd+=(--wrap "cd $(printf '%q' "$HERE") && $(printf '%q ' "$@")")
    if [[ $DRY -eq 1 ]]; then echo "+ ${cmd[*]}" >&2; return 0; fi
    local job
    job=$("${cmd[@]}")
    [[ "$job" =~ ^[0-9]+([\;].*)?$ ]] || { echo "Invalid sbatch response" >&2; exit 1; }
    job=${job%%;*}
    [[ -z "$JOURNAL" ]] || printf '%s\t%s\n' "$name" "$job" >> "$JOURNAL"
    echo "Submitted ${name} as job ${job}; output ${SLURMOUT}/${name}_${job}.out"
}
need_venv() { [[ -x "${PY_VENV}/bin/python" ]] || { echo "no venv, run './run.sh pyenv' first" >&2; exit 1; }; }

mkdir -p "${SLURMOUT}" 2>/dev/null || true

case "$STEP" in
    scan)   do_scan ;;
    check)  do_check ;;
    record) record_run ;;

    references) sb "0_Setup/1_fetch_references.slurm" ;;
    pyenv)      sb "0_Setup/2_python_env.slurm" ;;
    md5)        sb "0_Setup/3_verify_source_md5.slurm" "--array=1-${#LANE_DIRS[@]}" ;;
    samples)    loc "1_Pre_Preprocessing/1_make_samples_txt.sh" ;;
    concat)     n=$(n_samples); sb "1_Pre_Preprocessing/2_concatenate_rnaseq.slurm" "--array=1-${n}%${THROTTLE_CONCAT}" ;;
    hts)        n=$(n_samples); sb "2_HTStream_Preprocessing/1_hts_preproc.slurm" "--array=1-${n}%${THROTTLE_HTS}" ;;
    index)      sb "3_STAR_Alignment/1_star_index.slurm" ;;
    star)       n=$(n_samples); sb "3_STAR_Alignment/2_star.slurm" "--array=1-${n}%${THROTTLE_STAR}" ;;
    stats)      short_job star_stats 2 8G 01:00:00 bash 3_STAR_Alignment/3_star_stats.sh ;;
    counts)     short_job star_counts 2 16G 02:00:00 bash 3_STAR_Alignment/4_reads2counts.sh ;;
    finalize)   sb "3_STAR_Alignment/7_finalize.slurm" ;;
    htsqc)
        need_venv
        short_job hts_summary 2 8G 01:00:00 env HTS_OUT="${HTS_OUT}" SAMPLES_TXT="${SAMPLES_TXT}" \
            "${PY_VENV}/bin/python" 2_HTStream_Preprocessing/3_hts_summary.py \
            --out "${WORKDIR}/summary_htstream.txt"
        ;;

    all)
        if [[ $DRY -eq 0 ]]; then
            exec 8>"${WORKDIR}/submission.lock"
            flock -n 8 || { echo "A submission is already in progress" >&2; exit 1; }
            [[ ! -s "${WORKDIR}/jobs.tsv" ]] || { echo "jobs.tsv already exists; inspect existing jobs before submitting replacements" >&2; exit 1; }
            JOURNAL="${WORKDIR}/jobs.tsv"
        fi
        do_check || exit 1
        source "${WORKDIR}/resolved.sh"
        echo
        echo "=== submitting the whole chain ==="
        loc "1_Pre_Preprocessing/1_make_samples_txt.sh" >/dev/null
        [[ $DRY -eq 1 ]] || record_run
        j_ref=$(sb "0_Setup/1_fetch_references.slurm")
        j_env=$(sb "0_Setup/2_python_env.slurm")
        n=$(n_samples)
        echo "  ${n} samples"
        j_md5=$(sb "0_Setup/3_verify_source_md5.slurm" "--array=1-${#LANE_DIRS[@]}")
        j_cat=$(sb "1_Pre_Preprocessing/2_concatenate_rnaseq.slurm" "--array=1-${n}%${THROTTLE_CONCAT}" "--dependency=afterok:${j_md5}")
        j_hts=$(sb "2_HTStream_Preprocessing/1_hts_preproc.slurm"   "--array=1-${n}%${THROTTLE_HTS}"  "--dependency=afterok:${j_cat},afterok:${j_ref}")
        j_idx=$(sb "3_STAR_Alignment/1_star_index.slurm"            "--dependency=afterok:${j_ref}")
        j_str=$(sb "3_STAR_Alignment/2_star.slurm"                  "--array=1-${n}%${THROTTLE_STAR}" "--dependency=afterok:${j_hts},afterok:${j_idx}")
        j_fin=$(sb "3_STAR_Alignment/7_finalize.slurm"              "--dependency=afterok:${j_str},afterok:${j_env}")
        echo
        printf '  references %s\n  pyenv      %s\n  md5        %s\n  concat     %s\n  htstream   %s\n  star index %s\n  star align %s\n  stats+cnt  %s\n' \
            "$j_ref" "$j_env" "$j_md5" "$j_cat" "$j_hts" "$j_idx" "$j_str" "$j_fin"
        echo
        if [[ $DRY -eq 0 ]]; then
            printf 'references\t%s\npyenv\t%s\nmd5\t%s\nconcat\t%s\nhts\t%s\nindex\t%s\nstar\t%s\nfinalize\t%s\n' \
                "$j_ref" "$j_env" "$j_md5" "$j_cat" "$j_hts" "$j_idx" "$j_str" "$j_fin" > "${WORKDIR}/jobs.tsv"
        fi
        echo "  watch it:  ./run.sh status"
        echo "  output:    ${WORKDIR}"
        ;;

    confirm)
        need_venv
        short_job star_confirm 2 16G 02:00:00 env STAR_OUT="${STAR_OUT}" COUNTS_FILE="${COUNTS_FILE}" \
            "${PY_VENV}/bin/python" 3_STAR_Alignment/5_confirm_reads2counts.py \
            --star-dir "$STAR_OUT" --counts "$COUNTS_FILE" --samples "$SAMPLES_TXT" \
            --column "$COUNT_COLUMN" --report "$WORKDIR/count_validation.json"
        ;;
    kmerqc)
        need_venv
        short_job hts_kmer_qc 4 32G 12:00:00 bash -c 'cd "$1" && shift && exec "$@"' kmerqc "${WORKDIR}" \
            env CONCAT_DIR="${CONCAT_DIR}" HTS_OUT="${HTS_OUT}" SAMPLES_TXT="${SAMPLES_TXT}" \
            "${PY_VENV}/bin/python" "${HERE}/2_HTStream_Preprocessing/2_HTStream_Kmer_QC.py" \
            --layout "$LIBRARY_LAYOUT" "${EXTRA[@]+"${EXTRA[@]}"}"
        ;;
    compare)
        [[ ${#EXTRA[@]} -ge 1 ]] || { echo "usage: ./run.sh compare -- /path/to/other_counts.txt" >&2; exit 2; }
        need_venv
        short_job star_compare 2 16G 01:00:00 "${PY_VENV}/bin/python" 3_STAR_Alignment/6_compare_to_published.py \
            --mine "${COUNTS_FILE}" --published "${EXTRA[0]}" \
            --out "${WORKDIR}/comparison.txt"
        ;;
    status) squeue -u "$USER" -o "%.14i %.9P %.12j %.9T %.10M %.10L %R" ;;
    *) echo "unknown step: $STEP" >&2; exit 2 ;;
esac
