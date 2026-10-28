#!/usr/bin/env python3
"""Run the numbered pipeline stages sequentially without a scheduler."""

import argparse
import fcntl
import json
import os
import subprocess
import time
from pathlib import Path

STAGES = {
    "references": "0_Setup/1_fetch_references.slurm",
    "pyenv": "0_Setup/2_python_env.slurm",
    "md5": "0_Setup/3_verify_source_md5.slurm",
    "concat": "1_Pre_Preprocessing/2_concatenate_rnaseq.slurm",
    "hts": "2_HTStream_Preprocessing/1_hts_preproc.slurm",
    "index": "3_STAR_Alignment/1_star_index.slurm",
    "star": "3_STAR_Alignment/2_star.slurm",
    "finalize": "3_STAR_Alignment/7_finalize.slurm",
}


def save(path, value):
    partial = path.with_suffix(".part")
    partial.write_text(json.dumps(value, indent=2) + "\n")
    partial.replace(path)


def tasks(stage, work, lanes):
    if stage == "md5":
        return range(1, lanes + 1)
    if stage in {"concat", "hts", "star"}:
        samples = (work / "samples.txt").read_text().splitlines()
        if not samples or len(set(samples)) != len(samples):
            raise ValueError("Missing or duplicated sample list; run ./run.sh samples")
        return range(1, len(samples) + 1)
    return [0]


def run_stage(root, work, stage, task, threads, state, journal):
    key = f"{stage}/{task}"
    prior = sum(row["task"] == key for row in state["tasks"])
    log = work / "slurmout" / f"local_{stage}_{task}_{prior + 1}"
    env = os.environ.copy()
    env.update(PIPELINE_ROOT=str(root), SLURM_CPUS_PER_TASK=str(threads))
    if task:
        env["SLURM_ARRAY_TASK_ID"] = str(task)
    else:
        env.pop("SLURM_ARRAY_TASK_ID", None)
    row = {
        "task": key,
        "state": "RUNNING",
        "started_epoch": time.time(),
        "log": str(log),
    }
    state["tasks"].append(row)
    state["status"] = "running"
    save(journal, state)
    print(f"{key}: started; log {log}.out", flush=True)
    with (
        log.with_suffix(".out").open("x") as output,
        log.with_suffix(".err").open("x") as error,
    ):
        result = subprocess.run(
            ["bash", str(root / STAGES[stage])],
            cwd=work,
            env=env,
            stdout=output,
            stderr=error,
            check=False,
        )
    row.update(
        state="COMPLETED" if result.returncode == 0 else "FAILED",
        exit_code=result.returncode,
        finished_epoch=time.time(),
    )
    save(journal, state)
    if result.returncode:
        state["status"] = "failed"
        save(journal, state)
        raise RuntimeError(
            f"{key} failed with exit {result.returncode}; read {log}.err and .out"
        )
    print(f"{key}: completed", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("step", choices=["all", "status", *STAGES])
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--lanes", type=int, required=True)
    parser.add_argument("--threads", type=int, required=True)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    if args.lanes < 1 or args.threads < 1:
        parser.error("Lane and thread counts must be positive")
    root, work = args.root.resolve(), args.workdir.resolve()
    journal = work / "local_jobs.json"
    if args.step == "status":
        if not journal.exists():
            print("No local run recorded")
            return
        from collections import Counter

        state = json.loads(journal.read_text())
        latest = {row["task"]: row for row in state["tasks"]}
        print("Recorded local status: " + state["status"])
        for stage in STAGES:
            rows = [row for key, row in latest.items() if key.split("/")[0] == stage]
            if rows:
                print(
                    stage
                    + ": "
                    + ", ".join(
                        str(count) + " " + name
                        for name, count in sorted(
                            Counter(row["state"] for row in rows).items()
                        )
                    )
                )
        print("Full task and log record: " + str(journal))
        return
    os.umask(0o007)
    work.mkdir(parents=True, exist_ok=True)
    (work / "slurmout").mkdir(exist_ok=True)
    with (work / "submission.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if (work / "jobs.tsv").exists():
            raise ValueError(
                "Slurm jobs are recorded here; use a separate local run directory"
            )
        if args.step == "all" and journal.exists() and not args.dry_run:
            raise ValueError(
                "A local run is already recorded; inspect it and recover individual stages"
            )
        state = (
            json.loads(journal.read_text())
            if journal.exists()
            else {"executor": "local", "tasks": []}
        )
        if args.step == "all":
            for step in ["check", "samples"]:
                subprocess.run(
                    ["bash", str(root / "run.sh"), step], cwd=root, check=True
                )
        if args.step == "all" and not args.dry_run:
            subprocess.run(
                ["bash", str(root / "run.sh"), "record"], cwd=root, check=True
            )
        stages = list(STAGES) if args.step == "all" else [args.step]
        for stage in stages:
            for task in tasks(stage, work, args.lanes):
                if args.dry_run:
                    print(
                        f"{stage}/{task}: bash {STAGES[stage]} ({args.threads} STAR threads)"
                    )
                else:
                    run_stage(root, work, stage, task, args.threads, state, journal)
        if not args.dry_run:
            state["status"] = "completed"
            save(journal, state)
            print(
                "Completed requested stages. Final acceptance also requires pipeline_state.txt=validated and the count/BAM reports."
            )


if __name__ == "__main__":
    main()
