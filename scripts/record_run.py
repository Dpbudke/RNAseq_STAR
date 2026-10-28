#!/usr/bin/env python3
"""Record the source files and executable identities before processing starts."""

import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def record(root, work, settings):
    destination = work / "run_records"
    destination.mkdir(exist_ok=True)
    manifest = destination / "code_sha256.txt"
    if manifest.exists():
        raise ValueError("Source files are already recorded; inspect the existing run")
    files = [
        root / name
        for name in [
            "run.sh",
            "config.sh",
            "requirements.txt",
            "References/mouse_rrna.fasta",
        ]
    ]
    for name in [
        "0_Setup",
        "1_Pre_Preprocessing",
        "2_HTStream_Preprocessing",
        "3_STAR_Alignment",
        "scripts",
    ]:
        files.extend(
            p for p in (root / name).rglob("*") if p.suffix in {".py", ".sh", ".slurm"}
        )
    if (root / "tools").exists():
        files.extend(
            p
            for p in (root / "tools").rglob("*")
            if p.is_file()
            and (
                p.name in {"STAR", "LICENSE", "star_runtime.json", "star_runtime.sh"}
                or p.suffix == ".lua"
            )
        )
    hashes = []
    for source in sorted(set(files)):
        if not source.is_file() or source.is_symlink():
            raise ValueError("Source file is missing or linked: " + str(source))
        saved = destination / "code" / source.relative_to(root)
        saved.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, saved)
        sha = digest(source)
        if digest(saved) != sha:
            raise ValueError("Source snapshot differs")
        hashes.append((sha, source))
    for name in ["samples.txt", "input_manifest.tsv"]:
        source = work / name
        shutil.copy2(source, destination / name)
        hashes.append((digest(source), source))
    if settings:
        # Bind the chosen profile without copying arbitrary shell content.
        hashes.append((digest(Path(settings)), Path(settings)))
    if any("\n" in str(path) or "\\" in str(path) for _, path in hashes):
        raise ValueError("Unsupported source path for checksum manifest")
    programs = {}
    if os.environ.get("SOFTWARE_MODE") == "container":
        # The tools run inside the image; its checksum fixes their versions (container/Dockerfile).
        image = Path(os.environ["STAR_SIF"]).resolve()
        programs["container"] = {"path": str(image), "sha256": digest(image)}
    for name in [] if programs else ["STAR", "samtools", "hts_Stats", "python3"]:
        executable = shutil.which(name)
        if executable is None:
            raise ValueError("Required executable is unavailable: " + name)
        result = subprocess.run(
            [executable, "--version"], text=True, capture_output=True, check=False
        )
        programs[name] = {
            "path": executable,
            "sha256": digest(Path(executable)),
            "version_output": (result.stdout + result.stderr).strip(),
        }
    (destination / "software.json").write_text(json.dumps(programs, indent=2) + "\n")
    manifest.write_text("".join(sha + "  " + str(path) + "\n" for sha, path in hashes))
    print(
        f"Recorded {len(hashes)} source/configuration hashes and {len(programs)} executable identities"
    )


if __name__ == "__main__":
    os.umask(0o007)
    record(Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve(), sys.argv[3])
