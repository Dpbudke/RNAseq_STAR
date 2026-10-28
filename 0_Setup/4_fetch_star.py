#!/usr/bin/env python3
"""Install the tested Linux STAR release inside a pipeline copy."""

import argparse
import fcntl
import hashlib
import json
import os
import platform
import tempfile
import urllib.request
from pathlib import Path

VERSION = "2.7.11b"
COMMIT = "b1edc1208d91a53bf40ebae8669f71d50b994851"
BASE = f"https://raw.githubusercontent.com/alexdobin/STAR/{COMMIT}"
FILES = {
    f"star_{VERSION}/STAR": {
        "url": BASE + "/bin/Linux_x86_64_static/STAR",
        "sha256": "36e94b899a56b0ea5de5d65e722f55bc552b6713bd2d1e83a5c2abbd06c9881a",
    },
    f"star_{VERSION}/LICENSE": {
        "url": BASE + "/LICENSE",
        "sha256": "551e5bc4f91bf85b75c848aabba5acbcff4f4adf75e606cd17787038a155e69f",
    },
}
MODULE = f"""help("STAR {VERSION}: pinned Linux release")
local root = myFileName():match("^(.*)/modulefiles/[^/]+/[^/]+$")
prepend_path("PATH", pathJoin(root, "star_{VERSION}"))
"""
ACTIVATE = f'''# Source from run_settings.sh before submitting a new run.
_RNA_STAR_TOOLS="$(cd -- "$(dirname -- "${{BASH_SOURCE[0]}}")" && pwd)"
_RNA_STAR_HASH="$(sha256sum -- "$_RNA_STAR_TOOLS/star_{VERSION}/STAR")" || return 1
if [[ "${{_RNA_STAR_HASH%% *}}" != "{FILES[f"star_{VERSION}/STAR"]["sha256"]}" ]]; then
    echo "STAR executable differs from the recorded release" >&2
    return 1
fi
if [[ "${{SOFTWARE_MODE:-modules}}" == path ]]; then
    export PATH="$_RNA_STAR_TOOLS/star_{VERSION}:$PATH"
else
    module use "$_RNA_STAR_TOOLS/modulefiles" || return 1
    module unload rnaseq_star/{VERSION} || return 1
fi
MODULE_STAR="rnaseq_star/{VERSION}"
unset _RNA_STAR_HASH _RNA_STAR_TOOLS
'''


def digest(path):
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def make_directory(path):
    missing = []
    while not path.exists():
        missing.append(path)
        path = path.parent
    for path in reversed(missing):
        path.mkdir(mode=0o2770)
        path.chmod(0o2770)


def fetch(path, source, mode):
    if path.exists():
        if path.is_symlink() or digest(path) != source["sha256"]:
            raise ValueError("Existing runtime file differs: " + str(path))
        if path.stat().st_mode & 0o7777 != mode:
            path.chmod(mode)
        return
    make_directory(path.parent)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            dir=path.parent, suffix=".part", delete=False
        ) as stream:
            temporary = Path(stream.name)
            with urllib.request.urlopen(source["url"], timeout=30) as response:
                while block := response.read(1024 * 1024):
                    stream.write(block)
        if digest(temporary) != source["sha256"]:
            raise ValueError("Downloaded runtime checksum differs: " + str(path))
        temporary.chmod(mode)
        os.link(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def write_exact(path, text):
    make_directory(path.parent)
    if path.exists():
        if path.is_symlink() or path.read_text() != text:
            raise ValueError("Existing runtime configuration differs: " + str(path))
        return
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", dir=path.parent, suffix=".part", delete=False
        ) as stream:
            temporary = Path(stream.name)
            stream.write(text)
        temporary.chmod(0o660)
        os.link(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def install(root):
    if platform.system() != "Linux" or platform.machine() != "x86_64":
        raise ValueError("Run this setup on a Linux x86_64 host, such as Ceres")
    os.umask(0o007)
    root = root.resolve()
    make_directory(root)
    destinations = [
        *FILES,
        f"modulefiles/rnaseq_star/{VERSION}.lua",
        "star_runtime.sh",
        "star_runtime.json",
        ".star_setup.lock",
    ]
    if any(not (root / name).resolve().is_relative_to(root) for name in destinations):
        raise ValueError("Runtime destination resolves outside the tools directory")
    with (root / ".star_setup.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        for relative, source in FILES.items():
            fetch(
                root / relative, source, 0o770 if relative.endswith("/STAR") else 0o660
            )
        write_exact(root / f"modulefiles/rnaseq_star/{VERSION}.lua", MODULE)
        write_exact(root / "star_runtime.sh", ACTIVATE)
        write_exact(
            root / "star_runtime.json",
            json.dumps(
                {
                    "version": VERSION,
                    "git_commit": COMMIT,
                    "platform": "Linux x86_64",
                    "files": FILES,
                },
                indent=2,
            )
            + "\n",
        )
    print(f"Verified STAR {VERSION}: {root / f'star_{VERSION}/STAR'}")
    print("Setup records: " + str(root / "star_runtime.json"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--tools-dir", type=Path, default=Path(__file__).resolve().parents[1] / "tools"
    )
    arguments = parser.parse_args()
    try:
        install(arguments.tools_dir)
    except (ValueError, OSError) as error:
        parser.exit(2, f"ERROR: {error}\n")


if __name__ == "__main__":
    main()
