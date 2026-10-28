import importlib.util
import json
import os
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "local_run", ROOT / "scripts/local_run.py"
)
runner = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(runner)


class LocalExecutionChecks(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(
            prefix="local_run_codetestfile_", dir="/tmp"
        )
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name) / "pipeline"
        self.work = Path(self.temp.name) / "work"
        self.root.mkdir()
        self.work.mkdir()
        (self.work / "samples.txt").write_text("mouse_2\nmouse_10\n")
        (self.root / "run.sh").write_text("#!/bin/bash\nexit 0\n")
        for name in runner.STAGES.values():
            path = self.root / name
            path.parent.mkdir(exist_ok=True)
            path.write_text('test "$SLURM_CPUS_PER_TASK" = 3\n')
        self.args = [
            "local_run",
            "all",
            "--root",
            str(self.root),
            "--workdir",
            str(self.work),
            "--lanes",
            "3",
            "--threads",
            "3",
        ]

    def invoke(self, *extra):
        with patch("sys.argv", self.args + list(extra)):
            runner.main()

    def test_three_lanes_and_two_samples_run_in_dependency_order(self):
        self.invoke()
        report = json.loads((self.work / "local_jobs.json").read_text())
        self.assertEqual(
            [row["task"] for row in report["tasks"]],
            [
                "references/0",
                "pyenv/0",
                "md5/1",
                "md5/2",
                "md5/3",
                "concat/1",
                "concat/2",
                "hts/1",
                "hts/2",
                "index/0",
                "star/1",
                "star/2",
                "finalize/0",
            ],
        )
        self.assertTrue(
            all(
                row["state"] == "COMPLETED" and row["exit_code"] == 0
                for row in report["tasks"]
            )
        )
        original = (self.work / "local_jobs.json").read_bytes()
        with self.assertRaisesRegex(ValueError, "already recorded"):
            self.invoke()
        self.assertEqual((self.work / "local_jobs.json").read_bytes(), original)

    def test_failed_stage_prevents_downstream_work_and_keeps_stderr(self):
        (self.root / runner.STAGES["hts"]).write_text(
            "echo 'test failure' >&2\nexit 17\n"
        )
        with self.assertRaisesRegex(RuntimeError, "failed with exit 17"):
            self.invoke()
        report = json.loads((self.work / "local_jobs.json").read_text())
        self.assertEqual(report["status"], "failed")
        self.assertEqual(report["tasks"][-1]["task"], "hts/1")
        self.assertFalse(
            any(row["task"].startswith("star/") for row in report["tasks"])
        )
        self.assertEqual(
            (self.work / "slurmout/local_hts_1_1.err").read_text(), "test failure\n"
        )

    def test_dry_run_creates_no_job_record_or_stage_logs(self):
        self.invoke("--dry-run")
        self.assertFalse((self.work / "local_jobs.json").exists())
        self.assertEqual(list((self.work / "slurmout").iterdir()), [])

    def test_slurm_and_local_runs_cannot_share_an_output_directory(self):
        (self.work / "jobs.tsv").write_text("star\t123\n")
        with self.assertRaisesRegex(ValueError, "Slurm jobs are recorded"):
            self.invoke()

    def test_missing_configuration_is_an_error(self):
        env = os.environ.copy()
        env["PIPELINE_CONFIG"] = str(self.root / "missing.sh")
        result = subprocess.run(
            ["bash", "-c", 'source "$1"', "bash", str(ROOT / "config.sh")],
            env=env,
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("configuration file not found", result.stderr)

    def test_tool_module_cannot_replace_the_selected_python(self):
        approved, bundled = self.root / "approved", self.root / "bundled"
        approved.mkdir()
        bundled.mkdir()
        for directory, value in [
            (approved, "approved_python"),
            (bundled, "bundled_python"),
        ]:
            program = directory / "python3"
            program.write_text("#!/bin/bash\necho " + value + "\n")
            program.chmod(0o755)
        profile = self.root / "profile.sh"
        profile.write_text('SOFTWARE_MODE="modules"\n')
        env = os.environ.copy()
        env["PIPELINE_CONFIG"] = str(profile)
        command = "\n".join(
            [
                'approved="$2"; bundled="$3"',
                "module() {",
                '    if [[ "$1" == load ]]; then',
                '        case "$2" in',
                '            python_3/*) export PATH="$approved:$PATH";;',
                '            htstream/*) export PATH="$bundled:$PATH";;',
                "        esac",
                "    fi",
                "}",
                'source "$1"',
                'load_software "$MODULE_PYTHON"',
                'load_software "$MODULE_HTSTREAM"',
                "python3",
            ]
        )
        result = subprocess.run(
            [
                "bash",
                "-c",
                command,
                "bash",
                str(ROOT / "config.sh"),
                str(approved),
                str(bundled),
            ],
            env=env,
            text=True,
            capture_output=True,
            check=True,
        )
        self.assertEqual(result.stdout.strip(), "approved_python")


if __name__ == "__main__":
    unittest.main()
