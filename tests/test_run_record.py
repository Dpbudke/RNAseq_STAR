import hashlib
import importlib.util
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

SPEC = importlib.util.spec_from_file_location(
    "record_run", Path(__file__).resolve().parents[1] / "scripts/record_run.py"
)
record = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(record)


class RunRecordChecks(unittest.TestCase):
    def test_source_snapshot_and_config_hash_are_bound_without_copying_arbitrary_settings(
        self,
    ):
        with tempfile.TemporaryDirectory(
            prefix="run_record_codetestfile_", dir="/tmp"
        ) as temp:
            root, work = Path(temp) / "pipeline", Path(temp) / "work"
            root.mkdir()
            work.mkdir()
            (work / "samples.txt").write_text("mouse_1\n")
            (work / "input_manifest.tsv").write_text("sample\tlane\nmouse_1\tlane1\n")
            for name in [
                "run.sh",
                "config.sh",
                "requirements.txt",
                "References/mouse_rrna.fasta",
            ]:
                path = root / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text("source contents\n")
            settings = Path(temp) / "profile.sh"
            settings.write_text('WORKDIR="/study/work"\n')
            with (
                patch.object(record.shutil, "which", return_value="/bin/bash"),
                patch.object(
                    record.subprocess,
                    "run",
                    return_value=subprocess.CompletedProcess([], 0, "version\n", ""),
                ),
            ):
                record.record(root, work, str(settings))
            saved = work / "run_records"
            self.assertEqual(
                (saved / "code/config.sh").read_bytes(),
                (root / "config.sh").read_bytes(),
            )
            expected = hashlib.sha256(settings.read_bytes()).hexdigest()
            self.assertIn(
                expected + "  " + str(settings), (saved / "code_sha256.txt").read_text()
            )
            self.assertFalse((saved / "code/profile.sh").exists())
            self.assertFalse(
                any(
                    p.read_bytes() == settings.read_bytes()
                    for p in saved.rglob("*")
                    if p.is_file()
                )
            )
            with self.assertRaisesRegex(ValueError, "already recorded"):
                record.record(root, work, str(settings))


if __name__ == "__main__":
    unittest.main()
