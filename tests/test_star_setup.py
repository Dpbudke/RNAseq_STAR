import hashlib
import importlib.util
import io
import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

spec = importlib.util.spec_from_file_location(
    "star_setup", Path(__file__).resolve().parents[1] / "0_Setup/4_fetch_star.py"
)
setup = importlib.util.module_from_spec(spec)
spec.loader.exec_module(setup)


class StarSetupChecks(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(
            prefix="star_setup_codetestfile_", dir="/tmp"
        )
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.source = {
            "url": "https://example.invalid/STAR",
            "sha256": hashlib.sha256(b"valid").hexdigest(),
        }
        mask = os.umask(0o007)
        self.addCleanup(os.umask, mask)

    def test_corrupt_download_is_rejected_and_partial_removed(self):
        target = self.root / "STAR"
        with patch.object(
            setup.urllib.request, "urlopen", return_value=io.BytesIO(b"damaged")
        ):
            with self.assertRaisesRegex(ValueError, "checksum differs"):
                setup.fetch(target, self.source, 0o770)
        self.assertFalse(target.exists())
        self.assertEqual(list(self.root.iterdir()), [])

    def test_valid_download_and_cache_need_no_second_download(self):
        target = self.root / "STAR"
        with patch.object(
            setup.urllib.request, "urlopen", return_value=io.BytesIO(b"valid")
        ) as request:
            setup.fetch(target, self.source, 0o770)
            with patch.object(
                Path, "chmod", side_effect=PermissionError("different owner")
            ):
                setup.fetch(target, self.source, 0o770)
            self.assertEqual(request.call_count, 1)
        self.assertEqual(target.read_bytes(), b"valid")
        self.assertEqual(target.stat().st_mode & 0o777, 0o770)

    def test_existing_different_file_is_preserved(self):
        target = self.root / "STAR"
        target.write_bytes(b"other build")
        with patch.object(setup.urllib.request, "urlopen") as request:
            with self.assertRaisesRegex(ValueError, "Existing runtime file differs"):
                setup.fetch(target, self.source, 0o770)
            request.assert_not_called()
        self.assertEqual(target.read_bytes(), b"other build")

    def test_file_created_during_download_is_preserved(self):
        target = self.root / "STAR"

        def response(*args, **kwargs):
            target.write_bytes(b"another process")
            return io.BytesIO(b"valid")

        with patch.object(setup.urllib.request, "urlopen", side_effect=response):
            with self.assertRaises(FileExistsError):
                setup.fetch(target, self.source, 0o770)
        self.assertEqual(target.read_bytes(), b"another process")
        self.assertFalse(list(self.root.glob("*.part")))

    def test_existing_configuration_is_not_overwritten(self):
        target = self.root / "star_runtime.sh"
        setup.write_exact(target, "first\n")
        setup.write_exact(target, "first\n")
        with self.assertRaisesRegex(ValueError, "configuration differs"):
            setup.write_exact(target, "second\n")
        self.assertEqual(target.read_text(), "first\n")
        self.assertFalse(list(self.root.glob("*.part")))

    def test_symlink_cannot_redirect_installation_outside_tools(self):
        tools = self.root / "tools"
        outside = self.root / "outside"
        tools.mkdir()
        outside.mkdir()
        (tools / f"star_{setup.VERSION}").symlink_to(outside, target_is_directory=True)
        with (
            patch.object(setup.platform, "system", return_value="Linux"),
            patch.object(setup.platform, "machine", return_value="x86_64"),
        ):
            with self.assertRaisesRegex(ValueError, "outside the tools"):
                setup.install(tools)
        self.assertEqual(list(outside.iterdir()), [])

    def test_runtime_failure_stops_config_and_explicit_override_is_preserved(self):
        config = self.root / "config.sh"
        shutil.copyfile(Path(__file__).resolve().parents[1] / "config.sh", config)
        (self.root / "tools").mkdir()
        (self.root / "tools/star_runtime.sh").write_text("return 1\n")
        command = [
            "bash",
            "-c",
            'unset MODULE_STAR PIPELINE_CONFIG; source "$1"',
            "bash",
            str(config),
        ]
        self.assertEqual(subprocess.run(command, check=False).returncode, 1)
        (self.root / "run_settings.sh").write_text('MODULE_STAR="site_star/version"\n')
        command[2] += '; printf "%s\\n" "$MODULE_STAR"'
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        self.assertEqual(result.stdout.strip(), "site_star/version")


if __name__ == "__main__":
    unittest.main()
