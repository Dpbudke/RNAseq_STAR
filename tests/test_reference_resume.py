import hashlib
import importlib.util
import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

SPEC = importlib.util.spec_from_file_location(
    "fetch_refs", Path(__file__).resolve().parents[1] / "scripts/fetch_references.py"
)
refs = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(refs)


class ReferenceResumeChecks(unittest.TestCase):
    def test_interrupted_verification_keeps_other_cached_reference_records(self):
        with tempfile.TemporaryDirectory(
            prefix="reference_codetestfile_", dir="/tmp"
        ) as directory:
            root = Path(directory)
            names = ["GRCm39.primary_assembly.genome.fa", "gencode.vM35.annotation.gtf"]
            # A minimal genome and annotation that share a chromosome name (the fetcher checks that at the end).
            contents = {names[0]: b">chr1\nACGT\n", names[1]: b'chr1\tT\texon\t1\t4\t.\t+\t.\tgene_id "G";\n'}
            for name in names:
                (root / name).write_bytes(contents[name])
            rrna = root / "rrna_source.fasta"
            rrna.write_bytes(b">rRNA\nACGT\n")
            source_md5 = "a" * 32
            recorded = {
                name: {
                    "source_md5": source_md5,
                    "sha256": hashlib.sha256((root / name).read_bytes()).hexdigest(),
                }
                for name in names
            }
            report = root / "reference_manifest.json"
            report.write_text(json.dumps({"release": "M35", "files": recorded}))
            provider = "".join(
                source_md5 + "  " + name + ".gz\n" for name in names
            ).encode()
            original_digest = refs.digest

            def interrupt(path, algorithm="sha256"):
                if path.name == names[1]:
                    raise OSError("interrupted read")
                return original_digest(path, algorithm)

            args = [
                "fetch_references",
                "--directory",
                str(root),
                "--release",
                "M35",
                "--rrna",
                str(rrna),
            ]
            with (
                patch("sys.argv", args),
                patch.object(
                    refs.urllib.request, "urlopen", return_value=io.BytesIO(provider)
                ),
                patch.object(refs, "digest", side_effect=interrupt),
            ):
                with self.assertRaisesRegex(OSError, "interrupted read"):
                    refs.main()
            self.assertEqual(
                json.loads(report.read_text())["files"][names[1]], recorded[names[1]]
            )
            with (
                patch("sys.argv", args),
                patch.object(
                    refs.urllib.request, "urlopen", return_value=io.BytesIO(provider)
                ) as request,
            ):
                refs.main()
                self.assertEqual(request.call_count, 1)
            self.assertEqual(
                set(json.loads(report.read_text())["files"]),
                {*names, "mouse_rrna.fasta"},
            )


if __name__ == "__main__":
    unittest.main()
