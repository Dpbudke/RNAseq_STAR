import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path

SPEC = importlib.util.spec_from_file_location(
    "manifest_inputs",
    Path(__file__).resolve().parents[1] / "scripts/manifest_inputs.py",
)
inputs = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(inputs)


class ManifestInputChecks(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(
            prefix="manifest_inputs_codetestfile_", dir="/tmp"
        )
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.lanes = [self.root / "lane1", self.root / "lane2"]
        self.rows = []
        for lane in self.lanes:
            lane.mkdir()
            mates = [
                lane / ("prefix_mouse_10_" + mate + ".fastq.gz")
                for mate in ["R1", "R2"]
            ]
            for path in mates:
                path.write_bytes(b"input")
            self.rows.append(["mouse_10", str(lane.resolve()), *map(str, mates), 5, 5])
        self.manifest = self.root / "input_manifest.tsv"
        self.write()

    def write(self):
        with self.manifest.open("w") as stream:
            writer = csv.writer(stream, delimiter="\t")
            writer.writerow(["sample", "lane", "r1", "r2", "r1_bytes", "r2_bytes"])
            writer.writerows(self.rows)

    def test_custom_name_uses_exact_manifest_paths_in_configured_lane_order(self):
        paths = inputs.select(
            self.manifest, "mouse_10", "R1", list(reversed(self.lanes))
        )
        self.assertEqual(paths, [self.rows[1][2], self.rows[0][2]])

    def test_missing_lane_is_rejected(self):
        self.rows.pop()
        self.write()
        with self.assertRaisesRegex(ValueError, "one manifest row per lane"):
            inputs.select(self.manifest, "mouse_10", "R1", self.lanes)

    def test_changed_source_is_rejected(self):
        Path(self.rows[0][2]).write_bytes(b"changed")
        with self.assertRaisesRegex(ValueError, "missing, changed or outside"):
            inputs.select(self.manifest, "mouse_10", "R1", self.lanes)

    def test_duplicate_lanes_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "duplicated"):
            inputs.select(self.manifest, "mouse_10", "R1", [self.lanes[0]] * 2)


if __name__ == "__main__":
    unittest.main()
