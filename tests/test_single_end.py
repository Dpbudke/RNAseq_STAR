import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


manifest = load("sample_manifest_se", ROOT / "scripts/sample_manifest.py")
inputs = load("manifest_inputs_se", ROOT / "scripts/manifest_inputs.py")
summary = load("hts_summary_se", ROOT / "2_HTStream_Preprocessing/3_hts_summary.py")

REGEX = "s/_S[0-9]+_L[0-9]+$//"


class SingleEnd(unittest.TestCase):
    def lane(self, root, name, mates):
        lane = root / name
        lane.mkdir()
        for sample in ("A1", "B2"):
            for mate in mates:
                (lane / f"{sample}_S1_L001_{mate}_001.fastq.gz").write_bytes(b"x")
        return lane

    def test_single_end_collects_r1_only(self):
        with tempfile.TemporaryDirectory(prefix="single_end_codetestfile_") as d:
            lane = self.lane(Path(d), "lane1", ["R1"])
            samples, records = manifest.collect([str(lane)], REGEX, ["_R1_001.fastq.gz"], set())
            self.assertEqual(samples, ["A1", "B2"])
            self.assertTrue(all(r[3] == "" and r[5] == 0 for r in records))

    def test_paired_still_requires_r2(self):
        with tempfile.TemporaryDirectory(prefix="single_end_codetestfile_") as d:
            lane = self.lane(Path(d), "lane1", ["R1"])
            with self.assertRaises(ValueError):
                manifest.collect([str(lane)], REGEX, ["_R1_001.fastq.gz", "_R2_001.fastq.gz"], set())

    def test_no_r2_from_a_single_end_manifest(self):
        with tempfile.TemporaryDirectory(prefix="single_end_codetestfile_") as d:
            lane = self.lane(Path(d), "lane1", ["R1"])
            _, records = manifest.collect([str(lane)], REGEX, ["_R1_001.fastq.gz"], set())
            table = Path(d) / "input_manifest.tsv"
            with table.open("w", newline="") as f:
                w = csv.writer(f, delimiter="\t", lineterminator="\n")
                w.writerow(["sample", "lane", "r1", "r2", "r1_bytes", "r2_bytes"])
                w.writerows(records)
            self.assertEqual(len(inputs.select(str(table), "A1", "R1", [str(lane)])), 1)
            with self.assertRaises(ValueError):
                inputs.select(str(table), "A1", "R2", [str(lane)])

    def test_summary_without_duplicate_stage(self):
        import json
        def stage(note, **frag):
            return {"Program_details": {"options": {"notes": note}}, "Fragment": frag,
                    "Single_end": {"hits": 3}}
        entries = [stage("initial stats", **{"in": 100, "basepairs_in": 1000}), stage("screen phix"),
                   stage("count the number of rRNA reads"),
                   stage("trim adapters", basepairs_in=1000, basepairs_out=900),
                   stage("final stats", out=80, basepairs_out=800)]
        with tempfile.TemporaryDirectory(prefix="single_end_codetestfile_") as d:
            path = Path(d) / "S.json"
            path.write_text(json.dumps(entries))
            row = summary.summarise(str(path))
        self.assertEqual(row["dups"], 0)
        self.assertEqual(row["rrna"], 3)

    def test_screened_hits_for_both_layouts(self):
        self.assertEqual(summary.hits({"Paired_end": {"hits": 7}, "Single_end": {"hits": 0}}), 7)
        self.assertEqual(summary.hits({"Single_end": {"hits": 5}}), 5)
        self.assertEqual(summary.hits({}), 0)


if __name__ == "__main__":
    unittest.main()
