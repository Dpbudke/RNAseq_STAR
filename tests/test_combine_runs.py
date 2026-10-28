import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("combine_runs", ROOT / "scripts/combine_runs.py")
combine = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(combine)


class CombineRuns(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="combine_runs_codetestfile_", dir="/tmp")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def run_dir(self, name, samples, genes=("g1", "g2"), column=4, state="validated"):
        run = self.root / name
        (run / "03-Counts").mkdir(parents=True)
        lines = ["\t".join(["gene_id", *samples])]
        for i, gene in enumerate(genes):
            lines.append("\t".join([gene, *(str(i + j + 1) for j in range(len(samples)))]))
        counts = run / "03-Counts" / "counts.txt"
        counts.write_text("\n".join(lines) + "\n")
        sha = hashlib.sha256(counts.read_bytes()).hexdigest()
        (run / "count_validation.json").write_text(
            json.dumps({"status": "passed", "count_column": column, "matrix_sha256": sha}))
        (run / "pipeline_state.txt").write_text(state + "\n")
        return run

    def test_two_pools_join_in_order(self):
        a = self.run_dir("pool_a", ["S_1", "S_2"])
        b = self.run_dir("pool_b", ["S_10", "S_20"])
        result = combine.combine(self.root / "out", [a, b])
        self.assertEqual(result["samples"], 4)
        header = (self.root / "out/counts.txt").read_text().splitlines()[0].split("\t")
        self.assertEqual(header, ["gene_id", "S_1", "S_2", "S_10", "S_20"])
        self.assertTrue((self.root / "out/SHA256SUMS.txt").is_file())

    def test_sample_in_two_runs_is_refused(self):
        a = self.run_dir("pool_a", ["S_1"])
        b = self.run_dir("pool_b", ["S_1"])
        with self.assertRaisesRegex(ValueError, "more than one run"):
            combine.combine(self.root / "out", [a, b])

    def test_gene_order_must_match(self):
        a = self.run_dir("pool_a", ["S_1"])
        b = self.run_dir("pool_b", ["S_2"], genes=("g2", "g1"))
        with self.assertRaisesRegex(ValueError, "gene IDs"):
            combine.combine(self.root / "out", [a, b])

    def test_unfinished_or_changed_run_is_refused(self):
        a = self.run_dir("pool_a", ["S_1"])
        b = self.run_dir("pool_b", ["S_2"], state="running")
        with self.assertRaisesRegex(ValueError, "finalization"):
            combine.combine(self.root / "out", [a, b])
        b = self.run_dir("pool_c", ["S_3"])
        (b / "03-Counts/counts.txt").write_text("gene_id\tS_3\ng1\t9\ng2\t9\n")
        with self.assertRaisesRegex(ValueError, "changed"):
            combine.combine(self.root / "out", [a, b])

    def test_different_count_columns_are_refused(self):
        a = self.run_dir("pool_a", ["S_1"], column=4)
        b = self.run_dir("pool_b", ["S_2"], column=3)
        with self.assertRaisesRegex(ValueError, "count columns"):
            combine.combine(self.root / "out", [a, b])
        self.assertFalse((self.root / "out").exists())


if __name__ == "__main__":
    unittest.main()
