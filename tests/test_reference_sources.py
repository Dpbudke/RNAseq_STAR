import gzip
import hashlib
import importlib.util
import io
import json
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


refs = load("fetch_refs_sources", ROOT / "scripts/fetch_references.py")

GENOME = ">chr1 test\nAAAACCCCGGGGTTTT\nACGTACGTACGTACGT\n>chrM\nGGGGAAAACCCCTTTT\n"
GTF = (
    '#!genome-build test\n'
    'chr1\tT\texon\t5\t8\t.\t+\t.\tgene_id "G1"; gene_type "rRNA";\n'
    'chr1\tT\texon\t17\t20\t.\t-\t.\tgene_id "G2"; gene_biotype "rRNA";\n'
    'chr1\tT\texon\t1\t4\t.\t+\t.\tgene_id "G3"; gene_type "protein_coding";\n'
    'chrM\tT\texon\t1\t4\t.\t+\t.\tgene_id "G4"; gene_type "Mt_rRNA";\n'
)


def gz(data):
    return gzip.compress(data.encode(), mtime=0)


class MouseUnchanged(unittest.TestCase):
    """The GENCODE mouse path writes the same files and manifest bytes as the version before presets."""

    def run_fetch(self, script, directory, rrna, files):
        md5s = "".join(f"{hashlib.md5(data).hexdigest()}  {name}\n" for name, data in files.items())

        def urlopen(url, timeout=None):
            name = url.rsplit("/", 1)[1]
            return io.BytesIO(md5s.encode() if name == "MD5SUMS" else files[name])

        module = load(f"fetch_{script.stem}_{directory.name}", script)
        args = ["fetch", "--directory", str(directory), "--release", "M35", "--rrna", str(rrna)]
        with patch("sys.argv", args), patch.object(module.urllib.request, "urlopen", side_effect=urlopen):
            module.main()
        return {p.name: p.read_bytes() for p in sorted(directory.iterdir())}

    def test_same_outputs_as_previous_version(self):
        old = subprocess.run(["git", "-C", str(ROOT), "show", "bennettlab/main:scripts/fetch_references.py"],
                             capture_output=True, text=True)
        if old.returncode != 0:
            self.skipTest("previous version not available in this checkout")
        files = {"GRCm39.primary_assembly.genome.fa.gz": gz(GENOME),
                 "gencode.vM35.annotation.gtf.gz": gz(GTF)}
        with tempfile.TemporaryDirectory(prefix="refsources_codetestfile_") as d:
            root = Path(d)
            (root / "old_fetch.py").write_text(old.stdout)
            rrna = root / "mouse_rrna.fasta"
            rrna.write_text(">r\nACGT\n")
            before = self.run_fetch(root / "old_fetch.py", root / "old", rrna, files)
            after = self.run_fetch(ROOT / "scripts/fetch_references.py", root / "new", rrna, files)
            self.assertEqual(before, after)
            # Re-running on the folder the previous version wrote verifies it and changes nothing.
            again = self.run_fetch(ROOT / "scripts/fetch_references.py", root / "old", rrna, files)
            self.assertEqual(before, again)


class OtherSources(unittest.TestCase):
    def test_rrna_from_annotation(self):
        with tempfile.TemporaryDirectory(prefix="refsources_codetestfile_") as d:
            root = Path(d)
            (root / "genome.fa").write_text(GENOME)
            (root / "annotation.gtf").write_text(GTF)
            n = refs.rrna_from_annotation(root / "genome.fa", root / "annotation.gtf", root / "rrna.fasta")
            self.assertEqual(n, 3)
            text = (root / "rrna.fasta").read_text()
            self.assertIn(">G1|chr1:5-8(+)\nCCCC\n", text)
            self.assertIn(">G2|chr1:17-20(-)\nACGT\n", text)  # reverse complement of ACGT
            self.assertIn(">G4|chrM:1-4(+)\nGGGG\n", text)
            self.assertNotIn("G3", text)

    def test_no_rrna_genes_is_refused(self):
        with tempfile.TemporaryDirectory(prefix="refsources_codetestfile_") as d:
            root = Path(d)
            (root / "genome.fa").write_text(GENOME)
            (root / "annotation.gtf").write_text('chr1\tT\texon\t1\t4\t.\t+\t.\tgene_id "G3"; gene_type "lncRNA";\n')
            with self.assertRaises(ValueError):
                refs.rrna_from_annotation(root / "genome.fa", root / "annotation.gtf", root / "rrna.fasta")
            self.assertFalse((root / "rrna.fasta").exists())

    def test_contig_names_must_overlap(self):
        with tempfile.TemporaryDirectory(prefix="refsources_codetestfile_") as d:
            root = Path(d)
            (root / "genome.fa").write_text(">1\nACGT\n")
            (root / "annotation.gtf").write_text(GTF)
            with self.assertRaises(ValueError):
                refs.check_contigs(root / "genome.fa", root / "annotation.gtf")
            (root / "genome.fa").write_text(GENOME)
            refs.check_contigs(root / "genome.fa", root / "annotation.gtf")

    def test_ensembl_picks_primary_assembly_and_main_gtf(self):
        dna = ("1 2 Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz\n"
               "3 4 Rattus_norvegicus.mRatBN7.2.dna.primary_assembly.fa.gz\n"
               "5 6 Rattus_norvegicus.mRatBN7.2.dna_sm.primary_assembly.fa.gz\n")
        gtf = ("7 8 Rattus_norvegicus.mRatBN7.2.113.gtf.gz\n"
               "9 10 Rattus_norvegicus.mRatBN7.2.113.chr.gtf.gz\n"
               "11 12 Rattus_norvegicus.mRatBN7.2.113.abinitio.gtf.gz\n")

        class A:
            source, species, release = "ensembl", "rattus_norvegicus", "113"

        with patch.object(refs, "fetch_text", side_effect=lambda url: dna if "/fasta/" in url else gtf):
            ident, files, _ = refs.plan(A)
        self.assertEqual(ident["assembly"], "mRatBN7")
        self.assertTrue(files["genome.fa"]["source_url"].endswith(".dna.primary_assembly.fa.gz"))
        self.assertNotIn("dna_sm", files["genome.fa"]["source_url"])
        self.assertEqual(files["genome.fa"]["source_bsd_sum"], [3, 4])
        self.assertTrue(files["annotation.gtf"]["source_url"].endswith("mRatBN7.2.113.gtf.gz"))

    def test_custom_reference_and_refusal_of_another(self):
        with tempfile.TemporaryDirectory(prefix="refsources_codetestfile_") as d:
            root = Path(d)
            (root / "my.fa.gz").write_bytes(gz(GENOME))
            (root / "my.gtf").write_text(GTF)
            ref = root / "ref"
            args = ["fetch", "--directory", str(ref), "--source", "custom", "--name", "test",
                    "--genome-fa", str(root / "my.fa.gz"), "--gtf", str(root / "my.gtf")]
            with patch("sys.argv", args):
                refs.main()
            manifest = json.loads((ref / "reference_manifest.json").read_text())
            self.assertEqual(manifest["source"], "custom")
            self.assertEqual((ref / "genome.fa").read_text(), GENOME)
            self.assertEqual(manifest["files"]["rrna.fasta"]["records"], 3)
            with patch("sys.argv", args):
                refs.main()  # same reference again: verified, unchanged
            other = args[:6] + ["other"] + args[7:]
            with patch("sys.argv", other), self.assertRaises(ValueError):
                refs.main()


if __name__ == "__main__":
    unittest.main()
