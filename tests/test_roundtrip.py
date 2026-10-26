import gzip
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "examples"))
from fastq_check import alignment_stats, comparison, fastq_stats, load_stats
from roundtrip import make_inputs


class RoundtripTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for tool in ("samtools", "pigz", str(ROOT / "scripts/bam_merge")):
            if not shutil.which(tool):
                raise RuntimeError(f"Missing test dependency: {tool}; run make test")

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="cram test's-")
        self.directory = Path(self.tmp.name)
        make_inputs(self.directory)
        self.reference = self.directory / "reference.fa"
        self.f1 = self.directory / "original_R1.fastq"
        self.f2 = self.directory / "original_R2.fastq"
        subprocess.run(["samtools", "faidx", str(self.reference)], check=True)

    def tearDown(self):
        self.tmp.cleanup()

    def create(self, extra=(), success=True):
        cmd = [sys.executable, str(ROOT / "scripts/create_cram.py"), "--input", str(self.directory / "aligned.sam"),
               "--f1", str(self.f1), "--f2", str(self.f2), "--reference", str(self.reference),
               "--output-dir", str(self.directory / "stored"), "--threads", "1", "--quality-mode", "strict", *extra]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if success:
            self.assertEqual(result.returncode, 0, result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0, result.stderr)
        return result

    def test_create_cram_revert_and_compare_payloads(self):
        self.create()
        cram = self.directory / "stored/reads.cram"
        self.assertTrue(Path(str(cram) + ".crai").is_file())
        report = json.loads((self.directory / "stored/verification.json").read_text())
        self.assertTrue(report["match"])
        self.assertEqual(report["observed"]["mates"][0]["reads"], 3)
        restored = [self.directory / f"restored_R{i}.fastq.gz" for i in (1, 2)]
        result = subprocess.run([sys.executable, str(ROOT / "scripts/bam_revert.py"), "-i", str(cram),
                                 "-T", str(self.reference), "--f1", str(restored[0]), "--f2", str(restored[1]),
                                 "-s", str(self.directory / "revert.tsv"), "--threads", "1"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(fastq_stats(*restored, "strict"), fastq_stats(self.f1, self.f2, "strict"))
        sam = subprocess.check_output(["samtools", "view", "-T", str(self.reference), str(cram)], text=True)
        self.assertIn("ZB:Z:TT", sam)
        self.assertIn("YB:Z:CC", sam)
        self.assertIn("filtered\t589\t", sam)

    def test_gzip_sources_and_name_collated_input(self):
        for key in ("f1", "f2"):
            path = getattr(self, key)
            output = path.with_suffix(".fastq.gz")
            with gzip.open(output, "wb") as handle:
                handle.write(path.read_bytes())
            setattr(self, key, output)
        self.create(["--name-collated"])

    def test_reject_truncated_fastq_without_publishing(self):
        self.f1.write_text(self.f1.read_text().rsplit("\n", 2)[0] + "\n")
        self.create(success=False)
        self.assertFalse((self.directory / "stored").exists())

    def test_reject_mismatched_mate_counts(self):
        self.f2.write_text(self.f2.read_text() + "@extra/2\nAC\n+\nII\n")
        with self.assertRaisesRegex(ValueError, "different numbers"):
            fastq_stats(self.f1, self.f2)

    def test_quality_normalization_is_explicit(self):
        before = fastq_stats(self.f1, self.f2, "pipeline")
        strict = fastq_stats(self.f1, self.f2, "strict")
        self.f1.write_text(self.f1.read_text().replace("IHHG##", "IHHG!!"))
        self.assertEqual(before, fastq_stats(self.f1, self.f2, "pipeline"))
        self.assertNotEqual(strict, fastq_stats(self.f1, self.f2, "strict"))

    def test_name_changes_are_detected_and_cli_returns_failure(self):
        self.create()
        self.f1.write_text(self.f1.read_text().replace("pairA", "renamed"))
        self.f2.write_text(self.f2.read_text().replace("pairA", "renamed"))
        stats = self.directory / "expected.json"
        stats.write_text(json.dumps(fastq_stats(self.f1, self.f2, "strict")))
        result = subprocess.run([sys.executable, str(ROOT / "scripts/fastq_check.py"), "check", "--stats", str(stats),
                                 "-i", str(self.directory / "stored/reads.cram"), "-T", str(self.reference)], capture_output=True, text=True)
        self.assertEqual(result.returncode, 1, result.stderr)
        self.assertIn("Mismatch", result.stdout)

    def test_existing_output_directory_is_untouched(self):
        output = self.directory / "stored"
        output.mkdir()
        sentinel = output / "keep"
        sentinel.write_text("unchanged")
        self.create(success=False)
        self.assertEqual(sentinel.read_text(), "unchanged")

    def test_missing_restore_tag_fails(self):
        sam = self.directory / "broken.sam"
        sam.write_text("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:200\n"
                       "broken\t65\tchr1\t1\t60\t4M2H\t*\t0\t0\tACGT\tIIII\n")
        with self.assertRaisesRegex(RuntimeError, "hard clip"):
            alignment_stats(sam, self.reference)

    def test_optional_clipping_retains_the_read_payload(self):
        self.f1.write_text("@clip/1\n" + "A" * 40 + "TTGG\n+\n" + "I" * 44 + "\n")
        self.f2.write_text("@clip/2\nCCGG\n+\nIIII\n")
        sam = ("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:200\n"
               "clip\t73\tchr1\t1\t60\t40M4S\t*\t0\t0\t" + "A" * 40 + "TTGG\t" + "I" * 44 + "\n"
               "clip\t133\t*\t0\t0\t*\tchr1\t1\t0\tCCGG\tIIII\n")
        clipped = self.directory / "clipped.sam"
        result = subprocess.run([str(ROOT / "scripts/bam_dechimer"), "-i", "-", "-o", str(clipped),
                                 "-s", str(self.directory / "clipping.tsv")], input=sam, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("ZB:Z:TTGG", clipped.read_text())
        self.assertIn("40M4H", clipped.read_text())
        self.assertEqual(fastq_stats(self.f1, self.f2, "strict"), alignment_stats(clipped, self.reference, "strict"))

    def test_final_payload_mismatch_prevents_publication(self):
        self.f1.write_text(self.f1.read_text().replace("pairA", "different"))
        self.f2.write_text(self.f2.read_text().replace("pairA", "different"))
        result = self.create(success=False)
        self.assertIn("Restored CRAM differs", result.stderr)
        self.assertFalse((self.directory / "stored").exists())

    def test_corrupt_gzip_is_rejected(self):
        broken = self.directory / "broken.fastq.gz"
        data = bytearray(gzip.compress(self.f1.read_bytes()))
        data[-8] ^= 1
        broken.write_bytes(data)
        self.f1 = broken
        self.create(success=False)
        self.assertFalse((self.directory / "stored").exists())

    def test_truncated_cram_is_an_error(self):
        self.create()
        cram = self.directory / "stored/reads.cram"
        cram.write_bytes(cram.read_bytes()[:100])
        result = subprocess.run([sys.executable, str(ROOT / "scripts/fastq_check.py"), "check",
                                 "--stats", str(self.directory / "stored/fastq.stats.json"),
                                 "-i", str(cram), "-T", str(self.reference)], capture_output=True, text=True)
        self.assertEqual(result.returncode, 2, result.stderr)
        self.assertNotIn("Match", result.stdout)

    def test_invalid_statistics_are_rejected(self):
        path = self.directory / "bad.json"
        path.write_text('{}')
        with self.assertRaises(ValueError):
            load_stats(path)

    def test_duplicate_multiplicity_is_not_cancelled(self):
        before = fastq_stats(self.f1, self.f2, "strict")
        for path in (self.f1, self.f2):
            lines = path.read_text().splitlines(True)
            # Replace two copies with two different copies, keeping read counts.
            path.write_text("".join(lines[:4] * 2 + lines[8:]))
        after = fastq_stats(self.f1, self.f2, "strict")
        self.assertFalse(comparison(before, after)["match"])


if __name__ == "__main__":
    unittest.main()
