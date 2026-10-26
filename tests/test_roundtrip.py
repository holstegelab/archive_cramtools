from contextlib import redirect_stderr, redirect_stdout
import gzip
import io
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "examples"))
from fastq_check import alignment_stats, comparison, fastq_stats, load_stats
from roundtrip import make_inputs
from create_cram import publish_outputs
import create_cram


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

    def create(self, extra=(), success=True, prefix=None):
        if prefix is None:
            prefix = self.directory / "SAMPLE001"
        cmd = [sys.executable, str(ROOT / "scripts/create_cram.py"), "--input", str(self.directory / "aligned.sam"),
               "--f1", str(self.f1), "--f2", str(self.f2), "--reference", str(self.reference),
               "--output-prefix", str(prefix), "--threads", "1", "--quality-mode", "strict", *extra]
        result = subprocess.run(cmd, cwd=self.directory, capture_output=True, text=True)
        if success:
            self.assertEqual(result.returncode, 0, result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0, result.stderr)
        return result

    def test_create_cram_revert_and_compare_payloads(self):
        self.create()
        cram = self.directory / "SAMPLE001.cram"
        self.assertTrue(Path(str(cram) + ".crai").is_file())
        report = json.loads((self.directory / "SAMPLE001.verification.json").read_text())
        self.assertTrue(report["match"])
        creation = json.loads((self.directory / "SAMPLE001.creation.json").read_text())
        self.assertEqual(creation["verification"], "passed")
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

    def test_fast_mode_skips_payload_scans_and_can_be_verified_later(self):
        # Run the real encoder but fail if either expensive checksum pass starts.
        stdout, stderr = io.StringIO(), io.StringIO()
        with mock.patch.object(create_cram, "fastq_stats", side_effect=AssertionError("FASTQ scan")), \
             mock.patch.object(create_cram, "alignment_stats", side_effect=AssertionError("CRAM scan")), \
             redirect_stdout(stdout), redirect_stderr(stderr):
            status = create_cram.main([
                "--input", str(self.directory / "aligned.sam"),
                "--f1", str(self.f1), "--f2", str(self.f2),
                "--reference", str(self.reference), "--output-prefix", str(self.directory / "FAST"),
                "--threads", "1", "--fast",
            ])
        self.assertEqual(status, 0, stderr.getvalue())
        self.assertIn("read verification skipped", stdout.getvalue())
        self.assertNotIn("Created and verified", stdout.getvalue())
        self.assertEqual({p.name for p in self.directory.glob("FAST.*")}, {
            "FAST.cram", "FAST.cram.crai", "FAST.merge.stats.tsv",
            "FAST.badmap_R1.fastq.gz", "FAST.badmap_R2.fastq.gz", "FAST.creation.json",
        })
        creation = json.loads((self.directory / "FAST.creation.json").read_text())
        self.assertEqual(creation["verification"], "skipped")
        self.assertEqual(list(self.directory.glob(".FAST.*")), [])
        # Independent verification includes the pair dropped from the alignment.
        self.assertEqual(
            fastq_stats(self.f1, self.f2, "strict"),
            alignment_stats(self.directory / "FAST.cram", str(self.reference), "strict"),
        )

    def test_fast_mode_defers_payload_mismatch_to_later_check(self):
        self.f1.write_text(self.f1.read_text().replace("pairA", "different"))
        self.f2.write_text(self.f2.read_text().replace("pairA", "different"))
        self.create(["--fast"])
        stats = self.directory / "expected.json"
        result = subprocess.run([
            sys.executable, str(ROOT / "scripts/fastq_check.py"), "stats",
            "--f1", str(self.f1), "--f2", str(self.f2), "--quality-mode", "strict",
            "--output", str(stats),
        ], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        result = subprocess.run([
            sys.executable, str(ROOT / "scripts/fastq_check.py"), "check",
            "--stats", str(stats), "--input", str(self.directory / "SAMPLE001.cram"),
            "--reference", str(self.reference),
        ], capture_output=True, text=True)
        self.assertEqual(result.returncode, 1, result.stderr)
        self.assertIn("Mismatch", result.stdout)

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
        self.assertEqual(list(self.directory.glob("SAMPLE001.*")), [])
        self.assertEqual(list(self.directory.glob(".SAMPLE001.*")), [])

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
                                 "-i", str(self.directory / "SAMPLE001.cram"), "-T", str(self.reference)], capture_output=True, text=True)
        self.assertEqual(result.returncode, 1, result.stderr)
        self.assertIn("Mismatch", result.stdout)

    def test_existing_sidecar_is_untouched(self):
        sentinel = self.directory / "SAMPLE001.verification.json"
        sentinel.write_text("unchanged")
        for extra in ([], ["--fast"]):
            with self.subTest(extra=extra):
                self.create(extra, success=False)
                self.assertEqual(sentinel.read_text(), "unchanged")
                self.assertEqual(list(self.directory.glob("SAMPLE001.*")), [sentinel])

    def test_sample_prefix_is_used_for_every_output(self):
        self.create(prefix="SAMPLE.v1")
        expected = {
            "SAMPLE.v1.cram", "SAMPLE.v1.cram.crai", "SAMPLE.v1.fastq.stats.json",
            "SAMPLE.v1.verification.json", "SAMPLE.v1.merge.stats.tsv",
            "SAMPLE.v1.badmap_R1.fastq.gz", "SAMPLE.v1.badmap_R2.fastq.gz",
            "SAMPLE.v1.creation.json",
        }
        self.assertEqual({p.name for p in self.directory.glob("SAMPLE.v1.*")}, expected)
        creation = json.loads((self.directory / "SAMPLE.v1.creation.json").read_text())
        self.assertEqual(creation["sample"], "SAMPLE.v1")
        self.assertEqual(list(self.directory.glob(".SAMPLE.v1.*")), [])

    def test_two_samples_share_an_existing_results_directory(self):
        results = self.directory / "results"
        results.mkdir()
        self.create(prefix=results / "sampleA")
        first = {p.name: p.read_bytes() for p in results.iterdir()}
        self.create(prefix=results / "sampleB")
        self.assertEqual(len(list(results.iterdir())), 16)
        for name, data in first.items():
            self.assertEqual((results / name).read_bytes(), data)
        self.assertTrue((results / "sampleB.creation.json").is_file())

    def test_broken_symlink_output_is_not_overwritten(self):
        target = self.directory / "SAMPLE001.cram.crai"
        target.symlink_to(self.directory / "missing")
        self.create(success=False)
        self.assertTrue(target.is_symlink())
        self.assertFalse((self.directory / "SAMPLE001.cram").exists())

    def test_directory_is_rejected_as_a_sample_prefix(self):
        for prefix in (".", str(self.directory) + "/", str(self.directory)):
            with self.subTest(prefix=prefix):
                result = self.create(prefix=prefix, success=False)
                self.assertIn("sample filename prefix", result.stderr)

    def test_publication_collision_rolls_back_new_files(self):
        stage = self.directory / "stage"
        stage.mkdir()
        staged = {"cram": stage / "sample.cram", "index": stage / "sample.cram.crai"}
        destinations = {key: self.directory / path.name for key, path in staged.items()}
        for path in staged.values():
            path.write_bytes(b"new data")
        destinations["index"].write_bytes(b"concurrent writer")
        with self.assertRaises(FileExistsError):
            publish_outputs(staged, destinations)
        self.assertFalse(destinations["cram"].exists())
        self.assertEqual(destinations["index"].read_bytes(), b"concurrent writer")
        self.assertTrue(all(path.is_file() for path in staged.values()))

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
        self.assertEqual(list(self.directory.glob("SAMPLE001.*")), [])
        self.assertEqual(list(self.directory.glob(".SAMPLE001.*")), [])

    def test_corrupt_gzip_is_rejected(self):
        broken = self.directory / "broken.fastq.gz"
        data = bytearray(gzip.compress(self.f1.read_bytes()))
        data[-8] ^= 1
        broken.write_bytes(data)
        self.f1 = broken
        for extra in ([], ["--fast"]):
            with self.subTest(extra=extra):
                self.create(extra, success=False)
                self.assertEqual(list(self.directory.glob("SAMPLE001.*")), [])
                self.assertEqual(list(self.directory.glob(".SAMPLE001.*")), [])

    def test_truncated_cram_is_an_error(self):
        self.create()
        cram = self.directory / "SAMPLE001.cram"
        cram.write_bytes(cram.read_bytes()[:100])
        result = subprocess.run([sys.executable, str(ROOT / "scripts/fastq_check.py"), "check",
                                 "--stats", str(self.directory / "SAMPLE001.fastq.stats.json"),
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
