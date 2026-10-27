from contextlib import contextmanager
import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "examples"))
from fastq_check import alignment_stats, fastq_stats
from roundtrip import make_inputs


class StreamingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for tool in ("samtools", "pigz", str(ROOT / "scripts/bam_merge")):
            if not shutil.which(tool):
                raise RuntimeError(f"Missing test dependency: {tool}; run make test")

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="cram stream's-")
        self.directory = Path(self.tmp.name)
        make_inputs(self.directory)
        self.sam = self.directory / "aligned.sam"
        self.reference = self.directory / "reference.fa"
        self.f1 = self.directory / "original_R1.fastq"
        self.f2 = self.directory / "original_R2.fastq"
        subprocess.run(["samtools", "faidx", str(self.reference)], check=True)

    def tearDown(self):
        self.tmp.cleanup()

    def cli(self, script, args, data=None):
        return subprocess.run(
            [sys.executable, str(ROOT / "scripts" / script), *map(str, args)],
            input=data, capture_output=True, timeout=30, cwd=self.directory,
        )

    def create(self, prefix, extra=(), data=None):
        return self.cli("create_cram.py", [
            "--input", "-", "--f1", self.f1, "--f2", self.f2,
            "--reference", self.reference, "--output-prefix", self.directory / prefix,
            "--threads", "1", "--quality-mode", "strict", *extra,
        ], data)

    def alignment_bytes(self, path, fmt):
        return subprocess.check_output([
            "samtools", "view", "-h", "-O", fmt, "-T", str(self.reference), str(path),
        ])

    @contextmanager
    def pipes(self, payloads):
        """Independent writers; a reader that reopens or deadlocks hits a timeout."""
        writers = []
        paths = []
        with tempfile.TemporaryDirectory(dir=self.directory) as pipe_dir:
            try:
                for index, payload in enumerate(payloads):
                    source = Path(pipe_dir) / f"source{index}"
                    source.write_bytes(payload)
                    fifo = Path(pipe_dir) / f"input{index}.pipe"
                    os.mkfifo(fifo)
                    paths.append(fifo)
                    writers.append(subprocess.Popen([
                        sys.executable, "-c",
                        "import sys, time\n"
                        "with open(sys.argv[1], 'rb') as src, open(sys.argv[2], 'wb', buffering=0) as dst:\n"
                        "    dst.write(src.read(1))\n"
                        "    time.sleep(0.05)\n"
                        "    dst.write(src.read())\n",
                        str(source), str(fifo),
                    ], stderr=subprocess.PIPE))
                yield paths
                for writer in writers:
                    _, errors = writer.communicate(timeout=10)
                    self.assertEqual(writer.returncode, 0, errors)
            finally:
                for writer in writers:
                    if writer.poll() is None:
                        writer.kill()
                    writer.communicate(timeout=10)

    def assert_no_outputs(self, prefix):
        self.assertEqual(list(self.directory.glob(f"{prefix}.*")), [])
        self.assertEqual(list(self.directory.glob(f".{prefix}.*")), [])

    def test_creation_accepts_sam_bam_cram_stdin(self):
        expected = fastq_stats(self.f1, self.f2, "strict")
        for fmt in ("SAM", "BAM", "CRAM"):
            data = self.alignment_bytes(self.sam, fmt)
            for fast in (False, True):
                with self.subTest(format=fmt, fast=fast):
                    prefix = f"{fmt}_{fast}"
                    extra = ["--fast", "--name-collated"] if fast else []
                    result = self.create(prefix, extra, data)
                    self.assertEqual(result.returncode, 0, result.stderr)
                    cram = self.directory / f"{prefix}.cram"
                    self.assertEqual(alignment_stats(cram, str(self.reference), "strict"), expected)
                    creation = json.loads((self.directory / f"{prefix}.creation.json").read_text())
                    self.assertEqual(creation["verification"], "skipped" if fast else "passed")

    def test_creation_reads_alignment_fifo_once(self):
        with self.pipes([self.sam.read_bytes()]) as (fifo,):
            result = self.create("FIFO", ["--input", fifo])
        self.assertEqual(result.returncode, 0, result.stderr)
        report = json.loads((self.directory / "FIFO.verification.json").read_text())
        self.assertTrue(report["match"])

    def test_streamed_multiple_read_groups_are_rejected(self):
        data = self.sam.read_bytes().replace(b"@RG\t", b"@RG\tID:other\tSM:other\n@RG\t")
        for extra in ([], ["--fast", "--name-collated"]):
            with self.subTest(extra=extra):
                result = self.create("BAD", extra, data)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn(b"one read group", result.stderr)
                self.assert_no_outputs("BAD")

    def test_bad_alignment_stream_is_not_published(self):
        for extra in ([], ["--fast"]):
            with self.subTest(extra=extra):
                result = self.create("BAD", extra, b"not an alignment\n")
                self.assertNotEqual(result.returncode, 0)
                self.assert_no_outputs("BAD")

    def test_fast_creation_reads_separate_fastq_fifos(self):
        expected = fastq_stats(self.f1, self.f2, "strict")
        for compressed in (False, True):
            with self.subTest(compressed=compressed):
                payloads = [p.read_bytes() for p in (self.f1, self.f2)]
                if compressed:
                    payloads = list(map(gzip.compress, payloads))
                prefix = f"PIPED_{compressed}"
                with self.pipes(payloads) as (fifo1, fifo2):
                    result = self.create(prefix, ["--fast", "--f1", fifo1, "--f2", fifo2], self.sam.read_bytes())
                self.assertEqual(result.returncode, 0, result.stderr)
                cram = self.directory / f"{prefix}.cram"
                self.assertEqual(alignment_stats(cram, str(self.reference), "strict"), expected)

    def test_fastq_fifo_requires_fast_before_opening(self):
        fifo = self.directory / "reads.pipe"
        os.mkfifo(fifo)
        result = self.create("BAD", ["--f1", fifo], self.sam.read_bytes())
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(b"FASTQ pipes require --fast", result.stderr)
        self.assert_no_outputs("BAD")

    def test_fastq_checksums_from_separate_plain_and_gzip_fifos(self):
        # More data than a pipe buffer, with split gzip magic bytes.
        for source in (self.f1, self.f2):
            source.write_bytes(source.read_bytes() * 2000)
        expected = fastq_stats(self.f1, self.f2, "strict")
        for compressed in (False, True):
            with self.subTest(compressed=compressed):
                payloads = [p.read_bytes() for p in (self.f1, self.f2)]
                if compressed:
                    payloads = list(map(gzip.compress, payloads))
                output = self.directory / f"stats_{compressed}.json"
                with self.pipes(payloads) as (fifo1, fifo2):
                    result = self.cli("fastq_check.py", [
                        "stats", "--f1", fifo1, "--f2", fifo2,
                        "--quality-mode", "strict", "--output", output,
                    ])
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(json.loads(output.read_text()), expected)

    def test_either_fastq_mate_can_use_stdin(self):
        expected = fastq_stats(self.f1, self.f2, "strict")
        for mate in (1, 2):
            for compressed in (False, True):
                with self.subTest(mate=mate, compressed=compressed):
                    data = (self.f1 if mate == 1 else self.f2).read_bytes()
                    if compressed:
                        data = gzip.compress(data)
                    output = self.directory / "stdin.stats.json"
                    result = self.cli("fastq_check.py", [
                        "stats", "--f1", "-" if mate == 1 else self.f1,
                        "--f2", "-" if mate == 2 else self.f2,
                        "--quality-mode", "strict", "--output", output,
                    ], data)
                    self.assertEqual(result.returncode, 0, result.stderr)
                    self.assertEqual(json.loads(output.read_text()), expected)

    def test_invalid_fastq_stream_preserves_existing_stats(self):
        bad_gzip = bytearray(gzip.compress(self.f1.read_bytes()))
        bad_gzip[-8] ^= 1
        for data in (b"@broken/1\nAC\n+\n", bytes(bad_gzip)):
            with self.subTest(data=data[:10]):
                output = self.directory / "existing.json"
                output.write_bytes(b"previous statistics")
                result = self.cli("fastq_check.py", [
                    "stats", "--f1", "-", "--f2", self.f2, "--output", output,
                ], data)
                self.assertEqual(result.returncode, 2, result.stderr)
                self.assertEqual(output.read_bytes(), b"previous statistics")

    def test_reject_shared_fastq_pipe(self):
        fifo = self.directory / "shared.pipe"
        os.mkfifo(fifo)
        for sources in (("-", "-"), (fifo, fifo)):
            result = self.cli("fastq_check.py", [
                "stats", "--f1", sources[0], "--f2", sources[1], "--output", "bad.json",
            ], b"")
            self.assertEqual(result.returncode, 2, result.stderr)
            self.assertIn(b"separate", result.stderr)
            self.assertFalse((self.directory / "bad.json").exists())
        result = self.create("BAD", ["--fast", "--f1", fifo, "--f2", fifo], self.sam.read_bytes())
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(b"separate input pipes", result.stderr)
        self.assert_no_outputs("BAD")

    def test_check_accepts_alignment_stdin_and_reports_decode_errors(self):
        result = self.create("STORED", data=self.sam.read_bytes())
        self.assertEqual(result.returncode, 0, result.stderr)
        cram = self.directory / "STORED.cram"
        report = self.directory / "comparison.json"
        args = [
            "check", "--input", "-", "--reference", self.reference,
            "--stats", self.directory / "STORED.fastq.stats.json", "--report", report,
        ]
        for fmt in ("SAM", "BAM", "CRAM"):
            with self.subTest(format=fmt):
                result = self.cli("fastq_check.py", args, self.alignment_bytes(cram, fmt))
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertTrue(json.loads(report.read_text())["match"])
        report.write_bytes(b"previous comparison")
        result = self.cli("fastq_check.py", args, cram.read_bytes()[:100])
        self.assertEqual(result.returncode, 2, result.stderr)
        self.assertEqual(report.read_bytes(), b"previous comparison")


if __name__ == "__main__":
    unittest.main()
