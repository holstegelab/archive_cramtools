#!/usr/bin/env python3
"""Run a complete round trip using only tiny, synthetic records."""
import argparse
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from fastq_check import fastq_stats


def make_inputs(directory):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "reference.fa").write_text(">chr1\n" + "A" * 200 + "\n", encoding="ascii")
    records = [("pairA", "ACGTTT", "IHHG##", "AACCGG", "IIII##"),
               ("pairB", "TGCA", "IJ!#", "NACC", "#IHH"),
               ("filtered", "CCGG", "HHHH", "AATT", "IIII")]
    for mate in (1, 2):
        with (directory / f"original_R{mate}.fastq").open("w", encoding="ascii") as handle:
            for name, s1, q1, s2, q2 in records:
                seq, qual = (s1, q1) if mate == 1 else (s2, q2)
                handle.write(f"@{name}/{mate}\n{seq}\n+\n{qual}\n")
    # Coordinate sorted; two primary pairs plus supplementary/secondary hits.
    # The third FASTQ pair was dropped by an upstream filter and is re-added.
    (directory / "aligned.sam").write_text(
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:200\n@RG\tID:demo\tSM:synthetic\n"
        "pairA\t99\tchr1\t5\t60\t4M\t=\t100\t99\tACGT\tIHHG\tRG:Z:demo\tAS:i:4\n"
        "pairA\t2113\tchr1\t20\t30\t2M2H\t=\t100\t0\tAC\tIH\tRG:Z:demo\tAS:i:2\n"
        "pairA\t321\tchr1\t30\t20\t4M\t=\t100\t0\tACGT\tIHHG\tRG:Z:demo\tAS:i:3\n"
        "pairA\t147\tchr1\t100\t60\t4M\t=\t5\t-99\tGGTT\tIIII\tRG:Z:demo\tAS:i:4\n"
        "pairB\t77\t*\t0\t0\t*\t*\t0\t0\tTG\tIJ\tRG:Z:demo\n"
        "pairB\t141\t*\t0\t0\t*\t*\t0\t0\tNA\t#I\tRG:Z:demo\n",
        encoding="ascii")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", default="demo-output", help="new directory for synthetic inputs and demo results")
    args = parser.parse_args()
    directory = Path(args.work_dir).resolve()
    directory.mkdir(parents=True, exist_ok=False)
    make_inputs(directory)
    subprocess.run(["samtools", "faidx", str(directory / "reference.fa")], check=True)
    subprocess.run([sys.executable, str(ROOT / "scripts/create_cram.py"),
                    "--input", str(directory / "aligned.sam"), "--f1", str(directory / "original_R1.fastq"),
                    "--f2", str(directory / "original_R2.fastq"), "--reference", str(directory / "reference.fa"),
                    "--quality-mode", "strict", "--threads", "1", "--output-prefix", str(directory / "DEMO")], check=True)
    subprocess.run([sys.executable, str(ROOT / "scripts/bam_revert.py"),
                    "--input", str(directory / "DEMO.cram"), "--reference", str(directory / "reference.fa"),
                    "--f1", str(directory / "DEMO.restored_R1.fastq.gz"), "--f2", str(directory / "DEMO.restored_R2.fastq.gz"),
                    "--stats", str(directory / "DEMO.revert.stats.tsv"), "--threads", "1"], check=True)
    expected = fastq_stats(directory / "original_R1.fastq", directory / "original_R2.fastq", "strict")
    restored = fastq_stats(directory / "DEMO.restored_R1.fastq.gz", directory / "DEMO.restored_R2.fastq.gz", "strict")
    if restored != expected:
        raise RuntimeError("Restored FASTQs differ from originals")
    print(f"PASS: all 3 synthetic pairs restored, including qualities. Files: {directory}")


if __name__ == "__main__":
    main()
