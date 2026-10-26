#!/usr/bin/env python3
"""Store an alignment and its original paired FASTQs in a verified tagged CRAM."""
import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

from fastq_check import alignment_stats, comparison, fastq_stats, write_json

SCRIPTS = Path(__file__).resolve().parent


def pipeline(commands):
    processes = []
    try:
        previous = None
        for index, command in enumerate(commands):
            process = subprocess.Popen(command, stdin=previous,
                                       stdout=subprocess.PIPE if index < len(commands) - 1 else None)
            if previous is not None:
                previous.close()
            processes.append(process)
            previous = process.stdout
        failures = []
        for command, process in reversed(list(zip(commands, processes))):
            status = process.wait()
            if status:
                failures.append(f"{Path(str(command[0])).name} exited with status {status}")
        if failures:
            raise RuntimeError("; ".join(failures))
    finally:
        for process in processes:
            if process.poll() is None:
                process.terminate()
        for process in processes:
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input", required=True, help="SAM/BAM/CRAM from one read group, before tag annotation")
    parser.add_argument("--f1", required=True, help="original untrimmed R1 FASTQ[.gz], Phred+33")
    parser.add_argument("--f2", required=True, help="original untrimmed R2 FASTQ[.gz], Phred+33")
    parser.add_argument("-T", "--reference", required=True, help="matching indexed FASTA")
    parser.add_argument("-o", "--output-dir", required=True, help="new directory for CRAM, index and verification reports")
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--sort-memory", default="768M", help="samtools sort memory per thread")
    parser.add_argument("--tmp-dir", help="scratch directory for sorting/collation")
    parser.add_argument("--name-collated", action="store_true", help="input already groups QNAMEs; preserves aligner order and reduces buffering")
    parser.add_argument("--quality-mode", choices=("pipeline", "strict"), default="pipeline")
    parser.add_argument("--samtools", default="samtools")
    args = parser.parse_args(argv)
    try:
        if args.threads < 1:
            raise ValueError("--threads must be positive")
        for key in ("input", "f1", "f2", "reference"):
            value = Path(getattr(args, key)).resolve(strict=True)
            if not value.is_file():
                raise ValueError(f"--{key} must be a regular file")
            setattr(args, key, str(value))
        if not Path(args.reference + ".fai").is_file():
            raise ValueError("Index the reference first: samtools faidx reference.fa")
        for tool in (args.samtools, "pigz", str(SCRIPTS / "bam_merge")):
            if not shutil.which(tool):
                raise ValueError(f"Executable missing: {tool}; run make and check PATH")
        destination = Path(args.output_dir).absolute()
        if os.path.lexists(destination):
            raise ValueError(f"Output directory already exists: {destination}")
        if not destination.parent.is_dir():
            raise ValueError("Output parent directory does not exist")
        header = subprocess.check_output([args.samtools, "view", "-H", "-T", args.reference, args.input], text=True)
        if sum(line.startswith("@RG\t") for line in header.splitlines()) > 1:
            raise ValueError("Creation accepts one read group at a time; split the input before annotation")
        print("Validating original FASTQs and calculating payload digests...", file=sys.stderr)
        expected = fastq_stats(args.f1, args.f2, args.quality_mode)
        # Staging beside the destination allows a final directory rename. Scratch
        # can be elsewhere; failures leave the destination unpublished.
        with tempfile.TemporaryDirectory(prefix=f".{destination.name}.", dir=destination.parent) as stage_name:
            stage = Path(stage_name)
            with tempfile.TemporaryDirectory(prefix="cram-roundtrip-", dir=args.tmp_dir or stage) as scratch:
                cram = stage / "reads.cram"
                if args.name_collated:
                    commands = [[args.samtools, "view", "-h", "-T", args.reference, args.input]]
                else:
                    commands = [[args.samtools, "collate", "-Ou", "-@", str(args.threads),
                                 "--reference", args.reference, args.input, str(Path(scratch) / "collate")],
                                [args.samtools, "view", "-h", "-"]]
                commands += [[str(SCRIPTS / "bam_merge"), "-a", args.f1, "-b", args.f2,
                              "-ua", str(stage / "badmap_R1.fastq.gz"), "-ub", str(stage / "badmap_R2.fastq.gz"),
                              "-s", str(stage / "merge.stats.tsv")],
                             [args.samtools, "sort", "-@", str(args.threads), "-m", args.sort_memory,
                              "-T", str(Path(scratch) / "sort"), "--reference", args.reference,
                              "-O", "CRAM", "-o", str(cram), "-"]]
                pipeline(commands)
                subprocess.run([args.samtools, "quickcheck", "-v", str(cram)], check=True)
                print("Decoding the completed CRAM and verifying restored reads...", file=sys.stderr)
                observed = alignment_stats(cram, args.reference, args.quality_mode, args.samtools, args.threads)
                report = comparison(expected, observed)
                if not report["match"]:
                    raise ValueError("Restored CRAM differs from source FASTQs: " + ", ".join(report["differences"]))
                subprocess.run([args.samtools, "index", "-@", str(args.threads), str(cram)], check=True)
                write_json(stage / "fastq.stats.json", expected)
                write_json(stage / "verification.json", report)
                version = subprocess.check_output([args.samtools, "--version"], text=True).splitlines()[0]
                write_json(stage / "creation.json", {"tool": "cram_roundtrip", "version": "0.1.0",
                           "samtools": version, "reference_name": Path(args.reference).name,
                           "quality_mode": args.quality_mode, "threads": args.threads})
            if os.path.lexists(destination):
                raise ValueError("Output directory appeared during processing")
            os.rename(stage, destination)
        print(f"Created and verified: {destination / 'reads.cram'}")
        return 0
    except (OSError, ValueError, RuntimeError, EOFError, subprocess.SubprocessError) as error:
        print(f"create_cram: error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
