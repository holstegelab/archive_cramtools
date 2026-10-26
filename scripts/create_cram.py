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


def output_paths(prefix):
    """Append extensions without replacing dots within the sample identifier."""
    return {key: prefix.parent / (prefix.name + suffix) for key, suffix in (
        ("cram", ".cram"),
        ("index", ".cram.crai"),
        ("fastq_stats", ".fastq.stats.json"),
        ("verification", ".verification.json"),
        ("merge_stats", ".merge.stats.tsv"),
        ("badmap1", ".badmap_R1.fastq.gz"),
        ("badmap2", ".badmap_R2.fastq.gz"),
        # Publish the completion record last.
        ("creation", ".creation.json"),
    )}


def publish_outputs(staged, destinations):
    """Publish without overwriting existing files, rolling back on failure."""
    published = []
    try:
        for key, destination in destinations.items():
            # Staging is on the same filesystem. Unlike rename/replace, link
            # fails atomically if a concurrent writer created the destination.
            os.link(staged[key], destination)
            published.append(key)
    except BaseException:
        for key in reversed(published):
            try:
                if os.path.samefile(staged[key], destinations[key]):
                    destinations[key].unlink()
            except FileNotFoundError:
                pass
        raise


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
    parser.add_argument("-o", "--output-prefix", required=True, metavar="SAMPLE",
                        help="sample filename prefix, optionally with a directory (e.g. results/SAMPLE001)")
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
        if os.path.basename(args.output_prefix) in ("", ".", ".."):
            raise ValueError("--output-prefix must include a sample filename prefix")
        prefix = Path(args.output_prefix).expanduser().absolute()
        if prefix.is_dir():
            raise ValueError("--output-prefix must include a sample filename prefix, not an existing directory")
        if not prefix.parent.is_dir():
            raise ValueError("Output parent directory does not exist")
        destinations = output_paths(prefix)
        for destination in destinations.values():
            if os.path.lexists(destination):
                raise ValueError(f"Output already exists: {destination}")
        header = subprocess.check_output([args.samtools, "view", "-H", "-T", args.reference, args.input], text=True)
        if sum(line.startswith("@RG\t") for line in header.splitlines()) > 1:
            raise ValueError("Creation accepts one read group at a time; split the input before annotation")
        print("Validating original FASTQs and calculating payload digests...", file=sys.stderr)
        expected = fastq_stats(args.f1, args.f2, args.quality_mode)
        # Stage on the destination filesystem for publication without overwrites.
        # All outputs are verified before any final filename becomes visible.
        with tempfile.TemporaryDirectory(prefix=f".{prefix.name}.", dir=prefix.parent) as stage_name:
            stage = Path(stage_name)
            staged = output_paths(stage / prefix.name)
            with tempfile.TemporaryDirectory(prefix="cram-roundtrip-", dir=args.tmp_dir or stage) as scratch:
                cram = staged["cram"]
                if args.name_collated:
                    commands = [[args.samtools, "view", "-h", "-T", args.reference, args.input]]
                else:
                    commands = [[args.samtools, "collate", "-Ou", "-@", str(args.threads),
                                 "--reference", args.reference, args.input, str(Path(scratch) / "collate")],
                                [args.samtools, "view", "-h", "-"]]
                commands += [[str(SCRIPTS / "bam_merge"), "-a", args.f1, "-b", args.f2,
                              "-ua", str(staged["badmap1"]), "-ub", str(staged["badmap2"]),
                              "-s", str(staged["merge_stats"])],
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
                write_json(staged["fastq_stats"], expected)
                write_json(staged["verification"], report)
                version = subprocess.check_output([args.samtools, "--version"], text=True).splitlines()[0]
                write_json(staged["creation"], {"tool": "cram_roundtrip", "version": "0.1.0",
                           "sample": prefix.name,
                           "samtools": version, "reference_name": Path(args.reference).name,
                           "quality_mode": args.quality_mode, "threads": args.threads})
            publish_outputs(staged, destinations)
        print(f"Created and verified: {destinations['cram']}")
        return 0
    except (OSError, ValueError, RuntimeError, EOFError, subprocess.SubprocessError) as error:
        print(f"create_cram: error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
