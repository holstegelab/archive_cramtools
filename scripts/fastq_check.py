#!/usr/bin/env python3
"""Validate FASTQ payloads against restored SAM/BAM/CRAM primary reads."""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile
from itertools import zip_longest
from types import SimpleNamespace

from bam_revert import BamRevertError, SamtoolsSamReader, parse_record

SCHEMA = "cram-roundtrip-sha256-sum-v1"
MODULUS = 1 << 256


def read_fastq(path):
    """Read complete four-line FASTQ records, checking gzip CRCs through EOF."""
    with open(path, "rb") as probe:
        compressed = probe.read(2) == b"\x1f\x8b"
    opener = gzip.open if compressed else open
    with opener(path, "rt", encoding="ascii", newline=None) as handle:
        number = 0
        while True:
            header = handle.readline()
            if not header:
                break
            number += 1
            lines = [handle.readline() for _ in range(3)]
            if not all(lines):
                raise ValueError(f"{path}: truncated FASTQ record {number}")
            seq, plus, qual = [line.rstrip("\r\n") for line in lines]
            fields = header[1:].split()
            if not header.startswith("@") or not fields or not plus.startswith("+"):
                raise ValueError(f"{path}: malformed FASTQ record {number}")
            name = fields[0]
            if name.endswith(("/1", "/2")):
                name = name[:-2]
            if not name or not seq or len(seq) != len(qual):
                raise ValueError(f"{path}: invalid sequence/quality lengths at record {number}")
            if any(b not in "ACGTN" for b in seq):
                raise ValueError(f"{path}: expected uppercase A/C/G/T/N at record {number}")
            if any(not 33 <= ord(q) <= 126 for q in qual):
                raise ValueError(f"{path}: invalid Phred+33 quality at record {number}")
            yield name, seq, qual


def empty_stats(mode):
    return {"schema": SCHEMA, "quality_mode": mode,
            "mates": [{"reads": 0, "bases": 0, "sha256_sum": "0" * 64} for _ in range(2)]}


def accumulate(stats, mate, name, seq, qual):
    if not seq or "*" in (seq, qual) or len(seq) != len(qual):
        raise ValueError(f"{name}: missing or inconsistent sequence/quality")
    if stats["quality_mode"] == "pipeline":
        qual = "".join("!" if b in "Nn" or q == "#" else q for b, q in zip(seq, qual))
    # Length framing prevents ambiguity between fields. Include identity and mate,
    # and sum (rather than XOR) to retain duplicate multiplicity.
    digest = hashlib.sha256()
    for field in (name, str(mate + 1), seq, qual):
        value = field.encode("ascii")
        digest.update(len(value).to_bytes(8, "big"))
        digest.update(value)
    entry = stats["mates"][mate]
    entry["reads"] += 1
    entry["bases"] += len(seq)
    total = (int(entry["sha256_sum"], 16) + int.from_bytes(digest.digest(), "big")) % MODULUS
    entry["sha256_sum"] = f"{total:064x}"


def fastq_stats(f1, f2, mode="pipeline"):
    stats = empty_stats(mode)
    for r1, r2 in zip_longest(read_fastq(f1), read_fastq(f2)):
        if r1 is None or r2 is None:
            raise ValueError("R1 and R2 have different numbers of records")
        if r1[0] != r2[0]:
            raise ValueError(f"FASTQ mates have different names: {r1[0]!r}, {r2[0]!r}")
        accumulate(stats, 0, *r1)
        accumulate(stats, 1, *r2)
    if not stats["mates"][0]["reads"]:
        raise ValueError("Empty FASTQ input")
    return stats


def alignment_stats(path, reference, mode="pipeline", samtools="samtools", threads=2):
    stats = empty_stats(mode)
    args = SimpleNamespace(input=str(path), reference=reference, threads=threads,
                           samtools=samtools, name_collated=True)
    with SamtoolsSamReader(args) as stream:
        for line in stream:
            if line.startswith("@"):
                continue
            rec = parse_record(line.rstrip("\r\n").split("\t"))
            if not rec.is_primary():
                continue
            if not rec.flag & 1 or bool(rec.flag & 64) == bool(rec.flag & 128):
                raise ValueError(f"{rec.qname}: expected paired primary with exactly one mate flag")
            seq, qual = rec.restored_payload()
            accumulate(stats, 0 if rec.flag & 64 else 1, rec.qname, seq, qual)
    return stats


def load_stats(path):
    with open(path, encoding="utf-8") as handle:
        stats = json.load(handle)
    if not isinstance(stats, dict) or stats.get("schema") != SCHEMA or stats.get("quality_mode") not in ("pipeline", "strict"):
        raise ValueError("Unsupported statistics schema or quality mode")
    mates = stats.get("mates")
    if not isinstance(mates, list) or len(mates) != 2:
        raise ValueError("Statistics must contain two mates")
    for mate in mates:
        if not isinstance(mate, dict):
            raise ValueError("Invalid mate statistics")
        for key in ("reads", "bases"):
            if type(mate.get(key)) is not int or mate[key] <= 0:
                raise ValueError(f"Invalid statistics field: {key}")
        digest = mate.get("sha256_sum", "")
        if not isinstance(digest, str) or len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest):
            raise ValueError("Invalid statistics digest")
    if mates[0]["reads"] != mates[1]["reads"]:
        raise ValueError("Statistics have different mate counts")
    return stats


def comparison(expected, observed):
    differences = []
    for index in range(2):
        for key in ("reads", "bases", "sha256_sum"):
            if expected["mates"][index][key] != observed["mates"][index][key]:
                differences.append(f"R{index + 1}.{key}")
    return {"match": not differences, "differences": differences,
            "expected": expected, "observed": observed}


def write_json(path, data):
    path = Path(path)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(data, handle, indent=2)
            handle.write("\n")
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def protect_output(output, inputs):
    if output and any(Path(output).resolve() == Path(p).resolve() or
                      (Path(output).exists() and Path(p).exists() and os.path.samefile(output, p))
                      for p in inputs if p and p != "-"):
        raise ValueError("Output must not overwrite an input or reference")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    stats_parser = commands.add_parser("stats", help="Validate paired FASTQs and save expected payload digests")
    stats_parser.add_argument("--f1", required=True)
    stats_parser.add_argument("--f2", required=True)
    stats_parser.add_argument("--quality-mode", choices=("pipeline", "strict"), default="pipeline")
    stats_parser.add_argument("-o", "--output", required=True)
    check_parser = commands.add_parser("check", help="Compare restored primary reads with expected digests")
    check_parser.add_argument("-i", "--input", required=True)
    check_parser.add_argument("-T", "--reference")
    check_parser.add_argument("--stats", required=True)
    check_parser.add_argument("--report")
    check_parser.add_argument("--samtools", default="samtools")
    check_parser.add_argument("--threads", type=int, default=2)
    args = parser.parse_args(argv)
    try:
        if args.command == "stats":
            protect_output(args.output, [args.f1, args.f2])
            stats = fastq_stats(args.f1, args.f2, args.quality_mode)
            write_json(args.output, stats)
            print(f"Validated {stats['mates'][0]['reads']} read pairs", file=sys.stderr)
            return 0
        if args.threads < 1:
            raise ValueError("--threads must be positive")
        protect_output(args.report, [args.input, args.reference, args.stats])
        expected = load_stats(args.stats)
        observed = alignment_stats(args.input, args.reference, expected["quality_mode"], args.samtools, args.threads)
        report = comparison(expected, observed)
        if args.report:
            write_json(args.report, report)
        print("Match" if report["match"] else "Mismatch: " + ", ".join(report["differences"]))
        return 0 if report["match"] else 1
    except (OSError, ValueError, EOFError, BamRevertError) as error:
        print(f"fastq_check: error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
