# Validation

## Default standalone check

`fastq_check.py` validates four-line paired FASTQ structure, mate names, nonempty
sequences, sequence/quality lengths, supported bases, printable qualities, equal
mate counts and gzip decoding through EOF. The source quality encoding must be
Phred+33; printable Phred+64 cannot always be distinguished automatically. Convert
legacy encodings before using these commands.

The same checks apply to streamed input. FASTQ checksum calculation opens each
input once and detects gzip without seeking or reopening. R1 and R2 can use
separate pipes, with at most one mate on stdin. Alignment checking accepts
SAM/BAM/CRAM on stdin through `--input -`. See [streaming examples](streaming.md).
For a single stream containing both mates, `stats --interleaved --input -`
checks consecutive R1/R2 records with the same rules and produces the same JSON
as `--f1`/`--f2` for the corresponding separate files.

Each read contributes a SHA-256 digest of its normalized name, mate number,
sequence and quality. Each field is prefixed with its length as an unsigned
eight-byte big-endian integer. Digests are interpreted as big-endian integers and
summed modulo 2^256 separately for each mate. The check compares these sums,
read counts and base totals after restoring primary SAM/BAM/CRAM records. It is
independent of read ordering and includes duplicate multiplicity. Counts and
hashes are probabilistic integrity evidence, not an exact comparison or an
authenticated cryptographic signature.

The JSON schema is `cram-roundtrip-sha256-sum-v1`. It records the quality mode so a
later comparison uses the same policy. The check does not include RG IDs; use it
per read group with the original names. `--include-read-group-in-name` during
reversion deliberately changes those names.

## Optional creation check

`create_cram.py` runs the FASTQ checksum pass and compares the finished CRAM by
default. With `--fast`, it skips both passes and omits `.fastq.stats.json` and
`.verification.json`. The `.creation.json` record sets `verification` to `skipped`
instead of `passed`. Sorting, indexing, `samtools quickcheck` and the annotator's
own checks still run, but they do not establish that the original read information
is preserved. `--quality-mode` only controls verification and has no effect with
`--fast`.

To check a CRAM created with `--fast`, run `fastq_check.py stats` on the original
FASTQs, then `fastq_check.py check` on the CRAM, as shown in the README.

## Quality modes

| Mode | Comparison rule |
| --- | --- |
| `pipeline` (default) | Treat `#` as `!`; ignore differences in quality at `N/n` bases by mapping them to `!` |
| `strict` | Compare the actual ASCII quality values without normalization |

The pipeline mode preserves the normalization used by the lab's native checker
for aligners that change low-quality values or qualities at N bases. Normalization
only affects the comparison, not what the reverter writes. A passing pipeline
comparison does **not** establish exact original quality recovery. Strict mode
rejects these differences, even if the pipeline mode accepts them. Neither mode
compares original header comments, plus lines, ordering or gzip metadata.

## Original native checker

`scripts/fastcheck.c`, `setup.py`, `fastcheck_loader.py` and `hash_utils.py` retain
the original pipeline's FNV-1a implementation. Build it with `make native-check`.
Its Python interface provides `fastq_stats(R1, R2)`,
`fastq_stats_interleaved(path_or_stdin)` and `sam_stats(path_or_stdin)`. This code is
included for compatibility, inspection and speed comparisons; the supported
creation wrapper uses the standalone check above.

The native algorithm XORs independent 64-bit hashes of sequences and qualities
per mate, and returns read/base counts. It omits names and does not bind sequence
and quality to one record. Equal duplicate hashes cancel in pairs; equal XORs
and counts cannot prove read-for-read identity. Its historic FASTQ parser also
uses a heuristic for Phred+64 detection, and some truncated-input/subprocess
failures are not rejected. Use it only with known-valid source files, never as
the sole archival acceptance criterion. Its shell-based filenames must not
contain apostrophes. The standalone entry points handle such filenames safely.

The old native tuples/TSV fields and new JSON statistics are intentionally
different formats; do not substitute one for the other. The source tests for the
native quality-buffer boundary are retained.

## Verification coverage

The local verification on 2026-10-26 passed the test suite and the complete synthetic
demo with Python 3.12.3 and samtools/HTSlib 1.17.

`make test` runs integration tests using synthetic data. Tests cover actual CRAM
encoding/decoding, both strand orientations, restoring clipped and unmapped
payloads, re-adding fully filtered pairs, ignoring secondary/supplementary
records, duplicate names in separate RGs during reversion, compressed/plain
FASTQs, filenames with spaces/apostrophes, checks of read identity and quality
mode, malformed input, missing tags, failed samtools execution and preservation
of existing output. Sample-prefix tests cover multiple samples in a shared
directory, existing sidecars, broken symlinks and rollback on publication conflicts.
Fast-mode tests cover skipped checksum passes, recovery of dropped pairs,
verification afterwards, deferred mismatches and rejection of corrupt gzip input.
Streaming tests cover SAM/BAM/CRAM stdin, named alignment/FASTQ pipes, split gzip
headers, corrupt or truncated input, distinct mate streams and read-group checks.
Interleaved-input tests cover both quality modes, file/stdin/FIFO input, a `tee`
branch, missing mates, mismatched names and preservation of existing statistics.
`make test` requires the external tools and fails if they
are missing.

These tests demonstrate the included examples and contracts. They are not a
new whole-genome-scale performance benchmark or validation of every historical
aligner convention.
