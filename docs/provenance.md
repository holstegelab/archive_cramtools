# Source provenance

Extracted on 2026-10-26 from the working tree of
`short_read_analyzing_pipeline_Snakemake`, at Git HEAD
`a80d95c09231752069f209978e6a2f789db06a2e`.

The source checkout was not clean. In particular, `scripts/bam_merge.c` included
working-tree changes for interleaved FASTQ input. The extraction intentionally
used the current files and records their SHA-256 values and modification status
in [source-manifest.json](source-manifest.json). Source hashes identify the files
before the standalone changes described below. The parent repository was not
edited by this extraction.

## Retained implementations

- `bam_merge.c`: original tag annotation, supplementary CIGAR updates, badmap
  diagnostics and re-addition of original pairs absent from the alignment.
- `bam_dechimer.c`: original optional reversible clipping implementation.
- `bam_revert.py`: original standalone tag restoration, samtools collation,
  read-group-aware pairing and temporary output publication.
- `fastcheck.c`, `setup.py`, `fastcheck_loader.py`, `hash_utils.py`: original
  native checksum code and build/loading helpers.
- `test_bam_revert.py`, `test_fastcheck.py`: source regression tests.

## Standalone adaptations

Trailing whitespace in the extracted C annotation/clipping sources was normalized.

- `bam_merge.c`: shell quoting for source FASTQ paths, plain FASTQ support via
  pigz, propagation of decompression/write failures, rejection of malformed SAM
  records and existing restoration tags.
- `bam_revert.py`: reject incomplete/inconsistent tag pairs, absent payloads and
  mapped hard clips without matching tags; resolve output aliases and report
  malformed numeric fields as CLI errors.
- New `create_cram.py`: checks original FASTQs, groups input names if needed,
  invokes the C annotator and samtools, fully decodes/checks the finished CRAM,
  then publishes the CRAM, index and reports with a sample filename prefix.
- New `fastq_check.py`: streaming input validation and name/mate-aware SHA-256
  sums with explicit pipeline/strict quality modes. This is an addition to, not
  a byte-compatible replacement for, the historic FNV-1a XOR format.
- New standalone build/environment, integration tests, synthetic demonstration
  and English documentation.

The storage tag semantics remain compatible with the source tools. The new
quality checks can reject malformed legacy files that the original reverter
accepted. Source scripts not needed for this scope, cluster configuration,
encryption/delivery code, credentials and real datasets were not copied.

## Updating the extraction

Review upstream changes to the listed files, update the source manifest, reapply
the focused standalone adaptations and run `make test` and the demo. Do not copy
compiled `.so` files or executables from a different Python/platform environment.
