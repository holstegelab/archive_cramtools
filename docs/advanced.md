# Low-level use

The high-level creation command adds tags and produces a verified CRAM. It does
not realign, mark duplicates or modify the alignment with the optional
dechimerization algorithm.

## Annotation as a SAM stream

For custom workflows, the following is the extracted C annotator's interface.
The input must group all records for each QNAME; source FASTQs should be close to
that order to limit buffering. Work one read group at a time.

```bash
set -euo pipefail
samtools view -h aligned.name-grouped.bam \
  | scripts/bam_merge \
      -a original_R1.fastq.gz -b original_R2.fastq.gz \
      -ua badmap_R1.fastq.gz -ub badmap_R2.fastq.gz \
      -s merge.stats.tsv \
  | samtools sort --reference reference.fa -O CRAM -o tagged.cram -
```

`bam_merge` accepts plain paired FASTQs too. Its `--interleaved-fastq` option
accepts a plain file or named pipe, replacing `-a/-b`; SAM still arrives on stdin.
These low-level commands do not stage output files or perform the final check
automatically. Check the completed CRAM with `fastq_check.py` before accepting it.

## Optional reversible alignment clipping

`bam_dechimer.c` is the lab's additional overlap/chimeric-alignment clipping code.
It can move sequence and quality into the same four tags, so the original primary
payload remains recoverable. It changes alignment placement/CIGAR/flags and is
an optional processing decision, not necessary for CRAM encoding.

For example, start from an already annotated, name-grouped BAM, with mate fields
prepared by `samtools fixmate`:

```bash
set -euo pipefail
samtools fixmate -u -m tagged.name-grouped.bam - \
  | samtools view -h - \
  | scripts/bam_dechimer -i - -o - -s clipping.stats.tsv \
  | samtools sort --reference reference.fa -O CRAM -o clipped.cram -

python scripts/fastq_check.py check \
  --input clipped.cram --reference reference.fa --stats original.stats.json
```

The standalone wrapper does not reproduce the parent workflow's threshold-based
decision about when to run this step. The C tool supports `--loose_ends`,
`--min_align_length` (default 30) and `--max_read_dist` (default 100000).
Its historical `-d` switch has no effect: the C port does not implement the
internal validation it would disable. Always run the independent payload check
after clipping. Do not pass the result back through `bam_merge`.

## Encode existing tagged alignments

If the input already has the lab tags, skip annotation:

```bash
samtools sort --reference reference.fa -O CRAM -o stored.cram tagged.bam
samtools index stored.cram
python scripts/fastq_check.py check \
  --input stored.cram --reference reference.fa --stats original.stats.json
```

Use ordinary lossless CRAM settings; do not strip tags or discard quality scores.
Preserve unmapped primary reads, including QC-failed reads, in the archive.

## Native pipeline checksum API

After `make native-check`, this prints the original seven-value FASTQ tuple:

```bash
PYTHONPATH=scripts python -c \
  'import fastcheck; print(fastcheck.fastq_stats("R1.fastq.gz", "R2.fastq.gz"))'
```

The order is `(R1_sequence_xor, R1_quality_xor, R2_sequence_xor, R2_quality_xor,
pair_count, R1_base_count, R2_base_count)`. The SAM tuple is
`(R1_sequence_xor, R1_quality_xor, R1_read_count, R1_base_count,
R2_sequence_xor, R2_quality_xor, R2_read_count, R2_base_count)`.
See [the validation limitations](validation.md) before interpreting these values.
