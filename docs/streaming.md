# Streaming inputs

Use pipes to feed data directly from another command without first writing a
complete input file. FASTQ checksum calculation accepts either one stream with
alternating R1/R2 records, as used by the original pipeline, or two separate
mate streams. Outputs keep their usual filenames.

| Command | Streaming input |
| --- | --- |
| `create_cram.py --input -` | SAM, BAM or CRAM on stdin; verification remains enabled by default |
| `create_cram.py --input alignment.pipe` | SAM, BAM or CRAM from a named pipe |
| `create_cram.py --fast --f1 R1.pipe --f2 R2.pipe` | Separate original FASTQ streams through named pipes, plain or gzip compressed |
| `fastq_check.py stats --interleaved --input -` | Both mates on stdin, alternating R1/R2 records, plain or gzip compressed |
| `fastq_check.py stats --f1 R1.pipe --f2 R2.pipe` | Separate original FASTQ streams, plain or gzip compressed |
| `fastq_check.py stats --f1 - --f2 R2.fastq.gz` | R1 on stdin and R2 from a file or pipe; the roles can be reversed |
| `fastq_check.py check --input -` | SAM, BAM or CRAM on stdin, compared with a saved statistics JSON file |
| `bam_revert.py --input -` | SAM, BAM or CRAM on stdin; the usual two FASTQ.gz outputs remain files |

An input pipe is consumed once and cannot be reused for another command. With
`--f1`/`--f2`, each mate needs its own stream: never point them at the same pipe
or pass `--f1 - --f2 -`. Use `--interleaved` for a single paired FASTQ stream.
The reference FASTA and `.fai` index remain regular files.

## Pipe an alignment into creation

```bash
set -euo pipefail
samtools view -u SAMPLE001.aligned.bam \
  | python scripts/create_cram.py \
      --input - \
      --f1 SAMPLE001_R1.fastq.gz --f2 SAMPLE001_R2.fastq.gz \
      --reference reference.fa --output-prefix SAMPLE001 --threads 4
```

You can replace `samtools view` with an aligner or another producer of valid
SAM/BAM/CRAM. Creation reads the alignment once and checks the read-group headers
as they pass through. It still re-adds dropped pairs from the original FASTQs.

The `create_cram.py` wrapper runs `samtools collate` before annotation and
`samtools sort` afterwards to produce an indexed, position-sorted CRAM.
`--name-collated` skips collation when records are already grouped by read name.
These creation steps can use temporary files, controlled by `--tmp-dir` and
`--sort-memory`.

The `bam_merge` annotator itself does not sort or group records; it expects
records for each read name to be together already. Neither checksum command
sorts or groups its input. These are differences between the commands, not
requirements imposed by streaming.

Normal verification calculates FASTQ checksums before consuming the alignment,
then compares them with the finished CRAM. The upstream alignment producer may
wait while those checksums are calculated. `--fast` skips these checks as it does
for file input. Neither mode writes an extra copy of the input alignment.

## Calculate checksums from one paired FASTQ stream

The original pipeline combines the two mates into alternating R1/R2 records and
uses `tee` to supply that stream to both the checksum reader and adapter
processing. Its checksum command uses `--interleaved --input -`. The standalone
command now accepts the same input layout:

```bash
set -euo pipefail
pigz -dc SAMPLE001.interleaved.fastq.gz \
  | python scripts/fastq_check.py stats \
      --interleaved --input - --quality-mode strict \
      --output SAMPLE001.original.fastq.stats.json
```

Each R1 record must be immediately followed by its matching R2 record. Both
mates are checked, including their names, sequences, qualities, read counts and
base totals. Incomplete pairs and malformed records fail the command.
`--input` can also name a plain/gzip file or named pipe; omitting it selects
stdin. Do not combine `--interleaved` with `--f1` or `--f2`.

This changes only how reads enter the checker. The same standalone SHA-256 JSON
statistics and quality modes are used for both layouts; this does not switch to
the original pipeline's legacy checksum format.

### Check the same stream while another command processes it

This follows the pipeline's `tee`/named-pipe arrangement. Here the downstream
command simply recompresses the reads; replace that final `pigz` command with
your processing step. No extra uncompressed FASTQ file is created.

```bash
(
  set -euo pipefail
  stream_dir=$(mktemp -d)
  checker_pid=
  cleanup_streams() {
    if [[ -n "$checker_pid" ]]; then
      kill "$checker_pid" 2>/dev/null || true
      wait "$checker_pid" 2>/dev/null || true
    fi
    rm -rf "$stream_dir"
  }
  trap cleanup_streams EXIT
  mkfifo "$stream_dir/check.pipe"

  python scripts/fastq_check.py stats \
    --interleaved --input - --quality-mode strict \
    --output SAMPLE001.original.fastq.stats.json < "$stream_dir/check.pipe" &
  checker_pid=$!

  pigz -dc SAMPLE001.interleaved.fastq.gz \
    | tee "$stream_dir/check.pipe" \
    | pigz > SAMPLE001.forwarded.fastq.gz

  wait "$checker_pid"
  checker_pid=
)
```

`pipefail` checks the main pipeline, and `wait` checks the checksum reader. Both
must succeed before accepting the outputs.

## Calculate checksums from two separate FASTQ streams

This Bash example uses one named pipe per mate. The two producers run
concurrently, and their exit codes are checked along with the checksum command.
The subshell keeps its cleanup trap separate from your surrounding shell.

```bash
(
  set -euo pipefail
  stream_dir=$(mktemp -d)
  producer_pids=()
  cleanup_streams() {
    for producer_pid in "${producer_pids[@]}"; do
      kill "$producer_pid" 2>/dev/null || true
    done
    for producer_pid in "${producer_pids[@]}"; do
      wait "$producer_pid" 2>/dev/null || true
    done
    rm -rf "$stream_dir"
  }
  trap cleanup_streams EXIT
  mkfifo "$stream_dir/R1.pipe" "$stream_dir/R2.pipe"

  pigz -dc SAMPLE001_R1.fastq.gz > "$stream_dir/R1.pipe" &
  producer_pids+=("$!")
  pigz -dc SAMPLE001_R2.fastq.gz > "$stream_dir/R2.pipe" &
  producer_pids+=("$!")

  python scripts/fastq_check.py stats \
    --f1 "$stream_dir/R1.pipe" --f2 "$stream_dir/R2.pipe" \
    --quality-mode strict --output SAMPLE001.original.fastq.stats.json

  for producer_pid in "${producer_pids[@]}"; do
    wait "$producer_pid"
  done
  producer_pids=()
)
```

The producers can supply gzip-compressed bytes directly instead of decompressing
first; the checksum reader detects compression without reopening the pipe. It
checks the same FASTQ structure, mate names, qualities and gzip integrity as for
file inputs. Start an independent producer for each mate so both streams can
make progress.

Only accept the result if the command **and both producers** succeed. An
external producer can fail after emitting complete read pairs, which can look
like a valid end of input to the reader. `set -o pipefail` covers ordinary `|`
pipelines; background producers need the explicit `wait` calls above. On a
failure, discard any result from that run, even if a JSON report exists.

## Create using two FASTQ streams

Use the same two-pipe setup above, replacing the `fastq_check.py stats` command
with the following command. Creation uses separate FASTQ inputs; the
`--interleaved` option above belongs to checksum calculation.

```bash
python scripts/create_cram.py \
  --input SAMPLE001.aligned.bam \
  --f1 "$stream_dir/R1.pipe" --f2 "$stream_dir/R2.pipe" \
  --reference reference.fa --output-prefix SAMPLE001 --threads 4 --fast
```

The alignment can also arrive on stdin with `--input -`. For creation, stdin is
reserved for the alignment; original FASTQ streams use named pipes.

FASTQ pipes require `--fast` because normal creation reads each original FASTQ
twice: once for checksums and once for annotation. Without `--fast`, the command
rejects FASTQ pipes before opening them. Use regular FASTQ files for verified
creation, or run the independent checksum and comparison commands with fresh
input streams later. The fast-mode creation record marks verification as skipped.

## Check or restore a streamed alignment

To compare streamed read data with previously saved FASTQ checksums:

```bash
set -euo pipefail
cat SAMPLE001.cram \
  | python scripts/fastq_check.py check \
      --input - --reference reference.fa \
      --stats SAMPLE001.original.fastq.stats.json --report SAMPLE001.comparison.json
```

The check reads primary records directly in any order, without sorting, grouping
or an alignment index. The saved statistics determine the quality comparison
policy. Exit codes remain 0 for a match, 1 for a mismatch and 2 for an input/tool
error.

The restoration command already accepts alignment input on stdin:

```bash
set -euo pipefail
cat SAMPLE001.cram \
  | python scripts/bam_revert.py \
      --input - --reference reference.fa \
      --f1 SAMPLE001.restored_R1.fastq.gz --f2 SAMPLE001.restored_R2.fastq.gz \
      --stats SAMPLE001.revert.stats.tsv --threads 4
```

Restoration groups mates internally and writes its normal two FASTQ.gz files.
No changes to the output format are needed to use streamed input.
