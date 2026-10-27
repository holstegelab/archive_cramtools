# CRAM roundtrip

Keep alignments and recoverable read data together in one CRAM file. This tool
stores bases and quality scores removed by trimming in extra SAM tags, so paired
FASTQs can be reconstructed later.

It also re-adds read pairs that are present in the original FASTQs but were
dropped from the alignment, for example by an earlier filtering step.

The files use standard, reference-based CRAM with four extra tags:
`YB/YQ` and `ZB/ZQ`. The tools run independently of the original analysis pipeline.
See [how the tags work](docs/format.md) for details.

Recovery restores the read sequences and quality scores stored in the CRAM and
its tags. Original read order, FASTQ comments and gzip metadata are not preserved.
Quality scores changed by an aligner cannot always be recovered exactly. The
verification commands below let you check this against the original FASTQs.

## Motivation

FASTQs are often kept so reads can be analysed again with a different aligner or
reference. Keeping both the FASTQs and the alignments takes up substantial
storage. A CRAM that also preserves the trimmed read data can serve both purposes.

We compared original gzipped FASTQs with matching files from our **CRAM 3.1
`archive`** workflow. The CRAM totals include encryption and indexes.

| Data | Samples | Original FASTQ.gz | CRAM archives + indexes | Less storage |
| --- | ---: | ---: | ---: | ---: |
| Whole-exome (WES) | 1,000 | 3.312 TB | 0.819 TB | **75.3%** |
| Whole-genome (WGS) | 448 | 23.222 TB | 4.300 TB | **81.5%** |

This saves **2.493 TB** for the exomes and **18.922 TB** for the genomes, while
also keeping the alignments. Compression depends on the data. See the
[storage comparison](docs/storage-comparison.md) for the measurements and method.

## Outline

1. **Prepare the inputs:** original paired FASTQs, their alignments and the
   matching reference FASTA.
2. **[Create a CRAM](#create-a-cram):** add the trimmed read data and
   write files named after your sample.
3. **[Check against the original FASTQs](#verify-that-a-cram-preserves-the-original-fastq-reads):**
   compare the recovered read information with the originals. Creation runs this
   check by default; saved checksums let you repeat it later.
4. **[Restore FASTQs](#restore-paired-fastqs):** recover paired reads for a new
   analysis whenever you need them.

You can also [stream inputs](#stream-inputs) into creation and checksum commands.

## Install and test

The easiest setup is the supplied conda environment:

```bash
conda env create -f environment.yml
conda activate cram-roundtrip
make
make test
make demo
```

You can also install the dependencies yourself: Linux/POSIX, Python 3.10+,
samtools 1.17+, pigz, a C compiler, make and zlib development headers. Then run
`make`. The main Python commands need no additional Python packages.

`make demo` creates three small, synthetic read pairs and a reference in
`demo-output/`. It creates a CRAM, restores the FASTQs and checks that the read
names, sequences and qualities match. To run the demo in a different, new folder:

```bash
python examples/roundtrip.py --work-dir /tmp/my-cram-demo
```

The optional original C checksum extension can be built with `make native-check`.
It requires setuptools and Python development headers; see
[advanced use](docs/advanced.md) for its interface and limitations.

## Create a CRAM

Supply the alignment, the **original untrimmed FASTQs** and the matching reference.
If reads were trimmed before alignment, this command supports trimming from the
**3′ end** only.

The inputs must meet these requirements:

- Process one read group at a time, with a unique name for each read pair.
- R1 and R2 must contain matching read names in the same order.
- Sequences must use uppercase A/C/G/T/N and qualities must use Phred+33 encoding.
- The alignment must not already contain the restoration tags.

See [supported inputs](docs/format.md) for details, or
[advanced use](docs/advanced.md#encode-existing-tagged-alignments) to encode an
alignment that already has the tags.

```bash
samtools faidx reference.fa

python scripts/create_cram.py \
  --input SAMPLE001.aligned.bam \
  --f1 SAMPLE001_R1.fastq.gz --f2 SAMPLE001_R2.fastq.gz \
  --reference reference.fa \
  --output-prefix SAMPLE001 \
  --threads 4
```

The alignment can be SAM, BAM or CRAM. FASTQs can be plain or gzip compressed.
The creation command uses `samtools collate` to put records with the same read
name together before adding tags. It then uses `samtools sort` to write a
position-sorted CRAM that can be indexed. `--name-collated` skips the first step
when the input is already grouped by read name.

Encoding currently uses samtools' default CRAM settings; it does not select
the CRAM 3.1 `archive` profile used in the storage comparison above.

Read pairs missing entirely from the alignment are re-added from the original
FASTQs as unmapped reads marked as QC-failed. The restore command includes these
pairs, so reads dropped by an earlier filter can be recovered too. This does not
repair a pair with only one mate missing from the alignment.

`--output-prefix` (or `-o`) sets the name used for every output file:

- `--output-prefix SAMPLE001` writes files in the current directory.
- `--output-prefix results/SAMPLE001` writes them in an existing `results/` folder.

Multiple samples can share the same folder by using different prefixes. Existing
output files are never overwritten.

By default, the tool reads back the complete CRAM and checks its reconstructed
reads against the original FASTQs before making the outputs available. Add
`--quality-mode strict` to require exact quality matches; the default `pipeline`
mode allows the differences described in the verification section below.

Add **`--fast`** to the creation command to skip the separate FASTQ checksum pass
and the final comparison. The original FASTQs are still read to add the tags and
recover dropped pairs. Sorting, indexing and `samtools quickcheck` still run,
but they do not establish that the CRAM preserves the original reads.

With `--fast`, the tool omits `.fastq.stats.json` and `.verification.json` and
records `"verification": "skipped"` in `.creation.json`. Normal creation records
`"verification": "passed"`. `--quality-mode` has no effect in fast mode. To verify
later, use the [comparison commands below](#verify-that-a-cram-preserves-the-original-fastq-reads)
with the original FASTQs.

| Output | Contents |
| --- | --- |
| `SAMPLE001.cram`, `SAMPLE001.cram.crai` | Tagged alignment, sorted by position, and its index |
| `SAMPLE001.fastq.stats.json` | Checksums and counts from the original FASTQs; omitted with `--fast` |
| `SAMPLE001.verification.json` | Results of the comparison between the CRAM and original FASTQs; omitted with `--fast` |
| `SAMPLE001.merge.stats.tsv` | Statistics on the added read data |
| `SAMPLE001.badmap_R1.fastq.gz`, `SAMPLE001.badmap_R2.fastq.gz` | A diagnostic subset of poorly aligned reads |
| `SAMPLE001.creation.json` | Sample name, software versions, settings and verification status; written last to mark completion |

Keep the exact reference FASTA and its `.fai` index with your archive or in a
reliable reference store. They are needed to read the CRAM. The `badmap` files are
only a diagnostic subset; use the restore command below to recover all read pairs.

For large datasets, `--tmp-dir /scratch/path` chooses where sorting uses temporary
space. `--sort-memory 768M` sets memory per sorting thread. If the alignment is
already grouped by read name, `--name-collated` skips regrouping. Keeping that
order close to the FASTQ order can reduce memory use when adding tags; this step
has a separate memory buffer that `--sort-memory` does not limit.

## Stream inputs

Creation and checksum checking accept piped inputs. The CRAM, checksums and
reports are still written to named files.

To create a CRAM from an alignment stream, use **`--input -`**:

```bash
set -euo pipefail
samtools view -u SAMPLE001.aligned.bam \
  | python scripts/create_cram.py \
      --input - \
      --f1 SAMPLE001_R1.fastq.gz --f2 SAMPLE001_R2.fastq.gz \
      --reference reference.fa --output-prefix SAMPLE001 --threads 4
```

SAM, BAM and CRAM streams are supported. With `--input -`, creation reads the
incoming alignment once and runs the creation steps described above.
Verification stays enabled unless you add `--fast`.

The checksum commands read records directly and do **not** sort or group them.
The original pipeline supplies **both FASTQ mates in one input stream**, with
each R1 record immediately followed by its R2 record. The standalone checksum
command accepts the same layout with `--interleaved --input -`:

```bash
set -euo pipefail
pigz -dc SAMPLE001.interleaved.fastq.gz \
  | python scripts/fastq_check.py stats \
      --interleaved --input - \
      --quality-mode strict --output SAMPLE001.original.fastq.stats.json
```

Here the input order is `pairA/1`, `pairA/2`, `pairB/1`, `pairB/2`, and so on.
Any producer with this layout can feed the checksum command directly. As in the
original pipeline, `tee` can also send the same stream to the checksum reader
and another processing step; see the [tee example](docs/streaming.md#check-the-same-stream-while-another-command-processes-it).

To check a streamed CRAM against saved FASTQ checksums:

```bash
set -euo pipefail
cat SAMPLE001.cram \
  | python scripts/fastq_check.py check \
      --input - --reference reference.fa \
      --stats SAMPLE001.original.fastq.stats.json --report SAMPLE001.comparison.json
```

For separate R1/R2 sources, **two separate named pipes** remain supported through
`--f1` and `--f2`. Both input modes read once and detect plain or gzip-compressed data.
Creation accepts these FASTQ pipes with `--fast`; normal verified creation needs
FASTQ files it can read twice. The reference FASTA and its index remain files.
See [streaming examples](docs/streaming.md) for both FASTQ layouts, producer exit
checks and the supported inputs for each command.

## Verify that a CRAM preserves the original FASTQ reads

Use this check to confirm that the CRAM preserves the read information from the
**original paired FASTQs**. It reconstructs reads using the extra tags and compares their
names, mate numbers, sequences and qualities with the originals. It also checks
read counts and total numbers of bases. Checksums allow the comparison to work
even when the read order has changed.

First calculate checksums from the original FASTQs, then check the CRAM against
them:

```bash
python scripts/fastq_check.py stats \
  --f1 SAMPLE001_R1.fastq.gz --f2 SAMPLE001_R2.fastq.gz \
  --quality-mode strict --output SAMPLE001.original.fastq.stats.json

python scripts/fastq_check.py check \
  --input SAMPLE001.cram --reference reference.fa \
  --stats SAMPLE001.original.fastq.stats.json --report SAMPLE001.comparison.json
```

Choose the quality policy when calculating the FASTQ checksums:

| Mode | What it checks |
| --- | --- |
| `strict` | Quality values must match exactly. |
| `pipeline` (default) | Treats `#` and `!` qualities as equivalent and ignores quality differences at `N` bases. |

A match in `pipeline` mode can therefore still include quality differences. Both
modes check read sequences and names, but neither compares FASTQ comments,
original ordering or gzip bytes. See [validation](docs/validation.md) for the full
comparison rules and checksum limitations.

Creation without `--fast` already saves `SAMPLE001.fastq.stats.json` from the
original FASTQs. You can use it with `--stats` to check the CRAM again without rereading the
originals. The saved file also records which quality policy to apply.

The `check` command returns **0 for a match, 1 for a mismatch and 2 for an input
or tool error**. The JSON report describes a completed comparison, including any
mismatch. If the command fails before completing the comparison, an older report
may remain, so always check the exit status.

## Restore paired FASTQs

Use the supplied restore command to put the tagged bases and quality scores back
into the reads:

```bash
python scripts/bam_revert.py \
  --input SAMPLE001.cram --reference reference.fa \
  --f1 SAMPLE001.restored_R1.fastq.gz --f2 SAMPLE001.restored_R2.fastq.gz \
  --stats SAMPLE001.revert.stats.tsv --threads 4
```

The script groups mates by read name and read group, restores the trimmed pieces
and returns reads to their original orientation. It requires one primary R1 and
one primary R2 for each pair. Both output FASTQs are gzip compressed.

If read names repeat across read groups, add `--include-read-group-in-name`.
This writes names such as `readname:readgroup/1` and `readname:readgroup/2`, so the
groups remain distinguishable. Those added read-group names also affect later
comparisons with the original FASTQs.

## Further documentation

- [Tag format and supported inputs](docs/format.md)
- [Checksums, quality policies and limitations](docs/validation.md)
- [Streaming inputs and named-pipe examples](docs/streaming.md)
- [Advanced use and optional clipping](docs/advanced.md)
- [Source provenance and standalone changes](docs/provenance.md)
- [Troubleshooting](docs/troubleshooting.md)

Run the commands above with `--help` for their options. This repository includes synthetic
examples; it contains no human sequencing data or reference genome. Licensing
status is recorded in [NOTICE](NOTICE).
