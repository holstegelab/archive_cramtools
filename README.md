# CRAM roundtrip

Store aligned paired reads together with the bases and qualities removed during
trimming, then reconstruct paired FASTQs from the resulting CRAM. This repository
extracts the tag annotation, restoration and checksum code from the Holstege lab's
`short_read_analyzing_pipeline_Snakemake` project.

The storage format is ordinary, reference-based CRAM with four custom SAM tags:
`YB/YQ` for a removed prefix and `ZB/ZQ` for a removed suffix, expressed in the
alignment's orientation. The underlying C annotation code and optional clipping
code are included, as is the original native FASTQ checksum implementation.
No Snakemake, scheduler, storage service or lab configuration is required.

**Restoration recovers read payloads, not byte-identical original FASTQ files.**
Order, header comments, plus-line comments and gzip metadata are not retained.
The default validation mode accepts the lab pipeline's quality normalization;
use `--quality-mode strict` to require exact stored qualities. See
[validation](docs/validation.md) for the precise guarantees.

## Install and test

Linux/POSIX, Python 3.10+, samtools 1.17+, pigz, a C compiler, make and zlib headers
are required. Python's standard library is sufficient for the main commands.

```bash
conda env create -f environment.yml
conda activate cram-roundtrip
make
make test
make demo
```

Alternatively install those dependencies yourself and run `make`. If headers or
libraries are outside the compiler's default search path, pass `CPPFLAGS` and
`LDFLAGS` to make. The optional native checksum extension needs Python development
headers and setuptools: `make native-check`.

`make demo` generates three tiny synthetic read pairs and a synthetic reference
in `demo-output/`, stores a CRAM, restores FASTQs and compares their names, mate
numbers, sequences and exact qualities. It includes forward/reverse alignments,
unmapped reads, supplementary/secondary alignments, and a pair absent from the
alignment. It refuses to reuse an existing demo directory. To choose another:

```bash
python examples/roundtrip.py --output-dir /tmp/my-cram-demo
```

## Create a verified CRAM

Start with an alignment produced from 3′-trimmed reads, plus the **original,
untrimmed** paired FASTQs. Creation handles one read group per invocation.
Sequences must be uppercase A/C/G/T/N and qualities must be Phred+33. R1 and R2
must have matching names in the same order; each pair name must be unique within
that read group. Existing restoration tags and arbitrary internal/5′ trimming
are not supported inputs for annotation. See [the format contract](docs/format.md).

```bash
samtools faidx reference.fa

python scripts/create_cram.py \
  --input aligned.bam \
  --f1 original_R1.fastq.gz --f2 original_R2.fastq.gz \
  --reference reference.fa \
  --output-dir stored_reads \
  --threads 4
```

SAM, BAM or CRAM input is accepted, with an explicit reference. FASTQs can be
plain or gzip compressed. Input is grouped by QNAME automatically. For aligner
output already grouped by QNAME and close to FASTQ order, add `--name-collated`;
this can substantially reduce the annotation step's memory use.

The new output directory appears only after the final CRAM has been fully decoded,
its reconstructed payload checked against the source FASTQs, and its index built.
Existing output directories are refused. Input files are never edited.

| Output | Purpose |
| --- | --- |
| `reads.cram`, `reads.cram.crai` | Coordinate-sorted tagged alignment and index |
| `fastq.stats.json` | Expected read-payload digests; retain for later checks |
| `verification.json` | Expected/observed counts, base totals and digests |
| `merge.stats.tsv` | Original annotation statistics |
| `badmap_R1.fastq.gz`, `badmap_R2.fastq.gz` | Original annotator's diagnostic FASTQs for poorly aligned reads |
| `creation.json` | Tool version, samtools version and creation settings |

The `badmap` files are diagnostics, not complete restored FASTQs. Keep the exact
reference FASTA and its `.fai` with your archive or in a durable reference store;
the output directory does not copy them.

For large datasets use `--tmp-dir /scratch/path` and tune `--sort-memory 768M`
(per samtools sort thread). The Python checks are streaming and use bounded
memory. The inherited annotation lookup can buffer many original FASTQs when
alignment and FASTQ order differ; `--sort-memory` does not limit that buffer.

## Restore paired FASTQs

```bash
python scripts/bam_revert.py \
  --input stored_reads/reads.cram --reference reference.fa \
  --f1 restored_R1.fastq.gz --f2 restored_R2.fastq.gz \
  --stats restored.stats.tsv --threads 4
```

Reversion automatically collates records, uses primary alignments, restores the
tags and reverses strand orientation. It requires exactly one primary R1 and R2
per `(QNAME, RG)` group. `--include-read-group-in-name` writes `QNAME:RG/1` and
`QNAME:RG/2`, useful when names recur across read groups. That option changes the
names compared by the standalone checksum tool. Without it the names are
`QNAME/1` and `QNAME/2`. Output files are gzip compressed regardless of suffix.

## Verify an existing CRAM

Save expected statistics once, while original FASTQs are available:

```bash
python scripts/fastq_check.py stats \
  --f1 original_R1.fastq.gz --f2 original_R2.fastq.gz \
  --quality-mode strict --output original.stats.json

python scripts/fastq_check.py check \
  --input stored_reads/reads.cram --reference reference.fa \
  --stats original.stats.json --report comparison.json
```

Or use `stored_reads/fastq.stats.json`, which creation already saved. `check`
returns **0 for a match, 1 for a mismatch, 2 for an input/tool error**. A failed
read never produces a successful result. A JSON report is written for completed
comparisons, including mismatches; on an input error an older report is left
untouched, so always check the exit status. Reversion and creation return nonzero
on failure. Run `--help` on each script for all options.

## Further documentation

- [Tag format, orientation and supported inputs](docs/format.md)
- [Checksum algorithms, quality policies and limitations](docs/validation.md)
- [Low-level commands and optional clipping](docs/advanced.md)
- [Source provenance and standalone changes](docs/provenance.md)
- [Troubleshooting](docs/troubleshooting.md)

This repository contains source code and synthetic examples only. There is no
bundled human sequencing data, encryption key or reference genome. Licensing
status is recorded in [NOTICE](NOTICE).
