# Storage comparison

Measured for **1,000 whole-exome samples (WES)** and **448 whole-genome samples
(WGS)**. For each sample, we matched all original gzipped FASTQs to one archived
encrypted CRAM and its index. Every sample had a complete set of files; none
were excluded because of a missing file.

| Data | Samples | Original FASTQ.gz, TB | CRAM + index, TB | Reduction | Space saved, TB |
| --- | ---: | ---: | ---: | ---: | ---: |
| WES | 1,000 | 3.311765 | 0.818769 | 75.277% | 2.492996 |
| WGS | 448 | 23.221744 | 4.299658 | 81.484% | 18.922087 |

All sizes use decimal units: 1 TB = 10^12 bytes. CRAM totals include encryption
overhead and indexes. The shared reference FASTA, backups/replicas and other
analysis products are excluded. Counts refer to processed sample/library records,
not a count of unique participants. Compression depends on the data and sequencing
technology.

## Whole-exome data

Each of the 1,000 samples has two original FASTQ.gz files.

| Stored files | Files | Bytes | Decimal TB |
| --- | ---: | ---: | ---: |
| Original paired FASTQ.gz | 2,000 | 3,311,764,715,580 | 3.311765 |
| Tagged CRAM, Crypt4GH encrypted (`.cram.c4gh`) | 1,000 | 818,751,459,024 | 0.818751 |
| CRAM indexes (`.cram.crai`) | 1,000 | 17,572,801 | 0.000018 |
| **CRAM archives plus indexes** | **2,000** | **818,769,031,825** | **0.818769** |

The aggregate reduction is **75.277%**, or **4.045-fold**, saving
**2,492,995,683,755 bytes (2.493 TB)** relative to keeping the gzipped FASTQs alone.
Per-sample reductions range from 71.71% to 78.88%, with a median of 75.32%.

## Whole-genome data

The 448 samples have **7,024 original FASTQ.gz files** in total. Some samples were
sequenced across multiple lanes or runs, so the comparison includes every FASTQ
listed for each sample, rather than assuming two files per sample. Each sample
has between 2 and 48 source FASTQ files, counting both mates.

| Stored files | Files | Bytes | Decimal TB |
| --- | ---: | ---: | ---: |
| Original paired FASTQ.gz, all lanes | 7,024 | 23,221,744,190,236 | 23.221744 |
| Tagged CRAM, Crypt4GH encrypted (`.cram.c4gh`) | 448 | 4,299,593,318,855 | 4.299593 |
| CRAM indexes (`.cram.crai`) | 448 | 64,311,641 | 0.000064 |
| **CRAM archives plus indexes** | **896** | **4,299,657,630,496** | **4.299658** |

The aggregate reduction is **81.484%**, or **5.401-fold**, saving
**18,922,086,559,740 bytes (18.922 TB)** relative to keeping the gzipped FASTQs alone.
Per-sample reductions range from 78.97% to 84.56%, with a median of 81.57%.

## Method

1. Take the processing manifests and match each sample/library ID to its
   processed CRAM and CRAI filenames. Use one archived CRAM per record; older
   alternative processing outputs are not counted again.
2. Sum the actual byte lengths of every original FASTQ.gz listed for each sample,
   including both mates and all lanes, using current filesystem metadata. The
   WES manifests identify 2,000 distinct input files and the WGS manifests identify
   7,024. Check that no source file is reused across records. No uncompressed
   FASTQ estimates or scheduler storage estimates are used.
3. Obtain current byte lengths from directory listings of the archived files.
   Verify that each matched CRAM and index exists and has a positive size. The
   CRAM totals also agree with an independent processed-file inventory.
4. Calculate `reduction = 1 - sum(CRAM + CRAI bytes) / sum(FASTQ.gz bytes)` and
   `fold reduction = sum(FASTQ.gz bytes) / sum(CRAM + CRAI bytes)`. The headline
   uses summed file sizes, not an average of sample ratios. Here, 1 TB = 10^12 bytes.

Only aggregate measurements are included here.

## Encoding and interpretation

The archival encoding command explicitly selects:

```bash
samtools view --output-fmt cram,version=3.1,archive \
  --reference reference.fa --write-index -o archived.cram tagged.bam
```

Encoding settings were identified from the processing workflow. This size
audit did not decrypt the files, inspect
every CRAM header, rerun encoding or perform a new read-content comparison.

The standalone `create_cram.py` currently uses samtools' default CRAM settings;
matching the archival encoding above requires explicitly selecting
`version=3.1,archive`. The measured storage savings apply to the archived files.

File sizes do not establish exact FASTQ recoverability. The stored tags and the
[quality/identity checks](validation.md) determine what can be reconstructed.
Keep the matching reference and validate recovery before choosing a retention
policy. Exact original gzip bytes, read order and FASTQ comments are outside this
format's restoration contract.
