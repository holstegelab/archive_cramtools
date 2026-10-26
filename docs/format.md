# Storage and restoration contract

## Tags

| Tag | SAM type | Meaning in this implementation |
| --- | --- | --- |
| `YB` | `Z` | Bases preceding `SEQ` in alignment orientation |
| `YQ` | `Z` | Corresponding ASCII qualities |
| `ZB` | `Z` | Bases following `SEQ` in alignment orientation |
| `ZQ` | `Z` | Corresponding ASCII qualities |
| `XT` | `i` | Number of bases added back as trimming annotations by `bam_merge` |
| `RG` | `Z` | Read-group identity, retained during restoration |

These are local conventions. The SAM specification reserves tags beginning with
X, Y and Z for local use; another tool can assign different meanings to the same
names. They are not a universal CRAM extension. See the
[SAM optional-fields specification](https://samtools.github.io/hts-specs/SAMtags.pdf).
Do not remove, repurpose or annotate these restoration tags twice.

For a primary alignment:

```text
aligned_sequence = YB + SEQ  + ZB
aligned_quality  = YQ + QUAL + ZQ
```

Missing tag pairs contribute empty strings. If FLAG contains `0x10`, reverse
complement the combined sequence and reverse the combined quality string to get
the original read orientation. The lab convention treats untagged unmapped reads
as already in original orientation; tagged unmapped reads follow the restoration
path. This is an inherited convention, not a claim about every possible SAM
producer. The reverter checks that tag pairs and their lengths are consistent,
and that mapped hard clips have the corresponding tag payload.

For example, an original forward read `ACGTTT` with qualities `IIII##`, aligned
after trimming as `ACGT/IIII`, becomes:

```text
CIGAR: 4M2H    SEQ: ACGT    QUAL: IIII    ZB:Z:TT    ZQ:Z:##
```

An original reverse-aligned read `AACCGG/IIII##` aligned as `GGTT/IIII` becomes:

```text
CIGAR: 2H4M    SEQ: GGTT    QUAL: IIII    YB:Z:CC    YQ:Z:##
```

Restoration combines `CCGGTT/##IIII`, then reverses orientation to recover
`AACCGG/IIII##`. The side of a tag describes alignment orientation, not R1/R2.

## Annotation requirements

The extracted `bam_merge.c` associates primary alignments with original paired
FASTQs by name. After removing a terminal `/1` or `/2`, pair names must be unique
within an annotation invocation. FASTQ comments after whitespace are ignored.
The C code is specialized for reads shortened at their original **3′ end**.
The retained sequence is a prefix of the original read; on reverse alignments it
is a suffix after reverse complementation. Arbitrary 5′ trimming, middle excision,
primer processing or sequence correction needs a different annotation strategy.

For mapped records, newly recovered material is represented by hard clipping in
the CIGAR and stored on the primary record. Soft-clipped bases already reside in
`SEQ` and need no extra storage. Supplementary CIGAR hard clips are extended as
needed; restoration uses the primary payload only. Secondary/supplementary
records do not become additional FASTQ reads.

Paired FASTQs absent from the entire alignment are re-added as paired unmapped
records, with flags 589 and 653 (including the vendor-quality-failure bit).
They inherit the last primary read group's ID seen by the annotator; if there
was no such record they have no RG tag. This is why creation is limited to one
read group. This recovery does not repair an alignment missing just one mate;
the default final comparison will reject lost payloads.

The creation wrapper rejects pre-existing restoration tags and multiple declared
read groups. It expects primary pairs with sequence and quality present. The
final decoded comparison runs by default, including when the C annotator itself
does not detect a discrepancy. `--fast` skips the separate FASTQ checksum pass
and this comparison, while retaining annotation, sorting, quickcheck and indexing.
It records verification as skipped. Input requirements still apply: duplicate
names, unsupported trimming and missing source payloads should be corrected
upstream.

## What is retained

The CRAM contains alignment records, stored primary sequence/quality, restoration
tags and read-group headers. Ordinary samtools operations can read the alignment.
To reconstruct our hard-clipped data use `bam_revert.py`: ordinary `samtools
fastq` does not reconstruct these custom tags. Its documented handling of
hard-clipped data is described in the
[samtools FASTQ manual](https://www.htslib.org/doc/samtools-fasta.html).

Original FASTQ ordering, comments, plus-line annotations and compressed bytes are
not retained. No original pre-processing BAM is reconstructed. Quality values
already overwritten by an aligner cannot be recovered from these tags when that
base remained in `SEQ`. This package does not implement encryption, alignment,
duplicate marking or a general FASTQ repair pipeline.

The matching reference sequence is required to decode the reference-based CRAM.
Retain the exact FASTA and index, not just a genome-build label. CRAM encoding and
reference handling are supplied by [samtools/HTSlib](https://www.htslib.org/doc/samtools.html).
