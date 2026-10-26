# Troubleshooting

| Symptom | Action |
| --- | --- |
| Missing `samtools`, `pigz` or `bam_merge` | Activate the environment and run `make`; check `PATH` |
| `zlib.h` missing while building | Install zlib development headers, or use the supplied conda environment |
| `Python.h` missing for native checks | Install Python development headers; the primary CLI does not need this extension |
| Reference mismatch or CRAM decode failure | Supply the exact encoding reference, check contig names/lengths and run `samtools faidx` |
| Output directory already exists | Choose a new directory; creation deliberately refuses replacement |
| Very high annotation RAM | Use QNAME-grouped aligner output close to source FASTQ order with `--name-collated`; work per read group |
| Temporary storage exhausted | Set `--tmp-dir` to sufficient scratch space and tune samtools sort memory/threads |
| Quality mismatch in strict mode | Inspect aligner quality changes; pipeline mode permits only the documented `#`/`!` and N-base differences |
| Sequence/read-count mismatch | Check original FASTQs, trimming direction, missing mates, duplicate QNAMEs and input read group |
| Hard clip without matching tag payload | The removed bases are absent; use the correct original tagged file or regenerate tags from original FASTQs |
| Reversion reports multiple/missing primary mates | Inspect `(QNAME, RG)` identities and mate flags; singletons and single-end inputs are unsupported |
| Native checks differ on legacy qualities | Native Phred detection is heuristic; normalize encoding explicitly and use the standalone checker |

Check the process exit code for every command. `samtools quickcheck` alone does
not read all records; the creation wrapper additionally decodes the complete
CRAM. For hand-written shell pipelines enable `set -euo pipefail` and run the
final comparison explicitly.

Reversion collates internally unless `--name-collated` is supplied. Set `TMPDIR`
to a suitable directory for its samtools temporary files. The two FASTQ files and
stats are written to temporary files beside their destinations before replacing
them. Replacement is atomic per file, not a transaction spanning all three
files. Use new output names when retaining an earlier restoration.
