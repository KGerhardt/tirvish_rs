# testdata — identity & timing harness

Four ~5 Mb Pacific white shrimp fragments (a near-worst case for TIRvish) and a
script that runs the original `gt tirvish` and `tirvish_rs` on each, times both
with `/usr/bin/time -v`, and diffs the predictions to show they are identical.

## Contents

- `chunks.tar.gz` — `chunk0.fa` … `chunk3.fa` (multi-FASTA fragments; the script
  extracts these to a temp dir).
- `run_compare.sh` — runs both tools per chunk, times each call, diffs outputs.
- `gff_to_tsv.py` — converts `gt tirvish` GFF3 into the exact TSV shape `tirvish_rs`
  emits, so the two can be diffed directly. Each row carries the six `(start,stop)`
  coordinate pairs gt writes per element — `full`, `TSD1`, `body`, `TIR1`, `TIR2`,
  `TSD2` — plus `sim`; the diff compares **every** column.
- `expected_candidates.tar.gz` — `chunk{0..3}.tirvish.tsv`, the reference output in
  that shape (handy without running `gt`).

## Usage

```bash
# needs genometools `gt` on PATH and python3; builds tirvish_rs if not built
./run_compare.sh

# or point at specific binaries
GT=/opt/genometools/bin/gt RS=../target/release/tirvish ./run_compare.sh
```

Per chunk it prints `OUTPUTS IDENTICAL: N predictions match` (or a diff), plus the
wall time and peak RSS of `gt suffixerator`, `gt tirvish`, and `tirvish_rs`. Exit
status is non-zero if any chunk's predictions differ.

The `gt` flags are exactly those TIR-Learner's `tirvish_new.one_tirvish` uses
(`suffixerator -tis -suf -lcp -des -ssp -sds -dna -mirrored`, then `tirvish -seed
20 -mintirlen 10 -maxtirlen 1000 -mintirdist 10 -maxtirdist 5000 -similar 80
-mintsd 2 -maxtsd 11 -vic 13 -seqids yes`).
