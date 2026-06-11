# tirvish_rs

A faithful Rust port of genometools' `gt tirvish` (GenomeTools 1.6.5) — de-novo
terminal-inverted-repeat (TIR) transposon detection — built to replace the
pathologically slow `gt tirvish` step in the TIR-Learner pipeline on
repeat-rich genomes.

**Independent tool.** May copy code from `grf_rs` (same author) but does not
depend on it or ship with it.

## Status: complete & bit-exact

All five stages are ported and validated against an instrumented `gt tirvish`
(1.6.5) reference. On the committed oracle (four ~5 Mb Pacific white shrimp
chunks, multi-contig), the full pipeline reproduces gt's output **exactly**:

| chunk | elements |
|-------|----------|
| chunk0 | 270 |
| chunk1 | 267 |
| chunk2 | 295 |
| chunk3 | 321 |

**1,153 elements total, every field identical** (`start stop tir1 tir2 tsd1
tsd2 sim`). Each stage was also validated tuple-for-tuple in gt's own internal
coordinates (33,776 seeds → 32,061 extended/scored pairs → final elements).

This is the **faithful** port (single-threaded, no algorithmic shortcuts).
Optimization (banded edit distance, fragment parallelism, alloc reuse) is the
next phase and must stay validated against this same oracle.

## Pipeline (gt source → module)

1. `encode` / `sa` — gt's `-mirrored` encseq layout; SA via libsais over an
   integer alphabet (ACGT + unique sentinels), LCP via Kasai.
2. `maxpairs` — Abouelhoda–Kurtz maximal-pairs (`esa-maxpairs.c` +
   `esa-bottomup-maxpairs.inc`).
3. `seeds` — `gt_tir_store_seeds` (midpos / distance / contig / length filters).
4. `xdrop` — `gt_evalxdroparbitscoresextend` (Zhang–Schwartz–Miller greedy Xdrop).
5. `tsd` — `gt_tir_search_for_TSDs` + best-TSD (brute-force flank MEMs; gt builds
   a per-seed suffix array, this doesn't).
6. `similarity` — `greedyunitedist` (Myers O(nd) unit edit distance) + the gate.
7. `pipeline` — sort + `gt_tir_remove_overlaps` "best" + GFF coordinate output.

## Build / run

```
cargo build --release
./target/release/tirvish <genome.fa>        # gold-TSV: seqid start stop tir1 tir2 tsd1 tsd2 sim
```

Per-stage validators (compare against an instrumented gt's `TIRVISH_TRACE` /
`TIRVISH_XD` dumps): `seedcount`, `xdropcheck`, `tsdcheck`, `simcheck`.

Parameters are locked to the TIR-Learner invocation: `-seed 20 -mintirlen 10
-maxtirlen 1000 -mintirdist 10 -maxtirdist 5000 -similar 80 -mintsd 2 -maxtsd 11
-vic 13` (Xdrop scores at gt defaults: mat 2, mis -2, ins/del -3, xdrop 5).
