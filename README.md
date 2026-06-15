# tirvish_rs

A faithful Rust port of genometools' `gt tirvish` (GenomeTools 1.6.5) — de-novo
terminal-inverted-repeat (TIR) transposon detection — built to replace the
sometimes slow `gt tirvish` step in the TIR-Learner pipeline on
repeat-rich genomes.

## Complete & bit-exact

`tirvish-rs` has been validated against an instrumented `gt tirvish`
(1.6.5) reference. On the committed oracle (four ~5 Mb Pacific white shrimp
chunks, multi-contig), the full pipeline reproduces gt's output **exactly**:

**1,153 elements total, every field identical** (`start stop tir1 tir2 tsd1
tsd2 sim`). Each stage was also validated tuple-for-tuple in gt's own internal
coordinates (33,776 seeds → 32,061 extended/scored pairs → final elements).

## Performance (per-stage, chunk0, single-threaded)

Wall seconds for the faithful port vs the instrumented gt 1.6.5 reference, on
chunk0 (~5.25 Mb, ~5.1 M seeds). Already ~2.7× faster end-to-end **before any
optimization** — from the Rust constant factor plus a brute-force TSD search
(gt builds a suffix array per seed; this doesn't).

| stage            | gt                                     | tirvish_rs | speedup |
|------------------|----------------------------------------|------------|---------|
| stage 1 (seeds)  | ~32.6 (suffixerator 6 + maxpairs 26.6) | 10.6       | ~3×     |
| 2 — Xdrop        | 132.7                                  | 98.5       | 1.35×   |
| 3 — TSD          | 70.6                                   | 19.2       | 3.7×    |
| 4 — similarity   | 405.7                                  | 93.9       | 4.3×    |
| 5 — sort+overlap | 1.4                                    | 1.0        | 1.4×    |
| **total**        | **~611**                               | **224**    | **2.7×**|

Note the dominant stage moved: in gt, similarity is ~66%; in tirvish_rs, Xdrop
(98.5 s) and similarity (93.9 s) are co-dominant, so optimization weights both.
(`TIRVISH_RS_TIME=1 ./target/release/tirvish <fa>` prints this breakdown.)

## Performance (overall, all chunks, single-threaded)

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
