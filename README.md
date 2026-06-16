# tirvish_rs

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](https://bioconda.github.io/recipes/tirvish-rs/README.html)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/tirvish-rs/badges/version.svg)](https://anaconda.org/bioconda/tirvish-rs)

A faithful Rust port of genometools' `gt tirvish` (GenomeTools 1.6.5) — de-novo
terminal-inverted-repeat (TIR) transposon detection — built to replace the
sometimes slow `gt tirvish` step in the TIR-Learner pipeline on
repeat-rich genomes.

## Install

Available from [Bioconda](https://bioconda.github.io/recipes/tirvish-rs/README.html):

```
conda install -c bioconda tirvish-rs
```

This installs the `tirvish` binary (and a `tirvish_rs` alias). You can also build
from source with `cargo build --release` (see [Build / run](#build--run)).

## Complete & bit-exact

`tirvish-rs` has been validated against an instrumented `gt tirvish`
(1.6.5) reference. On the committed oracle (four ~5 Mb Pacific white shrimp
chunks, multi-contig), the full pipeline reproduces gt's output **exactly**:

**1,153 elements total, every field identical**. Each hit is emitted as one row
carrying the six `(start, stop)` coordinate pairs gt writes per element — full
element, TSD1, body, TIR1, TIR2, TSD2 — plus `tir_similarity` (the exact data gt
puts in its GFF, minus the line-bucketed/regex parsing). Each stage was also
validated tuple-for-tuple in gt's own internal coordinates (33,776 seeds →
32,061 extended/scored pairs → final elements).

The exact expected outputs are committed in `testdata/expected_candidates.tar.gz`
(the four per-chunk candidate TSVs); `testdata/run_compare.sh` regenerates them
with `gt tirvish` + `tirvish_rs` side by side, times both, and diffs them.

## Performance (per-stage, chunk0, single-threaded)

Wall seconds for the final tirvish_rs vs the instrumented gt 1.6.5 reference, on
chunk0 (~5.25 Mb, ~5.1 M seeds). See `OPTIMIZATIONS.md` for the change-by-change
account of how each stage got here.

| stage            | gt                                     | tirvish_rs | speedup |
|------------------|----------------------------------------|------------|---------|
| stage 1 (seeds)  | ~32.6 (suffixerator 6 + maxpairs 26.6) | ~10        | ~3×     |
| 2 — Xdrop        | 132.7                                  | 48         | 2.8×    |
| 3 — TSD          | 70.6                                   | 14         | 5×      |
| 4 — similarity   | 405.7                                  | 17         | 24×     |
| 5 — sort+overlap | 1.4                                    | 0.7        | 2×      |
| **total**        | **~618**                               | **~94**    | **~6.5×**|

The dominant stage moved twice: in gt, similarity is ~66%; a bit-parallel banded
edit distance cut it to ~17 s, leaving Xdrop (~48 s) as the largest stage. (These
per-stage figures were measured during development; the env-gated timing
instrumentation has since been removed from the production build.)

## Performance (overall, all chunks, single-threaded)

End-to-end wall time per oracle chunk — gt = `suffixerator` (~6 s) + `tirvish`
search, vs tirvish_rs single-threaded. Consistent ~6.2–6.8× across all four:

| chunk   | gt (s)    | tirvish_rs (s) | speedup |
|---------|----------:|---------------:|--------:|
| chunk0  | ~618      | 94.9           | 6.5×    |
| chunk1  | ~743      | 119.1          | 6.2×    |
| chunk2  | ~544      | 80.5           | 6.8×    |
| chunk3  | ~521      | 78.2           | 6.7×    |
| **all** | **~2425** | **372.7**      | **6.5×**|

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
7. `pipeline` — sort + `gt_tir_remove_overlaps` "best" + TSV coordinate output
   (the six per-element `(start, stop)` pairs gt would emit as GFF features).

## Build / run

```
cargo build --release

# single fragment -> stdout (TSV: seqid + the six per-element coord pairs + sim)
./target/release/tirvish <genome.fa>

# batch of pre-chunked fragments: one <basename>.tirvish.tsv per fragment,
# parallel at the fragment level (a straggler's seeds get stolen at the tail)
./target/release/tirvish --batch <outdir> [--threads N] <frag1.fa> <frag2.fa> ...
./target/release/tirvish --batch <outdir> [--threads N]      # paths from stdin
```

It operates directly on a FASTA — no pre-built mirrored index needed. All gt
tirvish parameters are accepted as options; defaults are the TIR-Learner
invocation: `-seed 20 -mintirlen 10 -maxtirlen 1000 -mintirdist 10 -maxtirdist
5000 -similar 80 -mintsd 2 -maxtsd 11 -vic 13` (Xdrop scores `-mat 2 -mis -2
-ins -3 -del -3 -xdrop 5`).

Per-stage validators (compare against an instrumented gt's `TIRVISH_TRACE` /
`TIRVISH_XD` dumps): `seedcount`, `xdropcheck`, `tsdcheck`, `simcheck`.
