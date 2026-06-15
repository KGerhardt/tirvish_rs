# tirvish_rs — performance & algorithmic notes

A change-by-change record of how `tirvish_rs` went from a faithful-but-slow port
of `gt tirvish` to ~6.5× faster single-threaded (and a fragment-parallel batch
driver), with identical output.

All numbers are chunk0 of the oracle fixture (a ~5 Mb multi-FASTA of Pacific white
shrimp contigs — a near-worst case for TIRvish), single-threaded unless noted.
Every change is **bit-exact** against the committed gold sets (4 chunks, 1153
elements, every field). To reproduce the timings and verify identity yourself, see
`testdata/` (`run_compare.sh` runs `gt tirvish` and `tirvish_rs` side by side with
`/usr/bin/time -v` and diffs the predictions).

---

## Preface: Faithfulness constraint and memory usage design

`tirvish_rs` exactly reproduces `gt tirvish`'s raw prediction
multiset. Even if more efficient, arguably better algorithms are available
to complete a given logical step, tirvish-rs does not use them unless they have
the exact same final output. There are therefore two kinds of optimizations in the code:

- **Canonical stages** compute a well-defined mathematical value. *Any* correct
  algorithm yields the same number, so the most performant algorithm is preferred,
  e.g. any algorithm computing Levenshtein edit distance.
- **Heuristic stages** compute an *algorithm-defined* result: the answer depends
  on the specific algorithm's choices (scores, tie-breaks, pruning order). A
  "better" algorithm gives a *different* answer e.g., GenomeTools' greedy X-drop
  extension whose stopping point *is* the output. Optimizations must preserve
  algorithmic choices exactly.

The acceptance test is the oracle in `TIR-Learner/tirlearner_run/tirvish_oracle/`
(committed chunks + `gt tirvish` gold TSVs); per-stage validator bins
(`seedcount`/`xdropcheck`/`tsdcheck`/`simcheck`) check each stage in isolation.

This code was designed for deployment in the Purdue ANVIL supercomputer environment where
thread count and memory are forcibly tied together at 2GB RAM per thread. We targeted
this RAM limit as the acceptable upper cap for total per-thread RAM usage. This code was also
developed specifically for TIR-Learner v4 (https://github.com/KGerhardt/TIR-Learner) 
which uses genomeSplitter (https://github.com/KGerhardt/genomesplitter) to divide a genome into
~5 million base pair fragments with overlaps that allow all TIRvish recoveries to be found 
faithfully, if by parts - any overlap size exceeding the max distance at which TIR elements may
be found will work.

As a result of the expected fixed input size, we pursued some more memory-greedy decisions compared 
to the GenomeTools version of TIRvish, most notably in that we use a fully in-memory suffix array 
structure. Were a user to apply tirvish-rs to a large genome sequence, this could cause a large use of 
RAM. We recommend doing exactly what TIR-Learner does: use genomeSplitter to pre-chunk your genome, THEN
use tirvish-rs.

---

## Summary

Per-stage wall time, `gt` vs `tirvish_rs` (chunk0, single-threaded):

| Stage | `gt tirvish` | `tirvish_rs` (final) | Speedup |
|---|---:|---:|---:|
| 1 — SA + maximal-pair seeding | 32.6 s\* | 10 s | 3.3× |
| 2 — xdrop arm extension | 132.7 s | 48 s | 2.8× |
| 3 — TSD search | 70.6 s | 14 s | 5.0× |
| 4 — similarity (edit distance) | 405.7 s | 17 s | 24× |
| 5 — sort + overlap removal | 1.4 s | 0.7 s | 2× |
| **total** | **~618 s** | **~94 s** | **~6.5×** |

\* `gt` stage 1 = `suffixerator` (~6 s) + `gt_enumeratemaxpairs` (~26.6 s).

Peak RAM use for tirvish-rs was 0.66 GB / thread in these tests, well under our 
targeted 2 GB/thread. Most of this is the suffix array.

---

## Stage 2 — Xdrop arm extension

Xdrop extends each exact seed outward into full TIR arms via a greedy X-drop
gapped alignment (`gt_evalxdroparbitscoresextend`: scores mat +2 / mis −2 / ins −3
/ del −3, drop threshold 5), run twice per seed (left + right) over ~5 M seeds =
~10 M extensions. It is a **heuristic** stage (§0): only constant-factor wins are
legal.

### 2.1 SWAR snake

- **Problem.** The greedy front's "snake" (extending along matching bases) used a
  scalar byte-by-byte loop over the encoded sequence — ~7.4 B character compares.
- **Fix.** Reuse the 2-bit packed genome's SWAR longest-common-extension (32 bases
  per XOR), already used by stage 4 and unit-tested identical to the scalar LCP.
- **Gain.** xdrop **70 s → 65.5 s** (~7%).
- **Why it's better.** LCP is canonical, so SWAR is bit-exact; and it's free
  (reuses existing infrastructure). **But the modest gain was itself the finding**:
  it proved xdrop is *not* snake-bound (the scalar compares were cheap and
  branch-predictable). That redirected effort to the real cost ↓.

### 2.2 Rolling 2-row front buffer

- **Problem.** `gt` indexes the DP front by `frontidx(d,k) = d² + d + k` — a
  triangular array that grows to ~7.7 MB per call, costs a multiply per access
  (~5 per cell × 2.6 B cells), and strides through memory cache-unfriendly.
- **Fix.** For TIRvish's fixed scores the cost-distances are all 1, so every front
  read is from the *previous* generation only. Keep **two rolling rows**
  (`f_cur`, `f_prev`); diagonal `k` of generation `g` lives at local index `g+k`,
  turning the `d²+d+k` multiply into an add and shrinking the hot data to two tiny
  L1/L2-resident rows. The `gt` per-generation fills overwrite every diagonal, so
  no stale slot is ever read. Also dropped the write-only direction array (`gt`'s
  traceback aid, never read — we only need the endpoint).
- **Gain.** xdrop **65.5 s → 48.2 s** (~26%).
- **Why it's better.** xdrop turned out to be **bookkeeping/cache-bound**, not
  compute- or snake-bound. Removing the per-cell multiply and collapsing a 7.7 MB
  array to a couple of cache-resident rows attacks exactly that. Pure
  representation change → bit-exact (guarded by a debug-assert on the
  cost-distance invariant).

---

## Stage 3 — TSD search

### 3.1 Brute-force flank matching

- **Problem.** `gt` builds a per-seed enhanced suffix array over the tiny ±13 bp
  TSD flanks (`gt_sarrquerysubstringmatch`) — heavy machinery for a ~26 bp window.
- **Fix.** A direct brute-force maximal-match scan over the flanks (done during the
  original port).
- **Gain.** ~3.7× on stage 3 vs the ESA build.
- **Why it's better.** On windows this small, the asymptotic structure of a suffix
  array is pure overhead; a linear scan wins on constants. Output-identical
  (validated tuple-for-tuple against `gt`'s trace, including the SA tie-break order).

### 3.2 Thread-local scratch

- **Problem.** The flank scan allocated a small `Vec` (plus a per-offset match
  list) on every one of ~4.8 M calls — heavy allocator churn.
- **Fix.** Reuse per-thread scratch buffers (cleared, not reallocated), the same
  pattern as the xdrop/greedy front buffers.
- **Gain.** stage 3 **19 s → 14.9 s** (~22%).
- **Why it's better.** Eliminates per-seed heap traffic at the source. (Notably it
  did **not** reduce peak RSS — confirming the loop's memory "creep" was live
  candidate data, not allocator fragmentation. See §M.)

---

## Stage 4 — Similarity (Levenshtein):

Stage 4 was 66% of `gt`'s runtime. It computes, for each candidate's two TIR
arms, `sim = 100·(1 − edist/max(ulen,vlen))` and drops anything below 80%.

### 4.1 Banded edit distance

- **Problem.** `gt`'s `greedyunitedist` computes the *full* edit distance for
  every candidate, including the ~20% that fail the 80% gate — grinding out a
  large distance just to reject it.
- **Fix.** Cap the DP band at `max/5 + 2`. A passer has `edist ≤ 0.2·max`, so it
  stays exact; a failer bails the moment it exceeds the band.
- **Gain.** Removed most of the wasted work on failers; the recorded similarity
  of passers is unchanged.
- **Why it's better.** The gate only needs the *exact* value for passers; failers
  only need "below threshold." The `+2` margin guarantees no true passer is ever
  cut. (Same trick as `grf_rs::seq_complexity` bailing once the boolean is decided.)

### 4.2 Bit-parallel Levenshtein via `rapidfuzz` — the big one

- **Problem.** Even banded, our hand-written greedy/diagonal front DP is *scalar*:
  one DP cell per loop iteration. Similarity remained the largest stage.
- **Fix.** We dumped the **actual** arm-pair corpus (4.34 M pairs from chunk0 via
  the `simdump` bin) and benchmarked it against `strsim`, `triple_accel`,
  `stringzilla`, and `rapidfuzz`. `rapidfuzz`'s block-based **bit-parallel Myers**
  with a `score_cutoff` (= our band) won decisively, and produced a *bit-identical*
  banded result (matching checksum over all 4.34 M pairs). Adopted it in
  `compute_similarity`. The other crates' unbounded SIMD/scalar paths were slower
  and did not finish within the benchmark deadline.
- **Gain.** Similarity **89 s → ~15 s (≈5.7×)**.
- **Why it's better.** Bit-parallel Myers evaluates 64 DP cells per machine word;
  our scalar front did one per iteration. Since similarity is a *canonical*
  quantity (§0), swapping the algorithm is faithful by construction — the matching
  checksum and the oracle confirm it. (This **corrected an earlier false
  conclusion** that "Myers is slower here": that was a naïve single-word Myers; a
  production block-banded bit-parallel implementation wins by a wide margin.)

### 4.3 Reused arm buffers

- **Problem.** Feeding `rapidfuzz` an iterator that extracts each base from the
  2-bit packed genome on the fly (`base_at`) compiled worse than handing it a
  contiguous slice.
- **Fix.** Materialize each arm once into a per-thread reused `Vec<u8>` (a tight
  monomorphic fill loop), then pass `rapidfuzz` a slice.
- **Gain.** Similarity **20 s → 17.4 s**.
- **Why it's better.** Same number of `base_at` calls, but the bit-extraction
  optimizes far better in a dedicated loop than threaded through a generic
  iterator. The buffer is reused (≈1 KB/thread) — no per-pair allocation and, by
  design, **no** materialization of the whole corpus (which would blow the RAM
  budget; fine for a flat benchmark, wrong for production).

---

## Memory

Goal: stay comfortably under ~2 GB/thread on the worst genomes so many fragments
pack onto one node.

### M.1 Prune skipped candidates before the overlap sort

- **Problem.** `gt` keeps *every* length-passing candidate in one array for the
  final overlap sort. On chunk0 that's ~4.6 M `TirPair`s (~500 MB), even though
  ~270 survive.
- **Fix.** Retain only the *live* (non-skipped) candidates, plus the single
  global-first one. In `gt_tir_remove_overlaps`, every skipped candidate at index
  ≥ 1 is `continue`d before it can affect anything — it's inert. The *only* skipped
  candidate that influences output is the one at index 0 (it seeds the reference
  range unconditionally), so we carry just that one.
- **Gain.** Part of the 1.23 → 1.04 GB peak drop; shrinks the stage-5 sort from
  millions to hundreds.
- **Why it's better.** Output-identical by analysis of the overlap algorithm, and
  confirmed by the oracle. (A subtlety: the live set is still large — only ~20% of
  length-passers are skipped, because ~80% *pass* similarity and it's *overlap
  removal*, not the similarity gate, that collapses 3.67 M → 270.)

### M.2 Lean candidate projection

- **Problem.** The retained candidate is a 104-byte `TirPair` carrying compute-only
  fields (seed positions, mirror coords) that stage 5 never reads, multiplied by
  ~3.67 M live candidates.
- **Fix.** Project each kept candidate to a 64-byte `LeanPair` with only the ~9
  fields the sort/overlap/emission need.
- **Gain.** ~176 MB off the peak.
- **Why it's better.** 64 vs 104 bytes × millions is the dominant memory term once
  pruning (§M.1) is in place. `process_seed` still computes the full pair (stages
  2–4 need the dropped fields); we project at collection.

### M.3 Free dead SA scaffolding eagerly

- **Problem.** The suffix array, LCP table, and per-position helper arrays live for
  the whole run, but nothing in the per-seed loop reads them — they're dead once
  seeds are built. (`gt`'s i32 SA/LCP were also kept as a redundant copy of the
  u64 versions.)
- **Fix.** Scope/`mem::take` each piece the moment its last reader finishes:
  `sa_input` after SA build, `suftab`/`lcptab` after seeding, `seqnum_of` after
  `store_seed`, and drop the i32 SA/LCP duplicate.
- **Gain.** ~164 MB off the loop's footprint; overall peak **1.23 → 0.66 GB**.
- **Why it's better.** Nothing reads them later, so freeing is invisible to output;
  it just stops them occupying the high-water mark during the long per-seed loop.

---

## Parallelism

### P.1 Fragment-level batch driver with tail steal

- **Problem.** Parallelism was at the wrong level: only an inner per-seed
  `par_iter` existed. The real workload is many pre-batched 5 Mb fragment *files*;
  the only ways to run them were sequential single-file (wastes cores on each
  fragment's serial stage 1) or N processes (oversubscription).
- **Fix.** Mirror `grf_rs`'s structure: one global rayon pool; an outer
  `par_iter` over fragments (1/thread) with per-fragment output emission; the inner
  per-seed `par_iter` nests in the same pool and is coarsened
  (`with_min_len(8192)`) so it's near-inert in bulk. Result: **bulk =
  inter-fragment** (each worker owns a fragment), **tail = intra-fragment** (when
  remaining + in-flight fragments < workers, idle workers steal a straggler's
  seed chunks). CLI: `tirvish --batch <outdir> [--threads N] <paths…|stdin>`.
- **Gain.** ~2.7× faster than sequential single-file for the multi-fragment case;
  ~0.6 GB per concurrent fragment; tail-steal validated (idle workers pick up a
  straggler's seed chunks).
- **Why it's better.** It matches the actual deployment of 5 Mbp genome chunks,
- removes the per-fragment serial-stage-1 bottleneck (fragments overlap each other's stage 1), avoids
  oversubscription, and keeps RAM per fragment bounded.

---

## Ingestion

### needletail

- **Problem.** Five hand-rolled FASTA readers across the bins; no gzip/multi-line
  robustness.
- **Fix.** One `read_fasta` in `lib.rs` using `needletail`; all bins call it.
- **Gain.** Faster parse, gzip + multi-line support, less duplicated code.
- **Why it's better.** A maintained, fast parser; byte-identical records on the
  oracle. Pinned to `grf_rs`'s version for a shared resolved dependency set.

---

## Investigations that did **not** pan out (and why that matters)

Negative results, kept because they justify *not* doing the obvious things.

### Hand-rolled bit-parallel Myers for similarity → reverted (was slower)
A single-word Myers lost to the output-sensitive greedy front, which led to the
(wrong) belief that bit-parallel doesn't help here. The crate benchmark (§4.2)
overturned this: a *production* block-banded bit-parallel implementation wins ~5.7×.
**Lesson:** benchmark real libraries on the real corpus before concluding.

### Seed / candidate redundancy → refuted
Hypothesis: `gt` over-enumerates (e.g. a 100-base run as ~80 length-20 seeds), or
downstream work is mostly duplicate, so seeds could be dropped or memoized.
Measured: seed lengths have **no** length-20 spike (only 7.7%, mean 32, max 257) —
`gt` emits genuine *maximal* pairs, so a run collapses to one long seed. And the
arm-pair corpus is only **1.6× redundant** (61% unique). So a content cache is a
wash (≈1 GB cache or a hash whose cost ≈ the 4 µs edit distance it replaces), and
the seeds are real repeat content, not waste. **No cheap seed-drop exists**; the
5 M seeds reflect genuine shrimp repeat density.

### 2-bit/XOR SWAR LCE inside the similarity front → no speedup (kept as foundation)
The matching runs in similarity are short, so the SWAR window overhead doesn't
amortize there. (It *is* used in xdrop §2.1 and as the LCE substrate generally.)

### Wider-SIMD / banding variants on the similarity front → lost to output-sensitivity
The pairs reaching similarity are moderate-length, moderate-similarity, which puts
the output-sensitive greedy front in its sweet spot; regular-DP × hardware-width
attacks did more total cell work than they saved. (Superseded anyway by §4.2.)

---

## Open directions

- **xdrop is a heuristic stage** (§0), so any speedup must reproduce its exact
  result — it can't be algorithm-swapped the way similarity was. SIMD on the cell
  eval and further cache/branch work are open.
- **Reducing the seed count** (a bounded MEM finder) would hit stages 2–4 at once.
  The seeds are real maximal pairs (above), so dropping any requires proving it
  never produces a surviving element.
- **PyO3 entry** into TIR-Learner's `tirvish_new.py` would reuse the `run_batch`
  logic over in-memory fragments (as `grf_rs::do_rusty_grf` does).
