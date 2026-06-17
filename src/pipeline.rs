//! Full de-novo TIRvish pipeline: stages 1-5 -> final element list (gold format).
//!
//! Stage 5 = sort by gt_tir_compare_TIRs (contignumber, left_tir_start,
//! right_transformed_start), gt_tir_remove_overlaps "best" (keep max-similarity
//! within each overlapping cluster), tir_compactboundaries (drop skipped), then
//! emit GFF coordinates (tir_stream.c:940-1001) in the gold-TSV shape.

use crate::encode::{encode, ALPHA};
use crate::maxpairs::enumerate_maxpairs;
use crate::params::Params;
use crate::sa::sa_lcp;
use crate::seeds::{store_seed, Seed};
use crate::similarity::{compute_similarity, double_smaller};
use crate::tsd::{build_pair, search_for_tsds, TirPair};
use crate::encode::Encoded;
use crate::xdrop::{calc_distances, extend_seed, ArbitraryDistances, ArbitraryScores};
use rayon::prelude::*;
use std::cmp::Ordering;

/// One TIRvish hit, carried as the SIX (start, stop) coordinate pairs that gt
/// tirvish emits per element and that TIR-Learner's one_tirvish ingests as
/// `next_result[0..5]` — full element, TSD1, body (no-TSD), TIR1, TIR2, TSD2, in
/// that order. We emit exactly these position data and nothing derived: the
/// consumer rebuilds `next_result` from the columns with plain splits and runs
/// the identical length/size calculations it ran on the GFF. Coordinates are
/// 1-based inclusive, local to the contig, exactly as gt writes them. TIR1/TIR2
/// are the two arms in gt's `(start, end)` GFF order (TIR1 = the smaller).
/// `sim` (gt's tir_similarity) is the lone non-position field, kept as a trailing
/// diagnostic for the gold diffs; the consumer ignores it.
#[derive(Debug, Clone)]
pub struct Element {
    pub seqid: String,
    pub full_start: u64,
    pub full_stop: u64,
    pub tsd1_start: u64,
    pub tsd1_stop: u64,
    pub body_start: u64,
    pub body_stop: u64,
    pub tir1_start: u64,
    pub tir1_stop: u64,
    pub tir2_start: u64,
    pub tir2_stop: u64,
    pub tsd2_start: u64,
    pub tsd2_stop: u64,
    pub sim: f64,
}

/// Final elements as the gold-shape TSV (header + rows), sorted by (seqid, start)
/// like parse_tirvish. Shared by the single-file CLI and run_batch's per-fragment
/// emission. Sorts `els` in place.
pub fn elements_tsv(els: &mut [Element]) -> String {
    els.sort_by(|a, b| a.seqid.cmp(&b.seqid).then(a.full_start.cmp(&b.full_start)));
    let mut s = String::with_capacity(64 + els.len() * 80);
    s.push_str(
        "seqid\tfull_start\tfull_stop\ttsd1_start\ttsd1_stop\tbody_start\tbody_stop\t\
         tir1_start\ttir1_stop\ttir2_start\ttir2_stop\ttsd2_start\ttsd2_stop\tsim\n",
    );
    for el in els.iter() {
        s.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\n",
            el.seqid,
            el.full_start,
            el.full_stop,
            el.tsd1_start,
            el.tsd1_stop,
            el.body_start,
            el.body_stop,
            el.tir1_start,
            el.tir1_stop,
            el.tir2_start,
            el.tir2_stop,
            el.tsd2_start,
            el.tsd2_stop,
            el.sim
        ));
    }
    s
}

/// Inner per-seed parallelism granularity floor. The seed loop is the FRAGMENT's
/// internal parallelism; with fragment-level (outer) parallelism saturating the
/// pool in bulk, this keeps the inner split coarse (near-inert) and only lets
/// idle workers steal seeds in ~SEED_PAR_BLOCK-sized chunks at the tail.
const SEED_PAR_BLOCK: usize = 8192;

/// Compact projection of a TirPair retained for stage 5 (sort + overlap removal +
/// emission). Drops the compute-only fields (pos1/pos2/seed_len and the mirror
/// coords right_tir_start/end) that stages 2-4 needed but stage 5 never reads.
/// The coordinate fields are u32: every position is into the mirrored text, which
/// is < 2^31 (the 32-bit suffix-array limit guarded in `run`), so they provably
/// fit. ~40 B (was 64); with ~3.67M live candidates on a shrimp chunk this is the
/// dominant memory term. `seed_idx` carries the stable-sort tie-break.
#[derive(Clone, Copy)]
struct LeanPair {
    seed_idx: u32,
    contignumber: u32,
    left_tir_start: u32,
    left_tir_end: u32,
    right_transformed_start: u32,
    right_transformed_end: u32,
    tsd_length: u32,
    similarity: f64,
    skip: bool,
}

impl LeanPair {
    fn from_pair(seed_idx: u32, p: &TirPair) -> Self {
        // Positions are < 2^31 (mirrored text length is guarded < i32::MAX), so the
        // u64 -> u32 narrowing is lossless.
        LeanPair {
            seed_idx,
            contignumber: p.contignumber,
            left_tir_start: p.left_tir_start as u32,
            left_tir_end: p.left_tir_end as u32,
            right_transformed_start: p.right_transformed_start as u32,
            right_transformed_end: p.right_transformed_end as u32,
            tsd_length: p.tsd_length as u32,
            similarity: p.similarity,
            skip: p.skip,
        }
    }
}

/// gt_tir_compare_TIRs.
fn compare_tirs(a: &LeanPair, b: &LeanPair) -> Ordering {
    a.contignumber
        .cmp(&b.contignumber)
        .then(a.left_tir_start.cmp(&b.left_tir_start))
        .then(a.right_transformed_start.cmp(&b.right_transformed_start))
}

/// Incremental, bounded-memory driver for stage 5 (sort + overlap-removal + emit).
/// Live candidates are PUSHED roughly in `pos1` order (bounded disorder, see `run`);
/// a sliding window holds the not-yet-safe tail, and `drain_to(watermark)` flushes
/// every buffered candidate with `left_tir_start < watermark - disorder`, in
/// `(compare_tirs, seed_idx)` order, through the same cluster sweep the old
/// `remove_overlaps_pruned` did. `p0` (the global argmin over ALL length-passers,
/// possibly skipped) is tracked as a running min via `observe_p0`; it is provably
/// settled before the first drain emits anything, so the sweep core is seeded
/// lazily on first emission. Output is identical to the whole-array sort +
/// remove_overlaps + emit; the candidate set never fully resides.
struct OverlapSweep<'a> {
    disorder: u32,
    win: Vec<LeanPair>,    // window buffer (unsorted between drains)
    p0: Option<LeanPair>,  // running global argmin over all length-passers
    seeded: bool,          // sweep core seeded (p0 settled) yet?
    ref_start: u32,
    ref_end: u32,
    best: Option<LeanPair>, // current cluster's running max-similarity survivor
    out: Vec<Element>,
    e: &'a Encoded,
    contigs: &'a [(String, Vec<u8>)],
}

/// Build the output Element for one survivor pair: the six (start, stop)
/// coordinate pairs gt emits per element, 1-based and local to the contig.
/// Extracted so the streaming overlap sweep can emit survivors directly.
fn make_element(pair: &LeanPair, e: &Encoded, contigs: &[(String, Vec<u8>)]) -> Element {
    let seqstart = e.fwd_seqstart[pair.contignumber as usize];
    // Emit the sequence id verbatim — never mutate identifiers the input carries.
    let name = contigs[pair.contignumber as usize].0.clone();
    // Widen the u32 coords back to u64 for the offset arithmetic against seqstart.
    let lts = pair.left_tir_start as u64;
    let lte = pair.left_tir_end as u64;
    let rts = pair.right_transformed_start as u64;
    let rte = pair.right_transformed_end as u64;
    let tsd = pair.tsd_length as u64;
    // The full element is anchored on the left-arm start and the right-arm end
    // (± one TSD); the body (no-TSD element) sits one TSD inside each end; the two
    // TSDs flank the body.
    let full_start = lts - seqstart - tsd + 1;
    let full_stop = rte - seqstart + tsd + 1;
    let body_start = full_start + tsd;
    let body_stop = full_stop - tsd;
    let tsd1_start = full_start;
    let tsd1_stop = body_start - 1;
    let tsd2_start = body_stop + 1;
    let tsd2_stop = full_stop;
    // The two TIR arms, each in its own coordinates. gt emits the two TIR features
    // sorted by GtRange (start, then end), so TIR1 = the (start,end)-first arm.
    // Post-TSD the transformed right arm can start at or before the left arm (the
    // arms cross over when long), so sort rather than assume left-then-right.
    let left_arm = (lts - seqstart + 1, lte - seqstart + 1);
    let right_arm = (rts - seqstart + 1, rte - seqstart + 1);
    let (tir1, tir2) = if left_arm <= right_arm {
        (left_arm, right_arm)
    } else {
        (right_arm, left_arm)
    };
    Element {
        seqid: name,
        full_start,
        full_stop,
        tsd1_start,
        tsd1_stop,
        body_start,
        body_stop,
        tir1_start: tir1.0,
        tir1_stop: tir1.1,
        tir2_start: tir2.0,
        tir2_stop: tir2.1,
        tsd2_start,
        tsd2_stop,
        sim: pair.similarity,
    }
}

impl<'a> OverlapSweep<'a> {
    fn new(e: &'a Encoded, contigs: &'a [(String, Vec<u8>)], disorder: u32) -> Self {
        OverlapSweep {
            disorder, win: Vec::new(), p0: None, seeded: false,
            ref_start: 0, ref_end: 0, best: None, out: Vec::new(), e, contigs,
        }
    }

    /// Fold a length-passer (skipped or not) into the running global p0 = gt's pairs[0].
    fn observe_p0(&mut self, c: LeanPair) {
        let beats = match &self.p0 {
            None => true,
            Some(pp) => compare_tirs(&c, pp).then(c.seed_idx.cmp(&pp.seed_idx)).is_lt(),
        };
        if beats {
            self.p0 = Some(c);
        }
    }

    /// Buffer a live candidate in the window.
    fn push(&mut self, c: LeanPair) {
        self.win.push(c);
    }

    /// Flush every buffered candidate with `left_tir_start < watermark - disorder`
    /// (safe: no future candidate, all with pos1 >= watermark, can have a smaller
    /// left_tir_start), in (compare_tirs, seed_idx) order, through the cluster sweep.
    fn drain_to(&mut self, watermark: u32) {
        if self.win.is_empty() {
            return;
        }
        self.win.sort_by(|a, b| compare_tirs(a, b).then(a.seed_idx.cmp(&b.seed_idx)));
        let thresh = watermark.saturating_sub(self.disorder);
        // left_tir_start is a global position, so it is monotonic in the
        // (contig, left_tir_start) sort order -> a binary partition is valid.
        let split = self.win.partition_point(|c| c.left_tir_start < thresh);
        let drained: Vec<LeanPair> = self.win.drain(..split).collect();
        for c in drained {
            self.feed_core(c);
        }
    }

    /// One candidate through the cluster sweep. Seeds the core lazily on the first
    /// call, when p0 is provably settled (it has the global-min left_tir_start, so it
    /// is final before anything below the first watermark is emittable). Identical
    /// logic to the former remove_overlaps_pruned + emit.
    fn feed_core(&mut self, cur: LeanPair) {
        if !self.seeded {
            let p0 = self.p0.expect("a candidate implies a length-passer => p0 set");
            self.ref_start = p0.left_tir_start;
            self.ref_end = p0.right_transformed_end;
            self.seeded = true;
            if !p0.skip {
                // p0 == this first live candidate (global min); consume it as the seed.
                self.best = Some(cur);
                return;
            }
            // p0 skipped: phantom best (None); fall through to process `cur`.
            self.best = None;
        }
        if self.ref_start <= cur.right_transformed_end && self.ref_end >= cur.left_tir_start {
            self.ref_end = self.ref_end.max(cur.right_transformed_end);
            self.best = match self.best {
                None => Some(cur), // phantom (skipped p0) always loses
                Some(b) => {
                    if double_smaller(b.similarity, cur.similarity) {
                        Some(cur) // cur strictly better; drop the old best
                    } else {
                        Some(b) // cur loses; drop it
                    }
                }
            };
        } else {
            // cluster closed: emit its survivor, then start a new cluster at cur.
            if let Some(b) = self.best {
                self.out.push(make_element(&b, self.e, self.contigs));
            }
            self.ref_start = cur.left_tir_start;
            self.ref_end = cur.right_transformed_end;
            self.best = Some(cur);
        }
    }

    /// Drain the remaining window and emit the final cluster's survivor.
    fn finish(mut self) -> Vec<Element> {
        self.drain_to(u32::MAX);
        if let Some(b) = self.best {
            self.out.push(make_element(&b, self.e, self.contigs));
        }
        self.out
    }
}

/// Stages 2-4 for a single seed: xdrop extension -> build_pair (length gate) ->
/// TSD search (mutates arm boundaries) -> similarity gate. Returns the resulting
/// TirPair (possibly skip=true; stage 5 still needs it in the array) or None if
/// the seed fails the length gate in build_pair. Reads only `e`/`scores`/`dist`/`p`
/// and the thread_local front buffers, so it is safe to call concurrently from
/// rayon workers.
fn process_seed(
    s: &Seed,
    e: &Encoded,
    scores: &ArbitraryScores,
    dist: &ArbitraryDistances,
    p: &Params,
) -> Option<TirPair> {
    let (s1, e1, s2, e2) = e.contig_bounds(s.contignumber);
    // Widen the u32 seed coords back to u64 for the extension/build arithmetic.
    let (pos1, pos2, slen) = (s.pos1 as u64, s.pos2 as u64, s.len as u64);
    let alilen = p.max_tir_len - slen;
    let (xl, xr) = extend_seed(
        &e.twobit, pos1, pos2, slen, s1, e1, s2, e2, alilen, scores, dist,
        p.xdrop_belowscore,
    );
    let mut pair = build_pair(
        pos1, pos2, slen, s.contignumber, xl.ivalue, xl.jvalue, xr.ivalue, xr.jvalue,
        e.total_logical, p.min_tir_len, p.max_tir_len,
    )?;
    let seq_start = e.fwd_seqstart[s.contignumber as usize];
    let seq_len = e.fwd_seqlen[s.contignumber as usize];
    search_for_tsds(
        &mut pair, &e.enc, seq_start, seq_len, p.vicinity, p.min_tsd_len, p.max_tsd_len,
    );
    if !pair.skip
        && (pair.left_tir_end <= pair.left_tir_start
            || pair.right_tir_end <= pair.right_tir_start)
    {
        pair.skip = true;
    }
    if !pair.skip {
        compute_similarity(&mut pair, &e.twobit, p.similar, p.sim_mult);
    }
    Some(pair)
}

pub fn run(contigs: &[(String, Vec<u8>)], p: &Params) -> Vec<Element> {
    let mut e = encode(contigs);
    let nsuf = e.num_suffixes();

    // Stage 1: SA/LCP -> maxpairs -> store_seed. ALL the SA scaffolding is dead
    // once `seeds` is built — the per-seed loop only reads e.enc/e.twobit + contig
    // metadata — so free each piece as soon as its last reader finishes: the i32
    // SA+LCP (inner scope), e.sa_input (after SA build), suftab/lcptab (after
    // maxpairs), e.seqnum_of (after store_seed). Keeps the loop's peak off them.
    let mut seeds: Vec<Seed> = Vec::new();
    {
        // tirvish_rs targets pre-fragmented input (genomeSplitter chunks, ~5 Mbp),
        // so the mirrored text is far under the 32-bit suffix-array limit. We use
        // the i32 libsais path and store the SA/LCP as u32 (not u64) — half the
        // floor, and the values provably fit since the text length is < 2^31.
        // Guard loudly rather than corrupt if someone runs a too-large sequence.
        assert!(
            e.sa_input.len() < i32::MAX as usize,
            "tirvish_rs: input too large ({} bp mirrored, limit {}). It expects \
             pre-fragmented genomes (~5 Mbp chunks); split larger sequences first.",
            e.sa_input.len(), i32::MAX
        );
        let (suftab, lcptab) = {
            let (sa, lcp) = sa_lcp(&e.sa_input, e.k);
            let suftab: Vec<u32> = sa[..nsuf].iter().map(|&x| x as u32).collect();
            drop(sa); // free the i32 SA before allocating lcptab (avoid a 4-array peak)
            let lcptab: Vec<u32> = lcp[..nsuf].iter().map(|&x| x as u32).collect();
            (suftab, lcptab)
        }; // i32 lcp freed
        e.sa_input = Vec::new(); // SA built; libsais text no longer needed (~T*4 B)
        enumerate_maxpairs(&suftab, &lcptab, p.seed, ALPHA, &e.enc, |len, p1, p2| {
            store_seed(
                &mut seeds, len, p1, p2, e.midpos, e.total_logical, e.num_contigs,
                &e.seqnum_of, p.min_tir_dist, p.max_tir_dist, p.max_tir_len,
            );
        });
    } // suftab/lcptab freed
    e.seqnum_of = Vec::new(); // only store_seed read it (~T*4 B)

    let scores = ArbitraryScores {
        mat: p.xdrop_mat, mis: p.xdrop_mis, ins: p.xdrop_ins, del: p.xdrop_del,
    };
    let dist = calc_distances(&scores);

    // Stages 2-4 are embarrassingly parallel: each seed maps independently to an
    // Option<TirPair> (None = length-gate fail), touching only immutable shared
    // state (e/scores/dist/p) + the thread_local front buffers.
    //
    // MEMORY: retain ONLY the live (non-skipped) length-passers, not all of them.
    // gt keeps every length-passer in first_pairs for the overlap sort, but in
    // gt_tir_remove_overlaps a skipped pair at index >= 1 is `continue`d before it
    // can touch refrng/maxsim — inert and droppable. The one skipped pair that can
    // affect output is the global-first (index-0 seed); keep it separately as `p0`.
    // The seed index (via .enumerate()) reproduces the validated stable-sort tie
    // -break. Each kept candidate is projected to a compact LeanPair.
    // Stages 2-5, STREAMED: process seeds in pos1 order via waves of parallel blocks,
    // window-sort the candidates (bounded disorder), and run the overlap sweep
    // incrementally so the full candidate set never resides. Faithful to the old
    // whole-array path: each candidate carries its ORIGINAL (SA-order) index as the
    // sort tie-break, and OverlapSweep reproduces remove_overlaps_pruned exactly.
    //
    // Disorder bound: left_tir_start in [pos1 - (maxtirlen+vicinity), pos1 +
    // (vicinity+maxtsd)] -- xdrop-left <= alilen <= maxtirlen, TSD shift <= vicinity
    // -- so a candidate is safe to emit once the watermark (min pos1 of unprocessed
    // seeds) passes its left_tir_start by that margin.
    let disorder = (p.max_tir_len + p.vicinity + p.max_tsd_len) as u32 + 8;

    // pos1 order over the seeds, carrying the SA-order index (the stage-5 tie-break).
    let mut order: Vec<u32> = (0..seeds.len() as u32).collect();
    order.sort_by_key(|&i| seeds[i as usize].pos1);
    let blocks: Vec<&[u32]> = order.chunks(SEED_PAR_BLOCK).collect();

    let mut sweep = OverlapSweep::new(&e, contigs, disorder);
    let wave = rayon::current_num_threads().max(1);
    let mut bi = 0;
    while bi < blocks.len() {
        let hi = (bi + wave).min(blocks.len());
        // Process this wave's contiguous blocks in parallel; rayon `collect`
        // preserves block (pos1) order. Each block is serial internally.
        let results: Vec<(Vec<LeanPair>, Option<LeanPair>)> = blocks[bi..hi]
            .par_iter()
            .map(|&blk| {
                let mut live: Vec<LeanPair> = Vec::new();
                let mut p0: Option<LeanPair> = None;
                for &i in blk.iter() {
                    if let Some(pair) = process_seed(&seeds[i as usize], &e, &scores, &dist, p) {
                        let lean = LeanPair::from_pair(i, &pair); // seed_idx = SA-order index
                        let beats = match &p0 {
                            None => true,
                            Some(pp) => compare_tirs(&lean, pp).then(i.cmp(&pp.seed_idx)).is_lt(),
                        };
                        if beats {
                            p0 = Some(lean);
                        }
                        if !pair.skip {
                            live.push(lean);
                        }
                    }
                }
                // Sort each block's candidates HERE (in the parallel worker) so the
                // serial drain only merges already-sorted runs (adaptive sort_by ->
                // ~linear) instead of full-sorting the window each time. Keeps the
                // sort work on the parallel phase, restoring CPU utilization.
                live.sort_by(|a, b| compare_tirs(a, b).then(a.seed_idx.cmp(&b.seed_idx)));
                (live, p0)
            })
            .collect();
        // Reassemble in block order; advance the window watermark to the next block's
        // smallest pos1 (every earlier seed is now processed).
        for (k, (live, blk_p0)) in results.into_iter().enumerate() {
            if let Some(pp) = blk_p0 {
                sweep.observe_p0(pp);
            }
            for c in live {
                sweep.push(c);
            }
            let next = bi + k + 1;
            let watermark = if next < blocks.len() {
                seeds[blocks[next][0] as usize].pos1
            } else {
                u32::MAX
            };
            sweep.drain_to(watermark);
        }
        bi = hi;
    }
    sweep.finish()
}
