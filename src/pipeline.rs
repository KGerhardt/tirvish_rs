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
/// 64 B vs TirPair's 104 B; with ~3.67M live candidates on a shrimp chunk this is
/// the dominant memory term. `seed_idx` carries the stable-sort tie-break.
#[derive(Clone, Copy)]
struct LeanPair {
    seed_idx: u32,
    contignumber: u32,
    left_tir_start: u64,
    left_tir_end: u64,
    right_transformed_start: u64,
    right_transformed_end: u64,
    tsd_length: u64,
    similarity: f64,
    skip: bool,
}

impl LeanPair {
    fn from_pair(seed_idx: u32, p: &TirPair) -> Self {
        LeanPair {
            seed_idx,
            contignumber: p.contignumber,
            left_tir_start: p.left_tir_start,
            left_tir_end: p.left_tir_end,
            right_transformed_start: p.right_transformed_start,
            right_transformed_end: p.right_transformed_end,
            tsd_length: p.tsd_length,
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

/// Result of the parallel per-seed pass: the live (non-skipped) length-passers,
/// each tagged with its seed index, plus `p0` = the single global-first
/// length-passer by (compare_tirs, seed_index) — gt's `pairs[0]`, which may
/// itself be skipped. See [`remove_overlaps_pruned`] for why only `p0` matters.
#[derive(Default)]
struct Collected {
    live: Vec<LeanPair>,
    p0: Option<LeanPair>,
}

/// gt_tir_remove_overlaps, "best" mode (keep max-similarity per overlapping
/// cluster), over the PRUNED candidate set.
///
/// `pairs` holds only the live (non-skipped) length-passers, sorted by
/// (compare_tirs, seed_index). `p0` is gt's `pairs[0]` — the global-first
/// length-passer, possibly skipped — which gt uses to seed the reference range
/// and `maxsimboundaries` at index 0 *unconditionally* (tir_stream.c:250-253).
///
/// This is output-identical to running gt's algorithm over the full first_pairs
/// array: every skipped pair at index >= 1 is `continue`d before it can touch
/// `refrng`/`maxsim`, so it is inert and droppable; only the index-0 seed
/// survives pruning, carried here as `p0`. When `p0` is skipped its similarity
/// is < 80 <= every live pair's, so it can extend the initial `ref_end` but can
/// never win the max-similarity comparison — modelled as a `None` (phantom)
/// maxsim that the first overlapping live pair always replaces. When `p0` is
/// live it equals `pairs[0]` (global min, present in the live set), so it seeds
/// from index 0 and the loop starts at 1, exactly as gt does.
fn remove_overlaps_pruned(pairs: &mut [LeanPair], p0: Option<&LeanPair>) {
    let p0 = match p0 {
        Some(p) => p,
        None => return, // no length-passers at all
    };
    let mut ref_start = p0.left_tir_start;
    let mut ref_end = p0.right_transformed_end;
    // maxsim: index into `pairs`, or None for the phantom skipped-p0 seed.
    let (mut maxsim, start) = if p0.skip {
        (None, 0usize) // p0 not in `pairs`; process every live pair
    } else {
        (Some(0usize), 1usize) // p0 == pairs[0]; seed from it, process the rest
    };
    for i in start..pairs.len() {
        // tirboundaries_overlap(refrng, pairs[i]); all of `pairs` is live here.
        if ref_start <= pairs[i].right_transformed_end && ref_end >= pairs[i].left_tir_start {
            ref_end = ref_end.max(pairs[i].right_transformed_end);
            match maxsim {
                None => maxsim = Some(i), // phantom (skipped p0) always loses
                Some(m) => {
                    if double_smaller(pairs[m].similarity, pairs[i].similarity) {
                        pairs[m].skip = true;
                        maxsim = Some(i);
                    } else {
                        pairs[i].skip = true;
                    }
                }
            }
        } else {
            ref_start = pairs[i].left_tir_start;
            ref_end = pairs[i].right_transformed_end;
            maxsim = Some(i);
        }
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
    let alilen = p.max_tir_len - s.len;
    let (xl, xr) = extend_seed(
        &e.twobit, s.pos1, s.pos2, s.len, s1, e1, s2, e2, alilen, scores, dist,
        p.xdrop_belowscore,
    );
    let mut pair = build_pair(
        s.pos1, s.pos2, s.len, s.contignumber, xl.ivalue, xl.jvalue, xr.ivalue, xr.jvalue,
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
            let lcptab: Vec<u32> = lcp[..nsuf].iter().map(|&x| x as u32).collect();
            (suftab, lcptab)
        }; // i32 sa/lcp freed
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
    let collected = seeds
        .par_iter()
        .with_min_len(SEED_PAR_BLOCK)
        .enumerate()
        .fold(Collected::default, |mut acc, (i, s)| {
            if let Some(pair) = process_seed(s, &e, &scores, &dist, p) {
                let idx = i as u32;
                let lean = LeanPair::from_pair(idx, &pair);
                // p0 = argmin over ALL length-passers (skipped or not) by
                // (compare_tirs, seed_index) = gt's pairs[0].
                let beats = match &acc.p0 {
                    None => true,
                    Some(pp) => compare_tirs(&lean, pp).then(idx.cmp(&pp.seed_idx)).is_lt(),
                };
                if beats {
                    acc.p0 = Some(lean);
                }
                if !pair.skip {
                    acc.live.push(lean);
                }
            }
            acc
        })
        .reduce(Collected::default, |mut a, mut b| {
            a.live.append(&mut b.live);
            a.p0 = match (a.p0, b.p0) {
                (None, x) | (x, None) => x,
                (Some(x), Some(y)) => {
                    if compare_tirs(&x, &y).then(x.seed_idx.cmp(&y.seed_idx)).is_lt() {
                        Some(x)
                    } else {
                        Some(y)
                    }
                }
            };
            a
        });

    // stage 5: sort the live survivors (compare_tirs, then seed index to match the
    // original stable full-array sort), then remove_overlaps seeded from p0.
    let mut live = collected.live;
    live.sort_by(|a, b| compare_tirs(a, b).then(a.seed_idx.cmp(&b.seed_idx)));
    remove_overlaps_pruned(&mut live, collected.p0.as_ref());

    let mut out = Vec::new();
    for pair in &live {
        if pair.skip {
            continue;
        }
        let seqstart = e.fwd_seqstart[pair.contignumber as usize];
        // Emit the sequence id verbatim — never mutate identifiers the input carries.
        let name = contigs[pair.contignumber as usize].0.clone();
        // The six (start, stop) pairs gt emits per element, 1-based and local to
        // the contig. The full element is anchored on the left-arm start and the
        // right-arm end (± one TSD); the body (no-TSD element) sits one TSD inside
        // each end; the two TSDs flank the body.
        let tsd = pair.tsd_length;
        let full_start = pair.left_tir_start - seqstart - tsd + 1;
        let full_stop = pair.right_transformed_end - seqstart + tsd + 1;
        let body_start = full_start + tsd;
        let body_stop = full_stop - tsd;
        let tsd1_start = full_start;
        let tsd1_stop = body_start - 1;
        let tsd2_start = body_stop + 1;
        let tsd2_stop = full_stop;
        // The two TIR arms, each in its own coordinates. gt emits the two TIR
        // features sorted by GtRange (start, then end), so TIR1 = the (start,end)-
        // first arm. Post-TSD the transformed right arm can start at or before the
        // left arm (the arms cross over when long), so sort rather than assume
        // left-then-right; when starts tie the smaller end wins.
        let left_arm = (
            pair.left_tir_start - seqstart + 1,
            pair.left_tir_end - seqstart + 1,
        );
        let right_arm = (
            pair.right_transformed_start - seqstart + 1,
            pair.right_transformed_end - seqstart + 1,
        );
        let (tir1, tir2) = if left_arm <= right_arm {
            (left_arm, right_arm)
        } else {
            (right_arm, left_arm)
        };
        out.push(Element {
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
        });
    }
    out
}
