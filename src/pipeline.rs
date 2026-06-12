//! Full de-novo TIRvish pipeline: stages 1-5 -> final element list (gold format).
//!
//! Stage 5 = sort by gt_tir_compare_TIRs (contignumber, left_tir_start,
//! right_transformed_start), gt_tir_remove_overlaps "best" (keep max-similarity
//! within each overlapping cluster), tir_compactboundaries (drop skipped), then
//! emit GFF coordinates (tir_stream.c:940-1001) in the gold-TSV shape.

use crate::encode::{encode, ALPHA};
use crate::maxpairs::enumerate_maxpairs;
use crate::params;
use crate::sa::sa_lcp;
use crate::seeds::{store_seed, Seed};
use crate::similarity::{compute_similarity, double_smaller};
use crate::tsd::{build_pair, search_for_tsds, TirPair};
use crate::encode::Encoded;
use crate::xdrop::{calc_distances, extend_seed, ArbitraryDistances, ArbitraryScores};
use rayon::prelude::*;
use std::cmp::Ordering;
use std::sync::atomic::{AtomicU64, Ordering as AtomicOrdering};

#[derive(Debug, Clone)]
pub struct Element {
    pub seqid: String,
    pub start: u64,
    pub stop: u64,
    pub tir1: u64,
    pub tir2: u64,
    pub tsd1: u64,
    pub tsd2: u64,
    pub sim: f64,
}

/// gt_tir_compare_TIRs.
fn compare_tirs(a: &TirPair, b: &TirPair) -> Ordering {
    a.contignumber
        .cmp(&b.contignumber)
        .then(a.left_tir_start.cmp(&b.left_tir_start))
        .then(a.right_transformed_start.cmp(&b.right_transformed_start))
}

/// gt_tir_remove_overlaps, "best" mode (keep max-similarity within each cluster).
fn remove_overlaps(pairs: &mut [TirPair]) {
    if pairs.is_empty() {
        return;
    }
    let mut maxsim_idx = 0usize; // maxsimboundaries
    let mut ref_start = pairs[0].left_tir_start;
    let mut ref_end = pairs[0].right_transformed_end;
    for i in 1..pairs.len() {
        if pairs[i].skip {
            continue;
        }
        // tirboundaries_overlap(refrng, pairs[i])
        if ref_start <= pairs[i].right_transformed_end && ref_end >= pairs[i].left_tir_start {
            ref_end = ref_end.max(pairs[i].right_transformed_end);
            if double_smaller(pairs[maxsim_idx].similarity, pairs[i].similarity) {
                pairs[maxsim_idx].skip = true;
                maxsim_idx = i;
            } else {
                pairs[i].skip = true;
            }
        } else {
            ref_start = pairs[i].left_tir_start;
            ref_end = pairs[i].right_transformed_end;
            maxsim_idx = i;
        }
    }
}

/// Stages 2-4 for a single seed: xdrop extension -> build_pair (length gate) ->
/// TSD search (mutates arm boundaries) -> similarity gate. Returns the resulting
/// TirPair (possibly skip=true; stage 5 still needs it in the array) or None if
/// the seed fails the length gate in build_pair. Pure w.r.t. shared state — only
/// reads `e`/`scores`/`dist` and the thread_local front buffers — so it is safe
/// to call concurrently from rayon workers. The atomics accumulate per-stage
/// CPU-time only when `timeit` is set.
#[allow(clippy::too_many_arguments)]
fn process_seed(
    s: &Seed,
    e: &Encoded,
    scores: &ArbitraryScores,
    dist: &ArbitraryDistances,
    timeit: bool,
    a_xdrop: &AtomicU64,
    a_tsd: &AtomicU64,
    a_sim: &AtomicU64,
) -> Option<TirPair> {
    let (s1, e1, s2, e2) = e.contig_bounds(s.contignumber);
    let alilen = params::MAX_TIR_LEN - s.len;
    let tx = if timeit { Some(std::time::Instant::now()) } else { None };
    let (xl, xr) = extend_seed(
        &e.enc, s.pos1, s.pos2, s.len, s1, e1, s2, e2, alilen, scores, dist,
        params::XDROP_BELOWSCORE,
    );
    if let Some(tx) = tx {
        a_xdrop.fetch_add(tx.elapsed().as_nanos() as u64, AtomicOrdering::Relaxed);
    }
    let mut pair = build_pair(
        s.pos1, s.pos2, s.len, s.contignumber, xl.ivalue, xl.jvalue, xr.ivalue, xr.jvalue,
        e.total_logical, params::MIN_TIR_LEN, params::MAX_TIR_LEN,
    )?;
    let seq_start = e.fwd_seqstart[s.contignumber as usize];
    let seq_len = e.fwd_seqlen[s.contignumber as usize];
    let tt = if timeit { Some(std::time::Instant::now()) } else { None };
    search_for_tsds(
        &mut pair, &e.enc, seq_start, seq_len, params::VICINITY,
        params::MIN_TSD_LEN, params::MAX_TSD_LEN,
    );
    if let Some(tt) = tt {
        a_tsd.fetch_add(tt.elapsed().as_nanos() as u64, AtomicOrdering::Relaxed);
    }
    if !pair.skip
        && (pair.left_tir_end <= pair.left_tir_start
            || pair.right_tir_end <= pair.right_tir_start)
    {
        pair.skip = true;
    }
    if !pair.skip {
        let ts = if timeit { Some(std::time::Instant::now()) } else { None };
        compute_similarity(&mut pair, &e.twobit, params::SIMILARITY_THRESHOLD);
        if let Some(ts) = ts {
            a_sim.fetch_add(ts.elapsed().as_nanos() as u64, AtomicOrdering::Relaxed);
        }
    }
    Some(pair)
}

pub fn run(contigs: &[(String, Vec<u8>)]) -> Vec<Element> {
    // env-gated per-stage timing (TIRVISH_RS_TIME), to compare vs gt's breakdown.
    let timeit = std::env::var("TIRVISH_RS_TIME").is_ok();
    let t_start = std::time::Instant::now();
    let mut d_xdrop = std::time::Duration::ZERO;
    let mut d_tsd = std::time::Duration::ZERO;
    let mut d_sim = std::time::Duration::ZERO;

    let e = encode(contigs);
    let nsuf = e.num_suffixes();
    let (sa, lcp) = sa_lcp(&e.sa_input, e.k);
    let suftab: Vec<u64> = sa[..nsuf].iter().map(|&x| x as u64).collect();
    let lcptab: Vec<u64> = lcp[..nsuf].iter().map(|&x| x as u64).collect();

    let mut seeds: Vec<Seed> = Vec::new();
    enumerate_maxpairs(&suftab, &lcptab, params::SEED, ALPHA, &e.enc, |len, p1, p2| {
        store_seed(
            &mut seeds, len, p1, p2, e.midpos, e.total_logical, e.num_contigs,
            &e.seqnum_of, params::MIN_TIR_DIST, params::MAX_TIR_DIST, params::MAX_TIR_LEN,
        );
    });

    let t_stage1 = t_start.elapsed();
    let scores = ArbitraryScores {
        mat: params::XDROP_MAT, mis: params::XDROP_MIS,
        ins: params::XDROP_INS, del: params::XDROP_DEL,
    };
    let dist = calc_distances(&scores);

    // first_pairs = ALL length-passing seeds (skip flags from stages 3+4).
    //
    // The per-seed body (stages 2-4: xdrop -> build_pair -> TSD -> similarity) is
    // embarrassingly parallel — each seed maps independently to an Option<TirPair>
    // and touches only immutable shared state (e.enc/e.twobit/scores/dist) plus
    // the thread_local xdrop/greedyedist front buffers. par_iter preserves input
    // order and stage 5 re-sorts, so the result is bit-exact for any worker count.
    //
    // Per-stage timing (TIRVISH_RS_TIME) accumulates CPU-time summed across all
    // rayon workers via relaxed atomics; the loop's wall time is t_after_loop -
    // t_stage1. With N workers the summed CPU-time can exceed the loop wall.
    let a_xdrop = AtomicU64::new(0);
    let a_tsd = AtomicU64::new(0);
    let a_sim = AtomicU64::new(0);
    let mut pairs: Vec<TirPair> = seeds
        .par_iter()
        .filter_map(|s| {
            process_seed(s, &e, &scores, &dist, timeit, &a_xdrop, &a_tsd, &a_sim)
        })
        .collect();
    if timeit {
        d_xdrop = std::time::Duration::from_nanos(a_xdrop.load(AtomicOrdering::Relaxed));
        d_tsd = std::time::Duration::from_nanos(a_tsd.load(AtomicOrdering::Relaxed));
        d_sim = std::time::Duration::from_nanos(a_sim.load(AtomicOrdering::Relaxed));
    }

    // stage 5
    let t_after_loop = t_start.elapsed();
    pairs.sort_by(compare_tirs);
    remove_overlaps(&mut pairs);

    let mut out = Vec::new();
    for pair in &pairs {
        if pair.skip {
            continue;
        }
        let seqstart = e.fwd_seqstart[pair.contignumber as usize];
        let name = contigs[pair.contignumber as usize].0.split(";;").next().unwrap().to_string();
        // gt emits the two TIR features sorted by GtRange (start, then end), so
        // parse_tirvish reads tir1 = the (start,end)-first arm. Normally the left
        // arm comes first, but post-TSD the transformed right arm can start at or
        // before the left arm; when they share a start, the smaller end wins.
        let left_len = pair.left_tir_end - pair.left_tir_start + 1;
        let right_len = pair.right_transformed_end - pair.right_transformed_start + 1;
        let left_key = (pair.left_tir_start, pair.left_tir_end);
        let right_key = (pair.right_transformed_start, pair.right_transformed_end);
        let (tir1, tir2) = if left_key <= right_key {
            (left_len, right_len)
        } else {
            (right_len, left_len)
        };
        out.push(Element {
            seqid: name,
            start: pair.left_tir_start - seqstart - pair.tsd_length + 1,
            stop: pair.right_transformed_end - seqstart + pair.tsd_length + 1,
            tir1,
            tir2,
            tsd1: pair.tsd_length,
            tsd2: pair.tsd_length,
            sim: pair.similarity,
        });
    }
    if timeit {
        let total = t_start.elapsed();
        // loop_wall is the wall time of the parallel per-seed loop; xdrop/tsd/sim
        // are CPU-seconds SUMMED across rayon workers (so their sum can exceed
        // loop_wall by ~Nthreads). nthreads from rayon's global pool.
        let loop_wall = t_after_loop.saturating_sub(t_stage1);
        eprintln!(
            "RSTIME stage1={:.1}s loop_wall={:.1}s [cpu: xdrop={:.1}s tsd={:.1}s sim={:.1}s] stage5={:.1}s total={:.1}s (nthreads={})",
            t_stage1.as_secs_f64(), loop_wall.as_secs_f64(), d_xdrop.as_secs_f64(),
            d_tsd.as_secs_f64(), d_sim.as_secs_f64(),
            (total - t_after_loop).as_secs_f64(), total.as_secs_f64(),
            rayon::current_num_threads()
        );
    }
    out
}
