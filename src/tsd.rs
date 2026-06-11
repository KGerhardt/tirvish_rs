//! Stage 3 — TSD search (target site duplication).
//!
//! Verbatim port of gt_tir_search_for_TSDs + gt_tir_store_TSDs + gt_tir_find_best_TSD
//! (tir_stream.c) and the matcher semantics of gt_sarrquerysubstringmatch /
//! gt_querysubstringmatch (esa-mmsearch.c): all LEFT-MAXIMAL exact matches of length
//! >= min_TSD between the left flank (db) and right flank (query), enumerated in
//! (query-offset ascending, then db-suffix SA order). The SA order matters: find_best
//! keeps the FIRST minimal-cost TSD (strict `<`).
//!
//! TSD search runs in FORWARD coordinates and MUTATES left_tir_start and
//! right_transformed_end BEFORE the similarity gate (so arms can shrink below the
//! seed). All flank bases are read from the forward `enc` (specials == ALPHA).

use crate::encode::ALPHA;
use crate::params::reversepos;
use std::cmp::Ordering;

#[derive(Debug, Clone)]
pub struct TirPair {
    pub pos1: u64,
    pub pos2: u64,
    pub seed_len: u64,
    pub contignumber: u32,
    pub left_tir_start: u64,
    pub left_tir_end: u64,
    pub right_tir_start: u64, // mirror coords
    pub right_tir_end: u64,   // mirror coords
    pub right_transformed_start: u64,
    pub right_transformed_end: u64,
    pub tsd_length: u64,
    pub skip: bool,
    pub similarity: f64,
}

/// One stored TSD (gt: Seed with pos1/offset/len reused as TSD fields).
struct Tsd {
    pos1: u64,   // absolute forward position of the TSD in the left flank
    offset: u64, // (right-flank absolute position) - pos1
    len: u64,
}

/// db-suffix lexicographic order == gt encseq SA order. End-of-string is the
/// separator, which sorts AFTER any base (special code), so a shorter suffix
/// (reaching the end first) sorts later.
fn cmp_suffix(db: &[u32], a: usize, b: usize) -> Ordering {
    let mut m = 0;
    loop {
        let ca = if a + m < db.len() { db[a + m] } else { u32::MAX };
        let cb = if b + m < db.len() { db[b + m] } else { u32::MAX };
        if ca != cb {
            return ca.cmp(&cb);
        }
        if a + m >= db.len() {
            return Ordering::Equal; // distinct positions => unreachable
        }
        m += 1;
    }
}

/// gt_querysubstringmatch: left-maximal exact matches (len >= minlen) of db vs
/// query, in (query-offset asc, db-suffix SA order). Stored as gt_tir_store_TSDs.
fn enumerate_tsds(db: &[u32], query: &[u32], minlen: usize, left_start: u64, right_start: u64) -> Vec<Tsd> {
    let mut out = Vec::new();
    let (dblen, qlen) = (db.len(), query.len());
    if dblen < minlen || qlen < minlen {
        return out;
    }
    for offset in 0..=(qlen - minlen) {
        // db positions whose minlen-prefix equals query[offset..offset+minlen]
        let mut matches: Vec<usize> = (0..=dblen - minlen)
            .filter(|&b| (0..minlen).all(|m| db[b + m] < ALPHA && db[b + m] == query[offset + m]))
            .collect();
        matches.sort_by(|&a, &b| cmp_suffix(db, a, b)); // SA order within the interval
        for dbstart in matches {
            // left-maximal (gt_mmsearch_isleftmaximal)
            let leftmax = dbstart == 0
                || offset == 0
                || db[dbstart - 1] >= ALPHA
                || db[dbstart - 1] != query[offset - 1];
            if !leftmax {
                continue;
            }
            // extend right (gt_mmsearch_extendright)
            let mut extend = 0usize;
            let (mut dbpos, mut qpos) = (dbstart + minlen, offset + minlen);
            while dbpos < dblen && qpos < qlen {
                if db[dbpos] >= ALPHA || db[dbpos] != query[qpos] {
                    break;
                }
                extend += 1;
                dbpos += 1;
                qpos += 1;
            }
            let matchlen = (minlen + extend) as u64;
            let pos1 = left_start + dbstart as u64;
            let right_abs = right_start + offset as u64;
            out.push(Tsd {
                pos1,
                offset: right_abs - pos1,
                len: matchlen,
            });
        }
    }
    out
}

/// gt_tir_find_best_TSD: pick the minimal-cost TSD (strict `<` => first wins),
/// mutate boundaries, set skip flags.
fn find_best_tsd(tsds: &[Tsd], pair: &mut TirPair, min_tsd: u64, max_tsd: u64) {
    let mut new_left_tir_start = pair.left_tir_start;
    let mut new_right_tir_end = pair.right_transformed_end;
    let mut best_cost = u64::MAX;
    let mut optimal_tsd_length = 0u64;

    for tsd in tsds {
        if tsd.len < min_tsd {
            continue;
        }
        let tsd_length = tsd.len;
        if tsd_length < max_tsd {
            let new_cost_left = if tsd.pos1 + tsd_length - 1 < pair.left_tir_start {
                pair.left_tir_start - (tsd.pos1 + tsd_length - 1)
            } else {
                (tsd.pos1 + tsd_length - 1) - pair.left_tir_start
            };
            let new_cost_right = if pair.right_transformed_end < tsd.pos1 + tsd.offset {
                (tsd.pos1 + tsd.offset) - pair.right_transformed_end
            } else {
                pair.right_transformed_end - (tsd.pos1 + tsd.offset)
            };
            let new_cost = new_cost_left + new_cost_right;
            if new_cost < best_cost {
                best_cost = new_cost;
                new_left_tir_start = tsd.pos1 + tsd_length;
                new_right_tir_end = tsd.pos1 + tsd.offset - 1;
                optimal_tsd_length = tsd_length;
            }
        }
    }

    if !tsds.is_empty() {
        pair.left_tir_start = new_left_tir_start;
        pair.right_transformed_end = new_right_tir_end;
        pair.tsd_length = optimal_tsd_length;
    } else {
        pair.skip = true;
    }
    if pair.right_transformed_end <= pair.right_transformed_start {
        pair.skip = true;
    }
    if pair.left_tir_end <= pair.left_tir_start {
        pair.skip = true;
    }
    if pair.tsd_length == 0 {
        pair.skip = true;
    }
}

/// gt_tir_search_for_TSDs: vicinity setup + enumerate + find_best.
#[allow(clippy::too_many_arguments)]
pub fn search_for_tsds(
    pair: &mut TirPair,
    enc: &[u32],
    seq_start_pos1: u64,
    seq_length: u64,
    vicinity: u64,
    min_tsd: u64,
    max_tsd: u64,
) {
    let seq_end_pos2 = seq_start_pos1 + seq_length - 1;

    let start_left = if pair.left_tir_start - seq_start_pos1 < vicinity {
        seq_start_pos1
    } else {
        pair.left_tir_start - vicinity
    };
    let end_left = if pair.left_tir_start + vicinity > pair.left_tir_end {
        pair.left_tir_end
    } else {
        pair.left_tir_start + vicinity
    };

    let start_right = if pair.right_transformed_start > pair.right_transformed_end - vicinity {
        pair.right_transformed_start
    } else {
        pair.right_transformed_end - vicinity
    };
    let end_right = if pair.right_transformed_end + vicinity > seq_end_pos2 {
        seq_end_pos2
    } else {
        pair.right_transformed_end + vicinity
    };

    if min_tsd > 1 {
        let db = &enc[start_left as usize..=end_left as usize];
        let query = &enc[start_right as usize..=end_right as usize];
        let tsds = enumerate_tsds(db, query, min_tsd as usize, start_left, start_right);
        find_best_tsd(&tsds, pair, min_tsd, max_tsd);
    }
}

/// Build the TIR pair from the seed + xdrop, apply gt's length check
/// (tir_stream.c:578-585). Returns None if the seed is dropped by the length
/// check (gt `continue`). `t` = total_logical (for reversepos).
#[allow(clippy::too_many_arguments)]
pub fn build_pair(
    pos1: u64,
    pos2: u64,
    seed_len: u64,
    contignumber: u32,
    li: u64, // xdrop_left.ivalue
    lj: u64, // xdrop_left.jvalue
    ri: u64, // xdrop_right.ivalue
    rj: u64, // xdrop_right.jvalue
    t: u64,
    min_tir_len: u64,
    max_tir_len: u64,
) -> Option<TirPair> {
    let tir_len = (pos1 + seed_len - 1 + ri) - (pos1 - lj + 1);
    if tir_len < min_tir_len || tir_len > max_tir_len {
        return None;
    }
    let left_tir_start = pos1 - li;
    let left_tir_end = pos1 + seed_len - 1 + ri;
    let right_tir_start = pos2 - lj;
    let right_tir_end = pos2 + seed_len - 1 + rj;
    let right_transformed_start = reversepos(t, right_tir_end);
    let right_transformed_end = reversepos(t, right_tir_start);
    Some(TirPair {
        pos1,
        pos2,
        seed_len,
        contignumber,
        left_tir_start,
        left_tir_end,
        right_tir_start,
        right_tir_end,
        right_transformed_start,
        right_transformed_end,
        tsd_length: 0,
        skip: false,
        similarity: 0.0,
    })
}
