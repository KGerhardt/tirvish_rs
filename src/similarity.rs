//! Stage 4 — TIR similarity.
//!
//! Verbatim port of genometools `greedyunitedist` (greedyedist.c — Myers O(nd)
//! greedy UNIT edit distance over a flat "front" array) + the similarity gate in
//! gt_tir_searchforTIRs (tir_stream.c:607-622) and gt_double_smaller_double /
//! gt_double_compare / gt_double_relative_equal (mathsupport.c).
//!
//! ulen/vlen are end-start (NOT +1) of the POST-TSD left arm (forward) and the
//! UNCHANGED mirror right arm. useq = forward enc slice, vseq = mirror enc slice
//! (revcomp already stored). sim = 100*(1 - edist/max(ulen,vlen)); the pair is
//! skipped iff sim is (strictly, epsilon-aware) below the threshold.

use crate::encode::ALPHA;
use crate::tsd::TirPair;

const ABS_ERR: f64 = 1.0e-100; // GT_DBL_MAX_ABS_ERROR
const REL_ERR: f64 = 1.0e-8; // GT_DBL_MAX_REL_ERROR

fn relative_equal(d1: f64, d2: f64) -> bool {
    if (d1 - d2).abs() < ABS_ERR {
        return true;
    }
    let relerr = if d2.abs() > d1.abs() {
        ((d1 - d2) / d2).abs()
    } else {
        ((d1 - d2) / d1).abs()
    };
    relerr <= REL_ERR
}

fn double_compare(d1: f64, d2: f64) -> i32 {
    if relative_equal(d1, d2) {
        0
    } else if d1 > d2 {
        1
    } else {
        -1
    }
}

/// gt_double_smaller_double.
pub fn double_smaller(d1: f64, d2: f64) -> bool {
    double_compare(d1, d2) < 0
}

/// gt_seqabstract_lcp(forward=true, ...): matching run, breaks at a special.
fn fwd_lcp(useq: &[u32], vseq: &[u32], u: usize, v: usize) -> usize {
    let maxlen = (useq.len() - u).min(vseq.len() - v);
    let mut l = 0;
    while l < maxlen {
        let uc = useq[u + l];
        if uc >= ALPHA {
            break;
        }
        let vc = vseq[v + l];
        if vc >= ALPHA {
            break;
        }
        if uc != vc {
            break;
        }
        l += 1;
    }
    l
}

// gt: frontspecparms — (left, width) for distance p, parameter r.
fn frontspecparms(ulen: i64, vlen: i64, p: i64, r: i64) -> (i64, i64) {
    if r <= 0 {
        (-p, p + p + 1)
    } else {
        let left = (-ulen).max(-p);
        let width = vlen.min(p) - left + 1;
        (left, width)
    }
}

// gt: accessfront — front row value at diagonal k, or integermin if out of band.
#[inline]
fn access(fs: &[i64], offset: i64, left: i64, width: i64, k: i64, imin: i64) -> i64 {
    if left <= k && k < left + width {
        fs[(offset - left + k) as usize]
    } else {
        imin
    }
}

/// gt: greedyunitedist — exact unit (Levenshtein) edit distance of useq vs vseq.
pub fn greedy_unit_edist(useq: &[u32], vseq: &[u32]) -> u64 {
    let ulen = useq.len() as i64;
    let vlen = vseq.len() as i64;
    let imin = -(ulen.max(vlen));
    let mut fs: Vec<i64> = vec![0];

    // firstfrontforward: front[0] at index 0
    fs[0] = if ulen == 0 || vlen == 0 {
        0
    } else {
        fwd_lcp(useq, vseq, 0, 0) as i64
    };
    if ulen == vlen && fs[0] == vlen {
        return 0;
    }

    let (mut prev_offset, mut prev_left, mut prev_width) = (0i64, 0i64, 1i64);
    let mut kval = 1i64;
    let mut r = 1 - ulen.min(vlen);
    loop {
        let offset = prev_offset + prev_width;
        let (left, width) = frontspecparms(ulen, vlen, kval, r);
        let need = (offset + width) as usize;
        if fs.len() < need {
            fs.resize(need, 0);
        }
        // evalfrontforward
        let mut k = left;
        while k < left + width {
            let stored = if r <= 0 || k <= -r || k >= r {
                // evalentryforward: max of three diagonal moves, then LCP extend
                let mut t = access(&fs, prev_offset, prev_left, prev_width, k, imin) + 1;
                let below = access(&fs, prev_offset, prev_left, prev_width, k - 1, imin);
                if t < below {
                    t = below;
                }
                let above = access(&fs, prev_offset, prev_left, prev_width, k + 1, imin) + 1;
                if t < above {
                    t = above;
                }
                if t < 0 || t + k < 0 {
                    imin
                } else {
                    let mut tt = t;
                    if ulen != 0 && vlen != 0 && tt < ulen && tt + k < vlen {
                        tt += fwd_lcp(useq, vseq, tt as usize, (tt + k) as usize) as i64;
                    }
                    if tt > ulen || tt + k > vlen {
                        imin
                    } else {
                        tt
                    }
                }
            } else {
                imin
            };
            fs[(offset - left + k) as usize] = stored;
            k += 1;
        }
        // alignment finished when the end diagonal reaches ulen
        if access(&fs, offset, left, width, vlen - ulen, imin) == ulen {
            return kval as u64;
        }
        prev_offset = offset;
        prev_left = left;
        prev_width = width;
        kval += 1;
        r += 1;
    }
}

/// gt: the similarity block. Computes ulen/vlen/edist/similarity on `pair`,
/// sets skip if below threshold. Returns (ulen, vlen, edist) for validation.
pub fn compute_similarity(pair: &mut TirPair, enc: &[u32], threshold: f64) -> (u64, u64, u64) {
    let ulen = pair.left_tir_end - pair.left_tir_start;
    let vlen = pair.right_tir_end - pair.right_tir_start;
    let useq = &enc[pair.left_tir_start as usize..(pair.left_tir_start + ulen) as usize];
    let vseq = &enc[pair.right_tir_start as usize..(pair.right_tir_start + vlen) as usize];
    let edist = greedy_unit_edist(useq, vseq);
    pair.similarity = 100.0 * (1.0 - edist as f64 / ulen.max(vlen) as f64);
    if double_smaller(pair.similarity, threshold) {
        pair.skip = true;
    }
    (ulen, vlen, edist)
}
