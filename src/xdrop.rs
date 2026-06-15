//! Stage 2 — Xdrop arm extension.
//!
//! Verbatim port of genometools 1.6.5 `xdrop.c` (`gt_evalxdroparbitscoresextend`,
//! the Zhang–Schwartz–Miller greedy Xdrop / Myers O(nd) with X-drop pruning) plus
//! `gt_seqabstract_lcp` and the per-seed left/right extension setup from
//! `gt_tir_searchforTIRs` (tir_stream.c:514-600).
//!
//! Sequences are slices of the mirror `enc` (ACGT 0..3, special == ALPHA). The
//! mirror already stores reverse-complement bases, so `dir_is_complement` is moot;
//! `GT_ISSPECIAL` becomes `code >= ALPHA`.

use crate::twobit::TwoBit;

#[derive(Clone, Copy)]
pub struct ArbitraryScores {
    pub mat: i32,
    pub mis: i32,
    pub ins: i32,
    pub del: i32,
}

#[derive(Clone, Copy)]
pub struct ArbitraryDistances {
    pub mis: i32,
    pub ins: i32,
    pub del: i32,
    pub gcd: i32,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct XdropBest {
    pub ivalue: u64,
    pub jvalue: u64,
    pub score: i64,
    pub best_d: i64,
    pub best_k: i64,
}

fn gcd_uint(mut a: u32, mut b: u32) -> u32 {
    while b != 0 {
        let t = b;
        b = a % b;
        a = t;
    }
    a
}

/// gt: gt_calculatedistancesfromscores (xdrop.c:128-156).
pub fn calc_distances(s: &ArbitraryScores) -> ArbitraryDistances {
    let (mat, mis, ins, del);
    if (s.mat as u32) % 2 > 0 {
        mat = s.mat * 2;
        mis = s.mis * 2;
        ins = s.ins * 2;
        del = s.del * 2;
    } else {
        mat = s.mat;
        mis = s.mis;
        ins = s.ins;
        del = s.del;
    }
    let gcd = gcd_uint(
        gcd_uint((mat - mis) as u32, (mat / 2 - ins) as u32),
        (mat / 2 - del) as u32,
    ) as i32;
    ArbitraryDistances {
        mis: (mat - mis) / gcd,
        ins: (mat / 2 - ins) / gcd,
        del: (mat / 2 - del) / gcd,
        gcd,
    }
}

const REPLACEMENTBIT: u8 = 1;
const DELETIONBIT: u8 = 2;
const INSERTIONBIT: u8 = 4;

/// gt: gt_evalxdroparbitscoresextend (xdrop.c:224-431). `forward` is gt's
/// `rightextension`. Returns the best-scoring extension; ivalue/jvalue are how
/// far it reached in useq/vseq.
thread_local! {
    // (f_cur, f_prev, big_t): two rolling front rows + the per-generation best
    // scores. Reused across calls (capacity retained).
    static XDROP_SCRATCH: std::cell::RefCell<(Vec<i64>, Vec<i64>, Vec<i64>)> =
        std::cell::RefCell::new((Vec::new(), Vec::new(), Vec::new()));
}

#[allow(clippy::too_many_arguments)]
pub fn eval_xdrop(
    forward: bool,
    ulen_u: usize,
    vlen_u: usize,
    tb: &TwoBit,
    ubase: usize,
    vbase: usize,
    xdropbelowscore: i64,
    scores: &ArbitraryScores,
    dist: &ArbitraryDistances,
) -> XdropBest {
    let ulen = ulen_u as i64;
    let vlen = vlen_u as i64;
    debug_assert!(ulen != 0 && vlen != 0);
    // SWAR LCE over the 2-bit mirror genome — identical value to the old scalar
    // seq_lcp (tb.lce is unit-tested == naive), just 32 bases/op. Forward (right
    // extension) ascends from ubase/vbase; reverse (left extension) descends.
    let lce = |u_start: i64, v_start: i64| -> i64 {
        let maxlen = (ulen - u_start).min(vlen - v_start);
        if maxlen <= 0 {
            return 0;
        }
        let m = maxlen as usize;
        if forward {
            tb.lce(true, ubase + u_start as usize, vbase + v_start as usize, m) as i64
        } else {
            tb.lce(false, ubase - u_start as usize, vbase - v_start as usize, m) as i64
        }
    };
    let end_k = ulen - vlen;
    let integermax = ulen.max(vlen);
    let integermin = -integermax;
    let mat_half = (scores.mat / 2) as i64;
    let gcd = dist.gcd as i64;
    let dmis = dist.mis as i64;
    let dins = dist.ins as i64;
    let ddel = dist.del as i64;
    let dback = (xdropbelowscore + mat_half) / gcd + 1;
    let allowed_mininf = dmis.max(dins).max(ddel) - 1;
    let eval = |k_anti: i64, d: i64| k_anti * mat_half - d * gcd;

    XDROP_SCRATCH.with(|cell| {
    let mut guard = cell.borrow_mut();
    let (f_cur, f_prev, big_t) = &mut *guard;
    big_t.clear();
    // ROLLING FRONT: for tirvish's fixed scores ddel=dmis=dins=1, so every front
    // read is from generation currd-1 -> keep just two rows (f_cur=currd,
    // f_prev=currd-1) instead of gt's d^2+d+k triangular array. Diagonal k of
    // generation g lives at local index (g+k), so the d^2+d+k multiply becomes an
    // add and the hot rows stay tiny (L1/L2). gt's per-generation fills overwrite
    // every k in [-currd,currd], so no stale slot is ever read. fd (gt traceback
    // direction) is never read here, so it is gone.
    debug_assert!(ddel == 1 && dmis == 1 && dins == 1);
    let mut best = XdropBest::default();

    let mut currd: i64 = 0;
    let lbound;
    let ubound;

    // phase 0 -> generation 0 lives in f_prev (the "previous" row for currd=1).
    let idx0 = lce(0, 0);
    if idx0 >= ulen || idx0 >= vlen {
        lbound = 1;
        ubound = -1;
    } else {
        lbound = 0;
        ubound = 0;
    }
    let mut tmp_dir: u8 = 0;
    if f_prev.is_empty() {
        f_prev.push(integermin);
    }
    f_prev[0] = idx0; // gen 0, k=0 -> local 0
    let mut bigt_tmp = eval(idx0 + idx0, 0);
    best.score = bigt_tmp;
    best.ivalue = idx0 as u64;
    best.jvalue = idx0 as u64;
    best.best_d = currd;
    best.best_k = 0;
    big_t.push(bigt_tmp);

    let mut lbound = lbound;
    let mut ubound = ubound;
    let mut current_mininf_gen: i64 = 0;
    let mut always_mininf = true;

    while lbound <= ubound {
        currd += 1;
        // current row holds k in [-currd, currd] -> local [0, 2*currd].
        let need = (2 * currd + 1) as usize;
        if f_cur.len() < need {
            f_cur.resize(need, integermin);
        }
        let gc = |k: i64| (currd + k) as usize; // f_cur (gen currd)
        let gp = |k: i64| (currd - 1 + k) as usize; // f_prev (gen currd-1)
        let mut k = lbound - 1;
        while k <= ubound + 1 {
            let mut i = integermin;
            // case 1: DELETION-EDGE (reads gen currd-1)
            if lbound < k && -(currd - 1) <= k - 1 && k - 1 <= currd - 1 {
                i = f_prev[gp(k - 1)] + 1;
                tmp_dir = DELETIONBIT;
            }
            // case 2: REPLACEMENT-EDGE
            if lbound <= k && k <= ubound && -(currd - 1) <= k && k <= currd - 1 {
                let row = f_prev[gp(k)] + 1;
                if (tmp_dir & DELETIONBIT) == 0 || row > i {
                    i = row;
                    tmp_dir = REPLACEMENTBIT;
                }
            }
            // case 3: INSERTION-EDGE
            if k < ubound && -(currd - 1) <= k + 1 && k + 1 <= currd - 1 {
                let row = f_prev[gp(k + 1)];
                if (tmp_dir & (DELETIONBIT | REPLACEMENTBIT)) == 0 || row > i {
                    i = row;
                    tmp_dir = INSERTIONBIT;
                }
            }
            let tmp_row;
            if i < 0 {
                if tmp_dir == 0 {
                    always_mininf = false;
                }
                tmp_row = integermin;
            } else {
                let mut i = i;
                let mut j = i - k;
                let previousd = currd - dback;
                if previousd > 0
                    && (previousd as usize) < big_t.len()
                    && eval(i + j, currd) < big_t[previousd as usize] - xdropbelowscore
                {
                    tmp_row = integermin;
                } else if k <= -currd
                    || k >= currd
                    || (f_prev[gp(k)] < i && i <= ulen.min(vlen + k))
                {
                    if ulen > i && vlen > j {
                        let lcp = lce(i, j);
                        i += lcp;
                        j += lcp;
                    }
                    always_mininf = false;
                    tmp_row = i;
                    if eval(i + j, currd) > bigt_tmp {
                        bigt_tmp = eval(i + j, currd);
                        best.score = bigt_tmp;
                        best.ivalue = i as u64;
                        best.jvalue = j as u64;
                        best.best_d = currd;
                        best.best_k = k;
                    }
                } else {
                    always_mininf = false;
                    tmp_row = f_prev[gp(k)];
                }
            }
            f_cur[gc(k)] = tmp_row;
            k += 1;
        }
        if always_mininf {
            current_mininf_gen += 1;
            if current_mininf_gen > allowed_mininf {
                break;
            }
        } else {
            current_mininf_gen = 0;
            always_mininf = true;
        }
        big_t.push(bigt_tmp);
        // fill out-of-bounds fronts with integermin (gt does this too)
        let mut kk = -currd;
        while kk < lbound - 1 {
            f_cur[gc(kk)] = integermin;
            kk += 1;
        }
        let mut kk = ubound + 2;
        while kk <= currd {
            f_cur[gc(kk)] = integermin;
            kk += 1;
        }
        // alignment finished
        if -currd <= end_k && end_k <= currd && f_cur[gc(end_k)] == ulen {
            break;
        }
        // pruning lower bound
        let mut kk = lbound - 1;
        while kk <= ubound + 1 {
            if f_cur[gc(kk)] > integermin {
                lbound = kk;
                break;
            }
            kk += 1;
        }
        // pruning upper bound
        let mut kk = ubound + 1;
        while kk >= lbound - 1 {
            if f_cur[gc(kk)] > integermin {
                ubound = kk;
                break;
            }
            kk -= 1;
        }
        // handling boundaries lower bound
        let mut kk = 0;
        while kk >= lbound {
            if f_cur[gc(kk)] == vlen + kk {
                lbound = kk;
                break;
            }
            kk -= 1;
        }
        // handling boundaries upper bound
        let mut kk = 0;
        while kk <= ubound {
            if f_cur[gc(kk)] == ulen {
                ubound = kk;
                break;
            }
            kk += 1;
        }
        // current generation becomes previous for the next.
        std::mem::swap(f_cur, f_prev);
    }
    best
    })
}

/// Per-seed left (reverse) + right (forward) extension, mirroring
/// gt_tir_searchforTIRs (tir_stream.c:514-600). Returns (xdrop_left, xdrop_right).
#[allow(clippy::too_many_arguments)]
pub fn extend_seed(
    tb: &TwoBit,
    pos1: u64,
    pos2: u64,
    len: u64,
    seqstart1: u64,
    seqend1: u64,
    seqstart2: u64,
    seqend2: u64,
    alilen: u64,
    scores: &ArbitraryScores,
    dist: &ArbitraryDistances,
    belowscore: i64,
) -> (XdropBest, XdropBest) {
    let mut xleft = XdropBest::default();
    let mut xright = XdropBest::default();

    // left (reverse) xdrop. Both arms have length l; the LCE descends from the
    // base just inside the seed (pos1-1 / pos2-1).
    if alilen != 0 && pos1 > seqstart1 && pos2 > seqstart2 {
        let l = if alilen <= pos1 - seqstart1 && alilen <= pos2 - seqstart2 {
            alilen
        } else {
            (pos1 - seqstart1).min(pos2 - seqstart2)
        };
        xleft = eval_xdrop(
            false, l as usize, l as usize, tb, (pos1 - 1) as usize, (pos2 - 1) as usize,
            belowscore, scores, dist,
        );
    }

    // right (forward) xdrop. LCE ascends from the base just past the seed.
    if alilen != 0 && pos1 + len < seqend1 && pos2 + len < seqend2 {
        let l = if alilen <= seqend1 - (pos1 + len) && alilen <= seqend2 - (pos2 + len) {
            alilen
        } else {
            (seqend1 - (pos1 + len)).min(seqend2 - (pos2 + len))
        };
        xright = eval_xdrop(
            true, l as usize, l as usize, tb, (pos1 + len) as usize, (pos2 + len) as usize,
            belowscore, scores, dist,
        );
    }
    (xleft, xright)
}
