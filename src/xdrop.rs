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

use crate::encode::ALPHA;

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

#[inline]
fn get_char(forward: bool, seq: &[u32], pos: usize) -> u32 {
    if forward {
        seq[pos]
    } else {
        seq[seq.len() - 1 - pos]
    }
}

/// gt: gt_seqabstract_lcp (seqabstract.c:205). Breaks at a special (>= ALPHA).
fn seq_lcp(forward: bool, useq: &[u32], vseq: &[u32], u_start: usize, v_start: usize) -> usize {
    let maxlen = (useq.len() - u_start).min(vseq.len() - v_start);
    let mut lcp = 0;
    while lcp < maxlen {
        let u_cc = get_char(forward, useq, u_start + lcp);
        if u_cc >= ALPHA {
            break;
        }
        let v_cc = get_char(forward, vseq, v_start + lcp);
        if v_cc >= ALPHA {
            break;
        }
        if u_cc != v_cc {
            break;
        }
        lcp += 1;
    }
    lcp
}

#[inline]
fn frontidx(d: i64, k: i64) -> usize {
    (d * d + d + k) as usize
}

/// gt: gt_evalxdroparbitscoresextend (xdrop.c:224-431). `forward` is gt's
/// `rightextension`. Returns the best-scoring extension; ivalue/jvalue are how
/// far it reached in useq/vseq.
thread_local! {
    static XDROP_SCRATCH: std::cell::RefCell<(Vec<i64>, Vec<u8>, Vec<i64>)> =
        std::cell::RefCell::new((Vec::new(), Vec::new(), Vec::new()));
}

pub fn eval_xdrop(
    forward: bool,
    useq: &[u32],
    vseq: &[u32],
    xdropbelowscore: i64,
    scores: &ArbitraryScores,
    dist: &ArbitraryDistances,
) -> XdropBest {
    let ulen = useq.len() as i64;
    let vlen = vseq.len() as i64;
    debug_assert!(ulen != 0 && vlen != 0);
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
    let (fr, fd, big_t) = &mut *guard;
    big_t.clear();
    // fronts indexed by frontidx(d,k); reused across calls (set-before-read, so
    // stale data is always overwritten before it could be read).
    let mut set_front = |d: i64, k: i64, row: i64, dir: u8, fr: &mut Vec<i64>, fd: &mut Vec<u8>| {
        let idx = frontidx(d, k);
        if idx >= fr.len() {
            fr.resize(idx + 1, integermin);
            fd.resize(idx + 1, 0);
        }
        fr[idx] = row;
        fd[idx] = dir;
    };
    let mut best = XdropBest::default();

    let mut currd: i64 = 0;
    let lbound;
    let ubound;

    // phase 0
    let idx0 = seq_lcp(forward, useq, vseq, 0, 0) as i64;
    if idx0 >= ulen || idx0 >= vlen {
        lbound = 1;
        ubound = -1;
    } else {
        lbound = 0;
        ubound = 0;
    }
    let mut tmp_dir: u8 = 0;
    set_front(0, 0, idx0, tmp_dir, fr, fd);
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
        let mut k = lbound - 1;
        while k <= ubound + 1 {
            let mut i = integermin;
            // case 1: DELETION-EDGE
            if lbound < k
                && currd - ddel >= 0
                && -(currd - ddel) <= k - 1
                && k - 1 <= currd - ddel
            {
                i = fr[frontidx(currd - ddel, k - 1)] + 1;
                tmp_dir = DELETIONBIT;
            }
            // case 2: REPLACEMENT-EDGE
            if lbound <= k
                && k <= ubound
                && currd - dmis >= 0
                && -(currd - dmis) <= k
                && k <= currd - dmis
            {
                let row = fr[frontidx(currd - dmis, k)] + 1;
                if (tmp_dir & DELETIONBIT) == 0 || row > i {
                    i = row;
                    tmp_dir = REPLACEMENTBIT;
                }
            }
            // case 3: INSERTION-EDGE
            if k < ubound
                && currd - dins >= 0
                && -(currd - dins) <= k + 1
                && k + 1 <= currd - dins
            {
                let row = fr[frontidx(currd - dins, k + 1)];
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
                    || (fr[frontidx(currd - 1, k)] < i && i <= ulen.min(vlen + k))
                {
                    if ulen > i && vlen > j {
                        let lcp = seq_lcp(forward, useq, vseq, i as usize, j as usize) as i64;
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
                    tmp_row = fr[frontidx(currd - 1, k)];
                }
            }
            set_front(currd, k, tmp_row, tmp_dir, fr, fd);
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
            set_front(currd, kk, integermin, tmp_dir, fr, fd);
            kk += 1;
        }
        let mut kk = ubound + 2;
        while kk <= currd {
            set_front(currd, kk, integermin, tmp_dir, fr, fd);
            kk += 1;
        }
        // alignment finished
        if -currd <= end_k && end_k <= currd && fr[frontidx(currd, end_k)] == ulen {
            break;
        }
        // pruning lower bound
        let mut kk = lbound - 1;
        while kk <= ubound + 1 {
            if fr[frontidx(currd, kk)] > integermin {
                lbound = kk;
                break;
            }
            kk += 1;
        }
        // pruning upper bound
        let mut kk = ubound + 1;
        while kk >= lbound - 1 {
            if fr[frontidx(currd, kk)] > integermin {
                ubound = kk;
                break;
            }
            kk -= 1;
        }
        // handling boundaries lower bound
        let mut kk = 0;
        while kk >= lbound {
            if fr[frontidx(currd, kk)] == vlen + kk {
                lbound = kk;
                break;
            }
            kk -= 1;
        }
        // handling boundaries upper bound
        let mut kk = 0;
        while kk <= ubound {
            if fr[frontidx(currd, kk)] == ulen {
                ubound = kk;
                break;
            }
            kk += 1;
        }
    }
    best
    })
}

/// Per-seed left (reverse) + right (forward) extension, mirroring
/// gt_tir_searchforTIRs (tir_stream.c:514-600). Returns (xdrop_left, xdrop_right).
#[allow(clippy::too_many_arguments)]
pub fn extend_seed(
    enc: &[u32],
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

    // left (reverse) xdrop
    if alilen != 0 && pos1 > seqstart1 && pos2 > seqstart2 {
        let l = if alilen <= pos1 - seqstart1 && alilen <= pos2 - seqstart2 {
            alilen
        } else {
            (pos1 - seqstart1).min(pos2 - seqstart2)
        };
        let useq = &enc[(pos1 - l) as usize..pos1 as usize];
        let vseq = &enc[(pos2 - l) as usize..pos2 as usize];
        xleft = eval_xdrop(false, useq, vseq, belowscore, scores, dist);
    }

    // right (forward) xdrop
    if alilen != 0 && pos1 + len < seqend1 && pos2 + len < seqend2 {
        let l = if alilen <= seqend1 - (pos1 + len) && alilen <= seqend2 - (pos2 + len) {
            alilen
        } else {
            (seqend1 - (pos1 + len)).min(seqend2 - (pos2 + len))
        };
        let useq = &enc[(pos1 + len) as usize..(pos1 + len + l) as usize];
        let vseq = &enc[(pos2 + len) as usize..(pos2 + len + l) as usize];
        xright = eval_xdrop(true, useq, vseq, belowscore, scores, dist);
    }
    (xleft, xright)
}
