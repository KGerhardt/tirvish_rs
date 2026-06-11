//! Stage 1 — maximal-pairs enumeration.
//!
//! Verbatim port of genometools 1.6.5 `esa-maxpairs.c` (the Abouelhoda–Kurtz /
//! REPuter callbacks) + `esa-bottomup-maxpairs.inc` (the stack-based LCP-interval
//! bottom-up traversal driver). Emits every maximal exact repeat
//! `(len, pos1, pos2)` with `pos1 < pos2` and `len >= searchlength`.
//!
//! Inputs are the suffix array and LCP table over the encoded (mirrored)
//! sequence, plus per-position left-character codes. The SA+LCP are canonical,
//! so this emits the identical pair set as gt's `gt_enumeratemaxpairs`.
//!
//! `enc[pos]`: left-character code at `pos`. Regular symbols are `0..alphabetsize`;
//! ANY special (separator / wildcard / contig boundary) is encoded as exactly
//! `alphabetsize` (== ISLEFTDIVERSE), which is all the left-maximality logic needs
//! (`add2poslist` routes it to the unique list; `CHECKCHAR` treats it as diverse).

/// `commonchar` sentinels, parameterized by alphabet size (gt: macros at
/// esa-maxpairs.c:27-28). ISLEFTDIVERSE marks a left-diverse (left-maximal) node.
#[inline]
fn is_left_diverse(alphabetsize: u32) -> u32 {
    alphabetsize
}
#[inline]
fn initial_char(alphabetsize: u32) -> u32 {
    alphabetsize + 1
}

#[derive(Clone)]
struct Listtype {
    start: usize,
    length: usize,
}

/// Per lcp-interval-tree node info (gt: `GtBUinfo_maxpairs`). Holds only offsets
/// into the shared `State` buffers, never the positions themselves.
#[derive(Clone)]
struct BuInfo {
    commonchar: u32,
    uniquecharposstart: usize,
    uniquecharposlength: usize,
    nodeposlist: Vec<Listtype>, // length == alphabetsize
}

impl BuInfo {
    fn new(alphabetsize: u32) -> Self {
        BuInfo {
            commonchar: 0,
            uniquecharposstart: 0,
            uniquecharposlength: 0,
            nodeposlist: vec![Listtype { start: 0, length: 0 }; alphabetsize as usize],
        }
    }
}

/// Global traversal state (gt: `GtBUstate_maxpairs`). The `poslist`/`uniquechar`
/// buffers are shared across all live nodes; nodes reference ranges via offsets.
struct State {
    searchlength: u64,
    alphabetsize: u32,
    poslist: Vec<Vec<u64>>, // one append-buffer per base
    uniquechar: Vec<u64>,
    initialized: bool,
}

impl State {
    fn new(searchlength: u64, alphabetsize: u32) -> Self {
        State {
            searchlength,
            alphabetsize,
            poslist: vec![Vec::new(); alphabetsize as usize],
            uniquechar: Vec::new(),
            initialized: false,
        }
    }

    /// gt: `setpostabto0_maxpairs` — reset the shared buffers, but only once per
    /// below-searchlength gap (guarded by `initialized`).
    fn setpostabto0(&mut self) {
        if !self.initialized {
            for base in 0..self.alphabetsize as usize {
                self.poslist[base].clear();
            }
            self.uniquechar.clear();
            self.initialized = true;
        }
    }

    /// gt: `add2poslist_maxpairs`. `base >= alphabetsize` => the unique list.
    fn add2poslist(&mut self, ninfo: &mut BuInfo, base: u32, leafnumber: u64) {
        if base >= self.alphabetsize {
            ninfo.uniquecharposlength += 1;
            self.uniquechar.push(leafnumber);
        } else {
            self.poslist[base as usize].push(leafnumber);
            ninfo.nodeposlist[base as usize].length += 1;
        }
    }
}

/// gt: `concatlists_maxpairs` — fold son's per-base lengths into father.
fn concatlists(alphabetsize: u32, father: &mut BuInfo, son: &BuInfo) {
    for base in 0..alphabetsize as usize {
        father.nodeposlist[base].length += son.nodeposlist[base].length;
    }
    father.uniquecharposlength += son.uniquecharposlength;
}

/// gt: `cartproduct1_maxpairs` — pair `leafnumber` with every position in
/// `ninfo`'s `base` list. Emits maximal pairs at depth `fatherdepth`.
fn cartproduct1(
    state: &State,
    fatherdepth: u64,
    ninfo: &BuInfo,
    base: u32,
    leafnumber: u64,
    process: &mut impl FnMut(u64, u64, u64),
) {
    let pl = &ninfo.nodeposlist[base as usize];
    let bucket = &state.poslist[base as usize];
    for &spval in &bucket[pl.start..pl.start + pl.length] {
        process(fatherdepth, leafnumber.min(spval), leafnumber.max(spval));
    }
}

/// gt: `cartproduct2_maxpairs` — cross product of two per-base lists.
fn cartproduct2(
    state: &State,
    fatherdepth: u64,
    ninfo1: &BuInfo,
    base1: u32,
    ninfo2: &BuInfo,
    base2: u32,
    process: &mut impl FnMut(u64, u64, u64),
) {
    let pl1 = &ninfo1.nodeposlist[base1 as usize];
    let pl2 = &ninfo2.nodeposlist[base2 as usize];
    let b1 = &state.poslist[base1 as usize];
    let b2 = &state.poslist[base2 as usize];
    for &v1 in &b1[pl1.start..pl1.start + pl1.length] {
        for &v2 in &b2[pl2.start..pl2.start + pl2.length] {
            process(fatherdepth, v1.min(v2), v1.max(v2));
        }
    }
}

/// gt: CHECKCHAR macro (esa-maxpairs.c:30) — left-maximality update.
#[inline]
fn checkchar(father: &mut BuInfo, cc: u32, alphabetsize: u32) {
    if father.commonchar != cc || cc >= is_left_diverse(alphabetsize) {
        father.commonchar = is_left_diverse(alphabetsize);
    }
}

/// gt: `processleafedge_maxpairs`.
fn processleafedge(
    firstsucc: bool,
    fatherdepth: u64,
    father: &mut BuInfo,
    leafnumber: u64,
    state: &mut State,
    enc: &[u32],
    process: &mut impl FnMut(u64, u64, u64),
) {
    if fatherdepth < state.searchlength {
        state.setpostabto0();
        return;
    }
    let leftchar = if leafnumber == 0 {
        initial_char(state.alphabetsize)
    } else {
        enc[(leafnumber - 1) as usize]
    };
    state.initialized = false;
    if firstsucc {
        father.commonchar = leftchar;
        father.uniquecharposlength = 0;
        father.uniquecharposstart = state.uniquechar.len();
        for base in 0..state.alphabetsize as usize {
            father.nodeposlist[base].start = state.poslist[base].len();
            father.nodeposlist[base].length = 0;
        }
        state.add2poslist(father, leftchar, leafnumber);
        return;
    }
    if father.commonchar != is_left_diverse(state.alphabetsize) {
        checkchar(father, leftchar, state.alphabetsize);
    }
    if father.commonchar == is_left_diverse(state.alphabetsize) {
        for base in 0..state.alphabetsize {
            if leftchar != base {
                cartproduct1(state, fatherdepth, father, base, leafnumber, process);
            }
        }
        let start = father.uniquecharposstart;
        for i in start..start + father.uniquecharposlength {
            let spval = state.uniquechar[i];
            process(fatherdepth, leafnumber.min(spval), leafnumber.max(spval));
        }
    }
    state.add2poslist(father, leftchar, leafnumber);
}

/// gt: `processbranchingedge_maxpairs`. `son` is None only when `firstsucc`
/// (gt passes NULL there and returns before dereferencing).
fn processbranchingedge(
    firstsucc: bool,
    fatherdepth: u64,
    father: &mut BuInfo,
    son: Option<&BuInfo>,
    state: &mut State,
    process: &mut impl FnMut(u64, u64, u64),
) {
    if fatherdepth < state.searchlength {
        state.setpostabto0();
        return;
    }
    // maxfreqcollect is NULL in the tirvish path (no -maxfreq) -> skip.
    state.initialized = false;
    if firstsucc {
        return;
    }
    let son = son.expect("non-first branching edge must have a son");
    let ld = is_left_diverse(state.alphabetsize);
    if father.commonchar != ld {
        if son.commonchar != ld {
            checkchar(father, son.commonchar, state.alphabetsize);
        } else {
            father.commonchar = ld;
        }
    }
    if father.commonchar == ld {
        let sstart = son.uniquecharposstart;
        let slen = son.uniquecharposlength;
        for chfather in 0..state.alphabetsize {
            for chson in 0..state.alphabetsize {
                if chson != chfather {
                    cartproduct2(state, fatherdepth, father, chfather, son, chson, process);
                }
            }
            for i in sstart..sstart + slen {
                let spval = state.uniquechar[i];
                cartproduct1(state, fatherdepth, father, chfather, spval, process);
            }
        }
        let fstart = father.uniquecharposstart;
        let flen = father.uniquecharposlength;
        for fi in fstart..fstart + flen {
            let fval = state.uniquechar[fi];
            for chson in 0..state.alphabetsize {
                cartproduct1(state, fatherdepth, son, chson, fval, process);
            }
            for i in sstart..sstart + slen {
                let spval = state.uniquechar[i];
                process(fatherdepth, fval.min(spval), fval.max(spval));
            }
        }
    }
    concatlists(state.alphabetsize, father, son);
}

/// One interval on the bottom-up stack (gt: `GtBUItvinfo_maxpairs`).
struct ItvInfo {
    lcp: u64,
    lb: u64,
    #[allow(dead_code)]
    rb: u64,
    info: BuInfo,
}

/// Driver — verbatim port of `gt_esa_bottomup_maxpairs` (esa-bottomup-maxpairs.inc).
///
/// `suftab[idx]` / `lcptab[idx]` are the suffix and LCP at SA rank `idx`
/// (`lcptab[0] == 0`); `n = suftab.len()` = number of non-special suffixes.
/// `enc` = per-position left-char codes (specials == `alphabetsize`).
pub fn enumerate_maxpairs(
    suftab: &[u64],
    lcptab: &[u64],
    searchlength: u64,
    alphabetsize: u32,
    enc: &[u32],
    mut process: impl FnMut(u64, u64, u64),
) {
    let n = suftab.len();
    let mut state = State::new(searchlength, alphabetsize);
    let mut stack: Vec<ItvInfo> = Vec::new();
    // PUSH_ESA_BOTTOMUP(0,0)
    stack.push(ItvInfo {
        lcp: 0,
        lb: 0,
        rb: u64::MAX,
        info: BuInfo::new(alphabetsize),
    });
    let mut firstedgefromroot = true;
    let mut lastinterval: Option<ItvInfo> = None;

    for idx in 0..n {
        let previoussuffix = suftab[idx];
        // LCP is LOOK-AHEAD: lcp(suftab[idx], suftab[idx+1]) == standard
        // look-behind lcptab[idx+1] (0 past the end). This is gt's streaming
        // convention; using look-behind lcptab[idx] drops the first leaf of
        // each interval. Verified against a brute-force maximal-pairs oracle.
        let lcpvalue = if idx + 1 < n { lcptab[idx + 1] } else { 0 };

        if lcpvalue <= stack.last().unwrap().lcp {
            let firstedge;
            if stack.last().unwrap().lcp > 0 || !firstedgefromroot {
                firstedge = false;
            } else {
                firstedge = true;
                firstedgefromroot = false;
            }
            let top = stack.last_mut().unwrap();
            processleafedge(
                firstedge, top.lcp, &mut top.info, previoussuffix, &mut state, enc,
                &mut process,
            );
        }
        debug_assert!(lastinterval.is_none());
        while lcpvalue < stack.last().unwrap().lcp {
            let mut li = stack.pop().unwrap();
            li.rb = idx as u64;
            lastinterval = Some(li);
            if lcpvalue <= stack.last().unwrap().lcp {
                let firstedge;
                if stack.last().unwrap().lcp > 0 || !firstedgefromroot {
                    firstedge = false;
                } else {
                    firstedge = true;
                    firstedgefromroot = false;
                }
                let son = lastinterval.take().unwrap();
                let top = stack.last_mut().unwrap();
                processbranchingedge(
                    firstedge, top.lcp, &mut top.info, Some(&son.info), &mut state,
                    &mut process,
                );
            }
        }
        if lcpvalue > stack.last().unwrap().lcp {
            if let Some(li) = lastinterval.take() {
                let lastintervallb = li.lb;
                // gt reuses the just-popped stack SLOT, so the new deeper node
                // INHERITS the popped interval's info (its accumulated poslist
                // ranges). processbranchingedge(firstsucc=true) returns without
                // touching it, so the inheritance is what carries those leaves.
                stack.push(ItvInfo {
                    lcp: lcpvalue,
                    lb: lastintervallb,
                    rb: u64::MAX,
                    info: li.info,
                });
                let top = stack.last_mut().unwrap();
                processbranchingedge(
                    true, top.lcp, &mut top.info, None, &mut state, &mut process,
                );
            } else {
                stack.push(ItvInfo {
                    lcp: lcpvalue,
                    lb: idx as u64,
                    rb: u64::MAX,
                    info: BuInfo::new(alphabetsize),
                });
                let top = stack.last_mut().unwrap();
                processleafedge(
                    true, top.lcp, &mut top.info, previoussuffix, &mut state, enc,
                    &mut process,
                );
            }
        }
    }
    // No post-loop final edge: the look-ahead form feeds lcpvalue=0 at idx=n-1,
    // which processes the last suffix and pops all open intervals in-loop.
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    // Longest common extension over `enc` (sentinel value `alphabetsize` is unique
    // per its single occurrence, so it terminates any match).
    fn lce(enc: &[u32], a: u64, b: u64) -> u64 {
        let (mut i, mut j) = (a as usize, b as usize);
        let mut l = 0u64;
        while i < enc.len() && j < enc.len() && enc[i] < ALPHA && enc[j] < ALPHA && enc[i] == enc[j] {
            i += 1; j += 1; l += 1;
        }
        l
    }
    const ALPHA: u32 = 4; // DNA

    // Build SA over the ACGT-start positions (0..m), ordered by their enc suffix
    // (sentinel at index m makes all distinct & terminates matches), + LCP table.
    fn naive_sa_lcp(enc: &[u32], m: usize) -> (Vec<u64>, Vec<u64>) {
        let mut sa: Vec<u64> = (0..m as u64).collect();
        sa.sort_by(|&a, &b| {
            let (mut i, mut j) = (a as usize, b as usize);
            loop {
                let (ca, cb) = (enc[i], enc[j]);
                if ca != cb { return ca.cmp(&cb); }
                i += 1; j += 1;
            }
        });
        let mut lcp = vec![0u64; m];
        for k in 1..m { lcp[k] = lce(enc, sa[k - 1], sa[k]); }
        (sa, lcp)
    }

    // Reference: every left- AND right-maximal exact pair (len, lo, hi), len>=sl.
    fn brute(enc: &[u32], m: usize, sl: u64) -> HashSet<(u64, u64, u64)> {
        let mut out = HashSet::new();
        for i in 0..m as u64 {
            for j in (i + 1)..m as u64 {
                let l = lce(enc, i, j);
                if l < sl { continue; }
                let leftmax = i == 0 || j == 0 || enc[i as usize - 1] != enc[j as usize - 1];
                if leftmax { out.insert((l, i, j)); }
            }
        }
        out
    }

    fn run_case(seq: &[u32], sl: u64) {
        let m = seq.len();
        let mut enc = seq.to_vec();
        enc.push(ALPHA); // trailing unique sentinel
        let (sa, lcp) = naive_sa_lcp(&enc, m);
        let mut got: HashSet<(u64, u64, u64)> = HashSet::new();
        enumerate_maxpairs(&sa, &lcp, sl, ALPHA, &enc, |len, p1, p2| {
            got.insert((len, p1, p2));
        });
        let want = brute(&enc, m, sl);
        assert_eq!(got, want, "mismatch on seq {:?} sl={}", seq, sl);
    }

    #[test]
    fn hand_examples() {
        run_case(&[0, 1, 2, 3, 0, 1, 2, 3], 2); // ACGTACGT direct repeat
        run_case(&[0, 0, 0, 0, 0], 2); // AAAAA (overlapping repeats)
        run_case(&[0, 1, 0, 1, 0, 1, 2], 2);
        run_case(&[3, 3, 0, 1, 3, 3, 0, 1, 3, 3], 3);
    }

    #[test]
    fn random_strings() {
        // Deterministic LCG; no Math.random/std rand dependency.
        let mut st: u64 = 0x9E3779B97F4A7C15;
        let mut next = || { st = st.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407); (st >> 33) as u32 };
        for _ in 0..300 {
            let n = 5 + (next() % 60) as usize;
            let seq: Vec<u32> = (0..n).map(|_| next() % 4).collect();
            for sl in [2u64, 3, 5] { run_case(&seq, sl); }
        }
    }
}
