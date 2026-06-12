//! 2-bit packed (mirrored) genome with SWAR/XOR longest-common-extension.
//!
//! The mechanical substrate for the comparison kernels (Xdrop, greedyunitedist,
//! later TSD). Each base is 2 bits (A=0 C=1 G=2 T=3); the mirror half already
//! stores the complement, so an inverted-repeat match is plain code equality.
//! `lce` XORs a u64 (32 bases) at a time and locates the first differing or
//! special code with `trailing_zeros`, replacing the per-base loop.
//!
//! Specials (separators / wildcards, code >= ALPHA) can't live in 2 bits, so a
//! parallel mask marks them at the same even-bit positions as the codes; a match
//! stops at the first special in either operand. A separate reverse packing makes
//! the left (reverse) Xdrop extension a plain forward LCE.

const MASK5555: u64 = 0x5555_5555_5555_5555; // one bit per 2-bit code field

pub struct TwoBit {
    fwd: Vec<u64>,
    fwd_spec: Vec<u64>,
    rev: Vec<u64>,      // codes in reverse order: rev[i] = code at (len-1-i)
    rev_spec: Vec<u64>,
    len: usize,
}

impl TwoBit {
    /// Pack from the per-position codes (`enc`: ACGT 0..3, special == alpha).
    pub fn from_enc(enc: &[u32], alpha: u32) -> Self {
        let len = enc.len();
        let nwords = len / 32 + 2; // +pad so window() can always read widx+1
        let mut fwd = vec![0u64; nwords];
        let mut fwd_spec = vec![0u64; nwords];
        let mut rev = vec![0u64; nwords];
        let mut rev_spec = vec![0u64; nwords];
        for (p, &c) in enc.iter().enumerate() {
            let (code, special) = if c < alpha { (c as u64 & 3, false) } else { (0u64, true) };
            let (w, sh) = (p / 32, 2 * (p % 32));
            fwd[w] |= code << sh;
            if special {
                fwd_spec[w] |= 1u64 << sh;
            }
            let rp = len - 1 - p;
            let (rw, rsh) = (rp / 32, 2 * (rp % 32));
            rev[rw] |= code << rsh;
            if special {
                rev_spec[rw] |= 1u64 << rsh;
            }
        }
        TwoBit { fwd, fwd_spec, rev, rev_spec, len }
    }

    /// 64-bit window = 32 codes starting at `pos`, normalized to bit 0.
    #[inline]
    fn window(words: &[u64], pos: usize) -> u64 {
        let bit = 2 * pos;
        let widx = bit / 64;
        let off = bit % 64;
        let lo = words[widx] >> off;
        let hi = if off == 0 { 0 } else { words[widx + 1] << (64 - off) };
        lo | hi
    }

    #[inline]
    fn lce_arr(words: &[u64], spec: &[u64], pa: usize, pb: usize, maxlen: usize) -> usize {
        let mut m = 0;
        while m < maxlen {
            let wa = Self::window(words, pa + m);
            let wb = Self::window(words, pb + m);
            let x = wa ^ wb;
            // bit at code k set iff codes differ
            let mism = (x | (x >> 1)) & MASK5555;
            let sa = Self::window(spec, pa + m);
            let sb = Self::window(spec, pb + m);
            let stop = mism | ((sa | sb) & MASK5555);
            if stop == 0 {
                m += 32;
            } else {
                let k = (stop.trailing_zeros() / 2) as usize;
                return (m + k).min(maxlen);
            }
        }
        maxlen
    }

    /// gt_seqabstract_lcp(forward, ...). Forward compares genome[pa+m] vs
    /// genome[pb+m]; reverse compares genome[pa-m] vs genome[pb-m] (pa/pb are the
    /// first positions of the run). Stops at the first mismatch or special.
    #[inline]
    pub fn lce(&self, forward: bool, pa: usize, pb: usize, maxlen: usize) -> usize {
        if forward {
            Self::lce_arr(&self.fwd, &self.fwd_spec, pa, pb, maxlen)
        } else {
            Self::lce_arr(
                &self.rev,
                &self.rev_spec,
                self.len - 1 - pa,
                self.len - 1 - pb,
                maxlen,
            )
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    const ALPHA: u32 = 4;

    // naive reference: gt_seqabstract_lcp over the code array.
    fn naive(enc: &[u32], forward: bool, pa: usize, pb: usize, maxlen: usize) -> usize {
        let n = enc.len();
        let mut m = 0;
        while m < maxlen {
            let (ia, ib) = if forward { (pa + m, pb + m) } else { (pa - m, pb - m) };
            if ia >= n || ib >= n { break; }
            if enc[ia] >= ALPHA || enc[ib] >= ALPHA || enc[ia] != enc[ib] { break; }
            m += 1;
        }
        m
    }

    #[test]
    fn swar_lce_matches_naive() {
        let mut st: u64 = 0xDEADBEEF12345678;
        let mut next = || { st = st.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407); (st >> 33) as u32 };
        for _ in 0..400 {
            let n = 40 + (next() % 300) as usize;
            // codes 0..3, with ~8% specials (ALPHA)
            let enc: Vec<u32> = (0..n).map(|_| { let r = next() % 50; if r < 4 { ALPHA } else { r % 4 } }).collect();
            let tb = TwoBit::from_enc(&enc, ALPHA);
            for _ in 0..40 {
                let pa = (next() as usize) % n;
                let pb = (next() as usize) % n;
                // forward
                let maxf = (n - pa).min(n - pb).min(1 + (next() as usize) % n);
                assert_eq!(tb.lce(true, pa, pb, maxf), naive(&enc, true, pa, pb, maxf),
                           "fwd pa={pa} pb={pb} maxf={maxf}");
                // reverse
                let maxr = (pa + 1).min(pb + 1).min(1 + (next() as usize) % n);
                assert_eq!(tb.lce(false, pa, pb, maxr), naive(&enc, false, pa, pb, maxr),
                           "rev pa={pa} pb={pb} maxr={maxr}");
            }
        }
    }
}
