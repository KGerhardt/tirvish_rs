//! Build gt's `-mirrored` encseq layout for SA construction + maxpairs.
//!
//! Layout (from genometools encseq.c:264 and tir_stream.c:1049-1051), for a
//! forward text F of length L (= concatenated contigs, internal separators
//! between them):
//!   positions [0, L)        = forward F
//!   position  L             = middle separator
//!   positions [L+1, 2L+1)   = reverse complement of F
//! so a mirror position p holds complement(F[2L - p]); total_logical = 2L+1;
//! midpos = (2L+1-1)/2 = L. A maximal exact repeat with pos1 in the forward
//! half and pos2 in the mirror half is an inverted repeat in F.
//!
//! Two encodings are produced:
//!   - `sa_input` (i32, for libsais): ACGT = 0..3, EVERY special (separator,
//!     wildcard/N, middle separator, their mirror images) gets a UNIQUE integer
//!     >= 4, so the suffix-array LCP caps at the first special and special-
//!     starting suffixes sort after all ACGT suffixes.
//!   - `enc` (u32, for maxpairs left-char): ACGT = 0..3, any special = ALPHA(4).

pub const ALPHA: u32 = 4; // DNA alphabet size (gt: gt_alphabet_num_of_chars == 4)

/// ACGT -> 0..3 (case-insensitive); anything else -> None (a special).
#[inline]
pub fn code(b: u8) -> Option<u32> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

pub struct Encoded {
    pub sa_input: Vec<i32>,   // libsais text: ACGT 0..3, specials unique >= 4
    pub enc: Vec<u32>,        // maxpairs left-char: ACGT 0..3, special == ALPHA
    pub k: i32,               // alphabet size for libsais (max symbol + 1)
    pub total_orig: u64,      // L (forward length incl. internal separators)
    pub total_logical: u64,   // T = 2L + 1
    pub midpos: u64,          // (T-1)/2 == L
    pub num_contigs: u64,     // mirrored count = 2 * input contigs
    pub seqnum_of: Vec<u32>,  // per-position seqnum; u32::MAX at separators
    pub fwd_seqstart: Vec<u64>, // forward contig start positions (len = input contigs)
    pub fwd_seqlen: Vec<u64>,   // forward contig lengths
    pub twobit: crate::twobit::TwoBit, // 2-bit packed mirror genome for SWAR LCE
}

impl Encoded {
    /// gt's seqstart1/seqend1/seqstart2/seqend2 for the forward contig `seqnum1`
    /// (tir_stream.c:493-498). seqend* are exclusive ends.
    pub fn contig_bounds(&self, seqnum1: u32) -> (u64, u64, u64, u64) {
        let seqstart1 = self.fwd_seqstart[seqnum1 as usize];
        let seqend1 = seqstart1 + self.fwd_seqlen[seqnum1 as usize];
        let seqstart2 = crate::params::reversepos(self.total_logical, seqend1);
        let seqend2 = crate::params::reversepos(self.total_logical, seqstart1);
        (seqstart1, seqend1, seqstart2, seqend2)
    }
}

impl Encoded {
    /// Number of non-special (ACGT) suffixes = the suftab/lcptab prefix length.
    pub fn num_suffixes(&self) -> usize {
        self.enc.iter().filter(|&&c| c < ALPHA).count()
    }
}

pub fn encode(contigs: &[(String, Vec<u8>)]) -> Encoded {
    let n = contigs.len() as u64;
    let mut sa_input: Vec<i32> = Vec::new();
    let mut enc: Vec<u32> = Vec::new();
    let mut seqnum_of: Vec<u32> = Vec::new();
    let mut fwd_seqstart: Vec<u64> = Vec::new();
    let mut fwd_seqlen: Vec<u64> = Vec::new();
    let mut next_sentinel: i32 = ALPHA as i32;

    // ---- forward half: contigs joined by unique-sentinel separators ----
    for (ci, (_, seq)) in contigs.iter().enumerate() {
        if ci > 0 {
            sa_input.push(next_sentinel);
            next_sentinel += 1;
            enc.push(ALPHA);
            seqnum_of.push(u32::MAX); // internal separator
        }
        fwd_seqstart.push(sa_input.len() as u64);
        fwd_seqlen.push(seq.len() as u64);
        for &b in seq {
            match code(b) {
                Some(c) => {
                    sa_input.push(c as i32);
                    enc.push(c);
                }
                None => {
                    sa_input.push(next_sentinel);
                    next_sentinel += 1;
                    enc.push(ALPHA);
                }
            }
            seqnum_of.push(ci as u32); // wildcards still belong to their contig
        }
    }
    let total_orig = sa_input.len() as u64;

    // ---- middle separator ----
    sa_input.push(next_sentinel);
    next_sentinel += 1;
    enc.push(ALPHA);
    seqnum_of.push(u32::MAX);

    // ---- mirror half: reverse complement of the forward half ----
    for p in (0..total_orig as usize).rev() {
        let e = enc[p];
        if e < ALPHA {
            let comp = 3 - e; // A<->T, C<->G via code complement
            sa_input.push(comp as i32);
            enc.push(comp);
        } else {
            sa_input.push(next_sentinel);
            next_sentinel += 1;
            enc.push(ALPHA);
        }
        let s = seqnum_of[p];
        // mirror contig of forward seqnum s is (2n-1 - s); separators stay MAX.
        seqnum_of.push(if s == u32::MAX { u32::MAX } else { (2 * n - 1 - s as u64) as u32 });
    }

    let total_logical = sa_input.len() as u64;
    let midpos = (total_logical - 1) / 2;
    let twobit = crate::twobit::TwoBit::from_enc(&enc, ALPHA);
    Encoded {
        sa_input,
        enc,
        k: next_sentinel,
        total_orig,
        total_logical,
        midpos,
        num_contigs: 2 * n,
        seqnum_of,
        fwd_seqstart,
        fwd_seqlen,
        twobit,
    }
}
