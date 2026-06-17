//! Stage 1 filter — verbatim port of `gt_tir_store_seeds` (tir_stream.c:135-178).
//!
//! For each maximal pair (len, pos1, pos2) from stage 1, keep it iff:
//!   1. it straddles the mirror midpoint: pos1 <= midpos <= pos2;
//!   2. distance = (REVERSEPOS(T, pos2) - len + 1) - pos1 is in [mintirdist, maxtirdist]
//!      (computed in unsigned wrap arithmetic, exactly as gt: an underflow wraps
//!      to a huge value and fails the upper bound);
//!   3. the two ends are on mirror-corresponding contigs:
//!      seqnum(pos2) == num_contigs - 1 - seqnum(pos1);
//!   4. len <= maxtirlen.

use crate::params::reversepos;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Seed {
    // Positions index the mirrored text (< 2^31) and len <= maxtirlen, so u32 fits
    // (16 B vs 40 B). `distance` is computed in store_seed for the filter but never
    // read afterward, so it is not stored.
    pub pos1: u32,
    pub pos2: u32,
    pub len: u32,
    pub contignumber: u32,
}

#[allow(clippy::too_many_arguments)]
pub fn store_seed(
    seeds: &mut Vec<Seed>,
    len: u64,
    pos1: u64,
    pos2: u64,
    midpos: u64,
    total_logical: u64,
    num_contigs: u64,
    seqnum_of: &[u32],
    min_tir_distance: u64,
    max_tir_distance: u64,
    max_tir_length: u64,
) {
    // gt asserts pos1 < pos2 (maxpairs emits min,max).
    // 1. mirrored vs. unmirrored
    if pos1 > midpos || pos2 < midpos {
        return;
    }
    // 2. distance constraints (unsigned, gt's evaluation order)
    let distance = reversepos(total_logical, pos2)
        .wrapping_sub(len)
        .wrapping_add(1)
        .wrapping_sub(pos1);
    if distance < min_tir_distance || distance > max_tir_distance {
        return;
    }
    // 3. same (mirror-corresponding) contig
    let seqnum1 = seqnum_of[pos1 as usize];
    let seqnum2 = seqnum_of[pos2 as usize];
    if seqnum2 as u64 != num_contigs - 1 - seqnum1 as u64 {
        return;
    }
    // 4. length constraint
    if len > max_tir_length {
        return;
    }
    seeds.push(Seed {
        pos1: pos1 as u32,
        pos2: pos2 as u32,
        len: len as u32,
        contignumber: seqnum1,
    });
}
