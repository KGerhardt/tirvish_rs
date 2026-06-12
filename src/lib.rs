//! tirvish_rs — a faithful Rust port of `gt tirvish` (genometools 1.6.5).
//!
//! INDEPENDENT tool: must not ship with grf_rs (which is destined to become a
//! full GRF C++→Rust port). May copy/vendor grf_rs code (the user's own), but
//! not depend on it. See project memory `tirvish-rust-port-plan`.
//!
//! # Contract
//! Reproduce `gt tirvish`'s RAW prediction multiset on identical input. The
//! acceptance test is the committed oracle in
//! `TIR-Learner/tirlearner_run/tirvish_oracle/` (4 shrimp chunks + gold TSVs).
//!
//! # Phase A — direct exact port (this crate's current goal)
//! Port all five stages faithfully, validating each against the instrumented
//! reference build (gt with TIRVISH_TRACE seed+pair dumps). Faithful, not yet
//! fast — the ESA seed enumeration is kept. Phase B later swaps stage 1 for a
//! bounded inverted-MEM finder, validated against Phase A.
//!
//! ## The five stages (tir_stream.c)
//! 1. SEED ENUM — `gt_enumeratemaxpairs` (Abouelhoda–Kurtz bottom-up lcp-interval
//!    maximal pairs) over the `-mirrored` suffix array → every maximal exact
//!    repeat of len >= SEED. `store_seeds` keeps those straddling midpos, in the
//!    distance band, on mirror-corresponding contigs, len <= MAX_TIR_LEN.
//!    Built here on a Rust SA/LCP (libsais); the SA+LCP are canonical, so a
//!    verbatim maxpairs port emits the identical pair set.
//! 2. XDROP — `gt_evalxdroparbitscoresextend` left+right from each seed →
//!    arm boundaries (left arm uses i-extensions, right arm j-extensions:
//!    arms can be asymmetric). Scores: see [`params`].
//! 3. TSD — search the ±VICINITY flanks, pick best TSD; MUTATES arm boundaries.
//!    Runs BEFORE stage 4. (Explains arms shorter than the seed.)
//! 4. SIMILARITY — greedyunitedist over the post-TSD arms; gate at
//!    [`params::SIMILARITY_THRESHOLD`].
//! 5. SORT + OVERLAP REMOVAL — qsort by (contig, left_tir_start,
//!    right_transformed_start) then `remove_overlaps` in "best" mode. The sort
//!    makes final output INDEPENDENT of seed enumeration order.
//!
//! ## Why faithfulness reduces to two checks
//! TIRvish sorts before overlap removal, so output is order-independent; and
//! every seed is a maximal exact match, so reproducing the exact seed set →
//! byte-identical seeds → deterministic stages 2–5 → identical output. Hence:
//! faithful ⟺ (1) exact seed set + (2) verbatim stages 2–5. Both oracle-checked.

pub mod encode;
pub mod maxpairs;
pub mod params;
pub mod pipeline;
pub mod sa;
pub mod seeds;
pub mod similarity;
pub mod tsd;
pub mod twobit;
pub mod xdrop;

// Stage modules land here as the port progresses (Phase A):
//   pub mod encode;     // ACGT + mirror string + contig/midpos coordinate system
//   pub mod sa;         // suffix array + LCP (libsais)
//   pub mod maxpairs;   // stage 1: bottom-up maximal-pairs (verbatim port)
//   pub mod seeds;      // stage 1 filter: store_seeds
//   pub mod xdrop;      // stage 2
//   pub mod tsd;        // stage 3
//   pub mod similarity; // stage 4 (greedyunitedist)
//   pub mod overlaps;   // stage 5
//   pub mod trace;      // parse the reference TIRVISH_TRACE dumps for validation
