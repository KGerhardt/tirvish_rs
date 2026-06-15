//! TIRvish parameters — exact values for gt 1.6.5 (verified: tir_stream.c,
//! esa-maxpairs.c, gt_tir.c, esa-bottomup-maxpairs.inc are byte-identical
//! between the v1.6.5 tag and the HEAD we read/instrumented).
//!
//! Values come from two places:
//!   - the TIR-Learner pipeline command (tirvish_new.one_tirvish), and
//!   - gt_tir.c option defaults, for any flag the pipeline does NOT pass.
//! The instrumented reference build is run with these same flags.

/// Min seed length = min length of an exact maximal repeat fed to stage 2.
/// gt_tir.c default 20 (min 5); pipeline `-seed 20`.
pub const SEED: u64 = 20;

/// Per-arm TIR length bounds. Pipeline `-mintirlen 10 -maxtirlen 1000`.
/// store_seeds rejects seeds with len > MAX_TIR_LEN (tir_stream.c:168); the
/// final per-arm length is re-checked after Xdrop extension (tir_stream.c:571).
pub const MIN_TIR_LEN: u64 = 10;
pub const MAX_TIR_LEN: u64 = 1000;

/// TIR distance (the spacer between the two arms). Pipeline
/// `-mintirdist 10 -maxtirdist 5000`. Checked on the SEED in store_seeds:
/// distance = (REVERSEPOS(T,pos2) - len + 1) - pos1  (tir_stream.c:157).
pub const MIN_TIR_DIST: u64 = 10;
pub const MAX_TIR_DIST: u64 = 5000;

/// TSD length bounds. Pipeline `-mintsd 2 -maxtsd 11`. TSD search only runs
/// when MIN_TSD_LEN > 1 (tir_stream.c:415).
pub const MIN_TSD_LEN: u64 = 2;
pub const MAX_TSD_LEN: u64 = 11;

/// Vicinity: ± bp around each 5'/3' boundary searched for the best TSD.
/// Pipeline `-vic 13`. The TSD search MUTATES the arm boundaries within this
/// window BEFORE the similarity gate (tir_stream.c:599 precedes :607).
pub const VICINITY: u64 = 13;

/// TIR similarity threshold (%). Pipeline `-similar 80` (gt_tir.c default 85).
/// sim = 100*(1 - edist/max(ulen,vlen)) via greedyunitedist over the POST-TSD
/// arms; pair dropped if sim < threshold (tir_stream.c:617-622).
pub const SIMILARITY_THRESHOLD: f64 = 80.0;

/// Xdrop arbitrary scores — gt_tir.c DEFAULTS (pipeline overrides none of these).
/// Constraints (xdrop.h:55): mat >= mis, mat >= 2*ins, mat >= 2*del.
pub const XDROP_MAT: i32 = 2; // gt_tir.c:183
pub const XDROP_MIS: i32 = -2; // gt_tir.c:190
pub const XDROP_INS: i32 = -3; // gt_tir.c:197
pub const XDROP_DEL: i32 = -3; // gt_tir.c:205

/// Xdrop below-score (greedy-extension drop threshold). gt_tir.c default 5.
pub const XDROP_BELOWSCORE: i64 = 5;

/// Overlap-removal mode. gt_tir.c default "best" (pipeline passes no -overlaps):
/// among overlapping pairs keep max-similarity, then max-length
/// (gt_tir_remove_overlaps, tir_stream.c). Output is fully non-overlapping.
pub const OVERLAPS_BEST: bool = true;

/// REVERSEPOS over the mirrored encseq (encseq.h:42): maps a mirror-half
/// position back to its forward counterpart. T = total length of the mirrored
/// encseq; midpos = (T-1)/2 (tir_stream.c:1051). A cross-mirror seed has
/// pos1 <= midpos <= pos2 (tir_stream.c:153).
#[inline]
pub const fn reversepos(total_len: u64, pos: u64) -> u64 {
    total_len - 1 - pos
}

/// Runtime TIRvish parameters (the gt tirvish options). `Default` is the locked
/// TIR-Learner value set (the consts above); the CLI overrides individual fields.
#[derive(Clone, Copy, Debug)]
pub struct Params {
    pub seed: u64,
    pub min_tir_len: u64,
    pub max_tir_len: u64,
    pub min_tir_dist: u64,
    pub max_tir_dist: u64,
    pub min_tsd_len: u64,
    pub max_tsd_len: u64,
    pub vicinity: u64,
    pub similar: f64,
    /// Precomputed (100 - similar)/100, so the similarity band in the hot loop is
    /// `max * sim_mult + 2` rather than recomputing the fraction per call.
    pub sim_mult: f64,
    pub xdrop_mat: i32,
    pub xdrop_mis: i32,
    pub xdrop_ins: i32,
    pub xdrop_del: i32,
    pub xdrop_belowscore: i64,
}

impl Default for Params {
    fn default() -> Self {
        Params {
            seed: SEED,
            min_tir_len: MIN_TIR_LEN,
            max_tir_len: MAX_TIR_LEN,
            min_tir_dist: MIN_TIR_DIST,
            max_tir_dist: MAX_TIR_DIST,
            min_tsd_len: MIN_TSD_LEN,
            max_tsd_len: MAX_TSD_LEN,
            vicinity: VICINITY,
            similar: SIMILARITY_THRESHOLD,
            sim_mult: (100.0 - SIMILARITY_THRESHOLD) / 100.0,
            xdrop_mat: XDROP_MAT,
            xdrop_mis: XDROP_MIS,
            xdrop_ins: XDROP_INS,
            xdrop_del: XDROP_DEL,
            xdrop_belowscore: XDROP_BELOWSCORE,
        }
    }
}

impl Params {
    /// Set the similarity threshold and keep `sim_mult` consistent.
    pub fn set_similar(&mut self, similar: f64) {
        self.similar = similar;
        self.sim_mult = (100.0 - similar) / 100.0;
    }
}
