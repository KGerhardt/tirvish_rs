//! tirvish_rs CLI (Phase A scaffold). Pipeline wiring lands as stages are ported.
use tirvish_rs::params;

fn main() {
    eprintln!("tirvish_rs {} — Phase A scaffold", env!("CARGO_PKG_VERSION"));
    eprintln!(
        "params: seed={} tirlen[{},{}] tirdist[{},{}] tsd[{},{}] vic={} similar={} \
         xdrop(mat={},mis={},ins={},del={},below={})",
        params::SEED,
        params::MIN_TIR_LEN, params::MAX_TIR_LEN,
        params::MIN_TIR_DIST, params::MAX_TIR_DIST,
        params::MIN_TSD_LEN, params::MAX_TSD_LEN,
        params::VICINITY, params::SIMILARITY_THRESHOLD,
        params::XDROP_MAT, params::XDROP_MIS, params::XDROP_INS,
        params::XDROP_DEL, params::XDROP_BELOWSCORE,
    );
}
