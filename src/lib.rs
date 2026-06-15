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

/// Read a FASTA into (header, uppercased-sequence) records via needletail.
///
/// `header` is the full record id (everything after `>` up to end of line, no
/// whitespace split) — pipeline.rs strips the TIR-Learner `;;<n>` suffix. The
/// sequence is uppercased to match gt's `-dna` handling; needletail joins
/// multi-line records and transparently handles gzip. Matches the previous
/// hand-rolled reader byte-for-byte on the oracle chunks while parsing faster.
pub fn read_fasta(path: &str) -> Vec<(String, Vec<u8>)> {
    let mut reader = needletail::parse_fastx_file(path).expect("open fasta");
    let mut out = Vec::new();
    while let Some(rec) = reader.next() {
        let rec = rec.expect("fasta record");
        let name = String::from_utf8_lossy(rec.id()).into_owned();
        let seq = rec.seq().iter().map(|b| b.to_ascii_uppercase()).collect();
        out.push((name, seq));
    }
    out
}

/// Size the global rayon pool to `threads` (call once at startup; later calls are
/// no-ops). The batch driver's outer per-fragment par_iter and pipeline::run's
/// inner per-seed par_iter both run in this one pool.
pub fn set_threads(threads: usize) {
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads.max(1))
        .build_global()
        .ok();
}

/// Batch entry: process pre-batched fragment FASTAs with parallelism at the
/// FRAGMENT level (1 fragment/worker), writing `<outdir>/<basename>.tirvish.tsv`
/// per fragment as it finishes. The per-seed inner par_iter in pipeline::run
/// nests in the same global pool, so it stays near-inert while every worker has
/// its own fragment (bulk = inter-fragment, 100% efficiency) and only splits a
/// fragment's seeds when workers go idle at the tail (intra-fragment steal).
/// Mirrors grf_rs::run_batch. Returns the number of fragments processed.
pub fn run_batch(paths: &[String], outdir: &str) -> usize {
    use rayon::prelude::*;
    paths.par_iter().for_each(|path| {
        let contigs = read_fasta(path);
        let mut els = pipeline::run(&contigs);
        let base = std::path::Path::new(path)
            .file_stem()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| "out".into());
        let out = format!("{}/{}.tirvish.tsv", outdir, base);
        std::fs::write(&out, pipeline::elements_tsv(&mut els))
            .unwrap_or_else(|e| panic!("write {out}: {e}"));
    });
    paths.len()
}

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
