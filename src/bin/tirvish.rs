//! tirvish_rs — full de-novo TIRvish. Runs all 5 stages, emits the gold-TSV shape
//! (seqid start stop tir1 tir2 tsd1 tsd2 sim), sorted like the oracle so it diffs
//! directly against gold/chunkN.gold.tsv.
//!
//! Single fragment (-> stdout):
//!     tirvish <fasta>
//! Batch of pre-batched fragments (-> <outdir>/<basename>.tirvish.tsv each),
//! parallel at the FRAGMENT level over the global pool:
//!     tirvish --batch <outdir> [--threads N] <frag1.fa> <frag2.fa> ...
//!     tirvish --batch <outdir> [--threads N]            # paths from stdin, one per line
//! Fragment-level parallelism runs each fragment 1/thread (100% efficient in
//! bulk); the per-seed inner par_iter only steals at the tail.

use std::io::Read;
use tirvish_rs::pipeline::{elements_tsv, run};
use tirvish_rs::{read_fasta, run_batch, set_threads};

fn main() {
    let args: Vec<String> = std::env::args().collect();

    if args.get(1).map(|s| s.as_str()) == Some("--batch") {
        let outdir = args.get(2).expect("usage: tirvish --batch <outdir> [--threads N] [paths...]");
        // optional --threads N
        let mut threads = 0usize; // 0 => rayon default (all logical cores)
        let mut rest_start = 3;
        if args.get(3).map(|s| s.as_str()) == Some("--threads") {
            threads = args.get(4).and_then(|s| s.parse().ok()).expect("--threads N");
            rest_start = 5;
        }
        let mut paths: Vec<String> = args[rest_start..].to_vec();
        if paths.is_empty() {
            // read newline-separated fragment paths from stdin (HPC manifest / xargs)
            let mut s = String::new();
            std::io::stdin().read_to_string(&mut s).expect("read stdin");
            paths = s.lines().map(str::trim).filter(|l| !l.is_empty()).map(String::from).collect();
        }
        if threads > 0 {
            set_threads(threads);
        }
        let t0 = std::time::Instant::now();
        let n = run_batch(&paths, outdir);
        eprintln!("{} fragments -> {} in {:.1}s", n, outdir, t0.elapsed().as_secs_f64());
        return;
    }

    // single-fragment mode -> stdout
    let path = args.get(1).expect("usage: tirvish <fasta>  |  tirvish --batch <outdir> [paths...]");
    let contigs = read_fasta(path);
    let mut els = run(&contigs);
    print!("{}", elements_tsv(&mut els));
    eprintln!("{} elements", els.len());
}
