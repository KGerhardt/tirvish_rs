//! tirvish_rs — full de-novo TIRvish. Runs all 5 stages, emits the gold-TSV shape
//! (seqid start stop tir1 tir2 tsd1 tsd2 sim), sorted like the oracle so it diffs
//! directly against gold/chunkN.gold.tsv.
//!
//! Accepts the gt tirvish parameters (it operates directly on a FASTA, so no
//! pre-built mirrored index is needed):
//!   -seed -mintirlen -maxtirlen -mintirdist -maxtirdist -similar
//!   -mintsd -maxtsd -vic -xdrop -mat -mis -ins -del
//! All default to the locked TIR-Learner values.
//!
//! Single fragment (-> stdout):
//!     tirvish [options] <fasta>
//! Batch of pre-batched fragments (-> <outdir>/<basename>.tirvish.tsv each),
//! parallel at the FRAGMENT level over the global pool:
//!     tirvish --batch <outdir> [--threads N] [options] <frag1.fa> ...
//!     tirvish --batch <outdir> [--threads N] [options]        # paths from stdin

use std::io::Read;
use tirvish_rs::params::Params;
use tirvish_rs::pipeline::{elements_tsv, run};
use tirvish_rs::{read_fasta, run_batch, set_threads};

fn nextv(args: &[String], i: &mut usize) -> String {
    *i += 1;
    args.get(*i).cloned().unwrap_or_else(|| {
        eprintln!("missing value for {}", args[*i - 1]);
        std::process::exit(2);
    })
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mut p = Params::default();
    let mut batch = false;
    let mut outdir = String::new();
    let mut threads = 0usize; // 0 => rayon default (all logical cores)
    let mut paths: Vec<String> = Vec::new();

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--batch" => { batch = true; outdir = nextv(&args, &mut i); }
            "--threads" => threads = nextv(&args, &mut i).parse().expect("--threads N"),
            "-seed" => p.seed = nextv(&args, &mut i).parse().expect("-seed int"),
            "-mintirlen" => p.min_tir_len = nextv(&args, &mut i).parse().expect("-mintirlen int"),
            "-maxtirlen" => p.max_tir_len = nextv(&args, &mut i).parse().expect("-maxtirlen int"),
            "-mintirdist" => p.min_tir_dist = nextv(&args, &mut i).parse().expect("-mintirdist int"),
            "-maxtirdist" => p.max_tir_dist = nextv(&args, &mut i).parse().expect("-maxtirdist int"),
            "-mintsd" => p.min_tsd_len = nextv(&args, &mut i).parse().expect("-mintsd int"),
            "-maxtsd" => p.max_tsd_len = nextv(&args, &mut i).parse().expect("-maxtsd int"),
            "-vic" => p.vicinity = nextv(&args, &mut i).parse().expect("-vic int"),
            "-similar" => p.set_similar(nextv(&args, &mut i).parse().expect("-similar float")),
            "-xdrop" => p.xdrop_belowscore = nextv(&args, &mut i).parse().expect("-xdrop int"),
            "-mat" => p.xdrop_mat = nextv(&args, &mut i).parse().expect("-mat int"),
            "-mis" => p.xdrop_mis = nextv(&args, &mut i).parse().expect("-mis int"),
            "-ins" => p.xdrop_ins = nextv(&args, &mut i).parse().expect("-ins int"),
            "-del" => p.xdrop_del = nextv(&args, &mut i).parse().expect("-del int"),
            other => paths.push(other.to_string()),
        }
        i += 1;
    }

    if batch {
        if paths.is_empty() {
            // newline-separated fragment paths from stdin (HPC manifest / xargs)
            let mut s = String::new();
            std::io::stdin().read_to_string(&mut s).expect("read stdin");
            paths = s.lines().map(str::trim).filter(|l| !l.is_empty()).map(String::from).collect();
        }
        if threads > 0 {
            set_threads(threads);
        }
        let t0 = std::time::Instant::now();
        let n = run_batch(&paths, &outdir, &p);
        eprintln!("{} fragments -> {} in {:.1}s", n, outdir, t0.elapsed().as_secs_f64());
        return;
    }

    let path = paths.first().expect("usage: tirvish [options] <fasta>  |  tirvish --batch <outdir> [options] [paths...]");
    let contigs = read_fasta(path);
    let mut els = run(&contigs, &p);
    print!("{}", elements_tsv(&mut els));
    eprintln!("{} elements", els.len());
}
