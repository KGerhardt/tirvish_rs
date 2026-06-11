//! Stage-1 validation harness: read a FASTA, build the mirror SA/LCP, run
//! maxpairs + store_seeds, and report the seed count (and optionally dump the
//! seed tuples for tuple-level diff against the gt TIRVISH_TRACE `SEED` lines).
//!
//! Usage: seedcount <fasta> [--dump]
//! Params are gt tirvish's: seed=20, mintirdist=10, maxtirdist=5000, maxtirlen=1000.

use std::fs;
use tirvish_rs::encode::{encode, ALPHA};
use tirvish_rs::maxpairs::enumerate_maxpairs;
use tirvish_rs::params;
use tirvish_rs::sa::sa_lcp;
use tirvish_rs::seeds::{store_seed, Seed};

fn read_fasta(path: &str) -> Vec<(String, Vec<u8>)> {
    let data = fs::read_to_string(path).expect("read fasta");
    let mut out = Vec::new();
    let mut name = String::new();
    let mut seq: Vec<u8> = Vec::new();
    for line in data.lines() {
        if let Some(h) = line.strip_prefix('>') {
            if !name.is_empty() {
                out.push((std::mem::take(&mut name), std::mem::take(&mut seq)));
            }
            name = h.to_string();
        } else {
            seq.extend(line.bytes().map(|b| b.to_ascii_uppercase()));
        }
    }
    if !name.is_empty() {
        out.push((name, seq));
    }
    out
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let path = &args[1];
    let dump = args.iter().any(|a| a == "--dump");

    let contigs = read_fasta(path);
    let total_bp: usize = contigs.iter().map(|(_, s)| s.len()).sum();
    eprintln!("contigs={} bp={}", contigs.len(), total_bp);

    let e = encode(&contigs);
    let nsuf = e.num_suffixes();
    eprintln!(
        "total_orig={} total_logical={} midpos={} num_contigs={} k={} nsuffixes={}",
        e.total_orig, e.total_logical, e.midpos, e.num_contigs, e.k, nsuf
    );

    let (sa, lcp) = sa_lcp(&e.sa_input, e.k);
    // ACGT suffixes sort first (codes < sentinels): take the nsuf-length prefix.
    debug_assert!(e.sa_input[sa[nsuf - 1] as usize] < ALPHA as i32);
    debug_assert!(nsuf == sa.len() || e.sa_input[sa[nsuf] as usize] >= ALPHA as i32);
    let suftab: Vec<u64> = sa[..nsuf].iter().map(|&x| x as u64).collect();
    let lcptab: Vec<u64> = lcp[..nsuf].iter().map(|&x| x as u64).collect();

    let mut seeds: Vec<Seed> = Vec::new();
    enumerate_maxpairs(&suftab, &lcptab, params::SEED, ALPHA, &e.enc, |len, p1, p2| {
        store_seed(
            &mut seeds, len, p1, p2, e.midpos, e.total_logical, e.num_contigs,
            &e.seqnum_of, params::MIN_TIR_DIST, params::MAX_TIR_DIST, params::MAX_TIR_LEN,
        );
    });

    println!("SEEDS\t{}", seeds.len());
    if dump {
        let mut s = seeds.clone();
        s.sort_by_key(|x| (x.pos1, x.pos2, x.len));
        for sd in &s {
            println!(
                "SEED\tpos1={}\tpos2={}\tlen={}\tdist={}\tseqnum1={}",
                sd.pos1, sd.pos2, sd.len, sd.distance, sd.contignumber
            );
        }
    }
}
