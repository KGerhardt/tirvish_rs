//! Dump the (left-arm, right-arm) string pairs that stage 4 (similarity) compares,
//! one TAB-separated pair per line (ACGT; 'N' for specials), as a benchmark corpus
//! for evaluating Levenshtein implementations against our greedy_unit_edist.
//!
//! These are EXACTLY the two inputs to greedy_unit_edist for every pair that
//! reaches the similarity gate: useq = mirror[left_tir_start..left_tir_end],
//! vseq = mirror[right_tir_start..right_tir_end] (post-TSD, end-start lengths,
//! revcomp already baked into the mirror). Mirrors simcheck's pipeline so the set
//! is identical to the pairs gt would dump at the PAIR trace.
//!
//! Usage: simdump <fasta> > pairs.tsv   (writes a count to stderr)

use std::io::{BufWriter, Write};
use tirvish_rs::encode::{encode, ALPHA};
use tirvish_rs::maxpairs::enumerate_maxpairs;
use tirvish_rs::params;
use tirvish_rs::sa::sa_lcp;
use tirvish_rs::seeds::{store_seed, Seed};
use tirvish_rs::tsd::{build_pair, search_for_tsds};
use tirvish_rs::xdrop::{calc_distances, extend_seed, ArbitraryScores};

fn main() {
    let path = std::env::args().nth(1).expect("usage: simdump <fasta> > pairs.tsv");
    let contigs = tirvish_rs::read_fasta(&path);
    let e = encode(&contigs);
    let nsuf = e.num_suffixes();
    let (sa, lcp) = sa_lcp(&e.sa_input, e.k);
    let suftab: Vec<u64> = sa[..nsuf].iter().map(|&x| x as u64).collect();
    let lcptab: Vec<u64> = lcp[..nsuf].iter().map(|&x| x as u64).collect();

    let mut seeds: Vec<Seed> = Vec::new();
    enumerate_maxpairs(&suftab, &lcptab, params::SEED, ALPHA, &e.enc, |len, p1, p2| {
        store_seed(
            &mut seeds, len, p1, p2, e.midpos, e.total_logical, e.num_contigs,
            &e.seqnum_of, params::MIN_TIR_DIST, params::MAX_TIR_DIST, params::MAX_TIR_LEN,
        );
    });

    let scores = ArbitraryScores {
        mat: params::XDROP_MAT, mis: params::XDROP_MIS,
        ins: params::XDROP_INS, del: params::XDROP_DEL,
    };
    let dist = calc_distances(&scores);

    let stdout = std::io::stdout();
    let mut w = BufWriter::with_capacity(1 << 20, stdout.lock());
    let mut line: Vec<u8> = Vec::with_capacity(2048);
    let mut dumped = 0usize;
    for s in &seeds {
        let (s1, e1, s2, e2) = e.contig_bounds(s.contignumber);
        let alilen = params::MAX_TIR_LEN - s.len;
        let (xl, xr) = extend_seed(
            &e.twobit, s.pos1, s.pos2, s.len, s1, e1, s2, e2, alilen, &scores, &dist,
            params::XDROP_BELOWSCORE,
        );
        let mut pair = match build_pair(
            s.pos1, s.pos2, s.len, s.contignumber, xl.ivalue, xl.jvalue, xr.ivalue, xr.jvalue,
            e.total_logical, params::MIN_TIR_LEN, params::MAX_TIR_LEN,
        ) {
            Some(p) => p,
            None => continue,
        };
        let seq_start = e.fwd_seqstart[s.contignumber as usize];
        let seq_len = e.fwd_seqlen[s.contignumber as usize];
        search_for_tsds(
            &mut pair, &e.enc, seq_start, seq_len, params::VICINITY,
            params::MIN_TSD_LEN, params::MAX_TSD_LEN,
        );
        // gt's pre-similarity skip (degenerate arm) — same condition as the pipeline.
        if !pair.skip
            && (pair.left_tir_end <= pair.left_tir_start
                || pair.right_tir_end <= pair.right_tir_start)
        {
            pair.skip = true;
        }
        if pair.skip {
            continue;
        }
        // exactly greedy_unit_edist's two inputs: [start, end) over the mirror.
        line.clear();
        for p in pair.left_tir_start..pair.left_tir_end {
            line.push(e.twobit.base_at(p as usize));
        }
        line.push(b'\t');
        for p in pair.right_tir_start..pair.right_tir_end {
            line.push(e.twobit.base_at(p as usize));
        }
        line.push(b'\n');
        w.write_all(&line).unwrap();
        dumped += 1;
    }
    w.flush().unwrap();
    eprintln!("dumped {dumped} arm pairs");
}
