//! Stage 4 validation: stages 1-4 per seed; for each pair that reaches the
//! similarity computation (gt's PAIR-trace set), dump ulen/vlen/edist/sim.
//! Compare against gt TIRVISH_TRACE PAIR lines. Usage: simcheck <fasta>

use tirvish_rs::encode::{encode, ALPHA};
use tirvish_rs::maxpairs::enumerate_maxpairs;
use tirvish_rs::params;
use tirvish_rs::sa::sa_lcp;
use tirvish_rs::seeds::{store_seed, Seed};
use tirvish_rs::similarity::compute_similarity;
use tirvish_rs::tsd::{build_pair, search_for_tsds};
use tirvish_rs::xdrop::{calc_distances, extend_seed, ArbitraryScores};


fn main() {
    let path = std::env::args().nth(1).expect("usage: simcheck <fasta>");
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

    let mut out = String::new();
    let mut reached = 0usize;
    for s in &seeds {
        let (s1, e1, s2, e2) = e.contig_bounds(s.contignumber);
        let alilen = params::MAX_TIR_LEN - s.len;
        let (xl, xr) = extend_seed(
            &e.enc, s.pos1, s.pos2, s.len, s1, e1, s2, e2, alilen, &scores, &dist,
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
        if !pair.skip
            && (pair.left_tir_end <= pair.left_tir_start
                || pair.right_tir_end <= pair.right_tir_start)
        {
            pair.skip = true;
        }
        if pair.skip {
            continue;
        }
        // reaches the similarity computation (gt dumps PAIR here)
        let (ulen, vlen, edist) = compute_similarity(&mut pair, &e.twobit, params::SIMILARITY_THRESHOLD);
        reached += 1;
        out.push_str(&format!(
            "P\tp1={}\tp2={}\tslen={}\tulen={}\tvlen={}\tedist={}\tsim={:.2}\n",
            s.pos1, s.pos2, s.len, ulen, vlen, edist, pair.similarity
        ));
    }
    print!("{out}");
    eprintln!("pairs reaching similarity: {reached}");
}
