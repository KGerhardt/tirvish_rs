//! Stage-2 validation: for every seed, dump (pos1,pos2,len, Li,Lj,Ri,Rj) where
//! L*/R* are the left/right Xdrop ivalue/jvalue. Compare against gt's TIRVISH_TRACE
//! PAIR lines (p1,p2,slen,Li,Lj,Ri,Rj). Usage: xdropcheck <fasta>

use tirvish_rs::encode::{encode, ALPHA};
use tirvish_rs::maxpairs::enumerate_maxpairs;
use tirvish_rs::params;
use tirvish_rs::sa::sa_lcp;
use tirvish_rs::seeds::{store_seed, Seed};
use tirvish_rs::xdrop::{calc_distances, extend_seed, ArbitraryScores};


fn main() {
    let path = std::env::args().nth(1).expect("usage: xdropcheck <fasta>");
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
        mat: params::XDROP_MAT,
        mis: params::XDROP_MIS,
        ins: params::XDROP_INS,
        del: params::XDROP_DEL,
    };
    let dist = calc_distances(&scores);
    let belowscore = params::XDROP_BELOWSCORE;

    let mut out = String::new();
    for s in &seeds {
        let (s1, e1, s2, e2) = e.contig_bounds(s.contignumber);
        let alilen = params::MAX_TIR_LEN - s.len; // gt: max_tir_length - seedptr->len
        let (xl, xr) = extend_seed(
            &e.twobit, s.pos1, s.pos2, s.len, s1, e1, s2, e2, alilen, &scores, &dist, belowscore,
        );
        out.push_str(&format!(
            "XD\tp1={}\tp2={}\tslen={}\tLi={}\tLj={}\tRi={}\tRj={}\n",
            s.pos1, s.pos2, s.len, xl.ivalue, xl.jvalue, xr.ivalue, xr.jvalue
        ));
    }
    print!("{out}");
    eprintln!("xdrop computed for {} seeds", seeds.len());
}
