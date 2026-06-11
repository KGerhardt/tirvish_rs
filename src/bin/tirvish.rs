//! tirvish_rs — full de-novo TIRvish. Reads a FASTA, runs all 5 stages, emits the
//! gold-TSV shape (seqid start stop tir1 tir2 tsd1 tsd2 sim), sorted like the
//! oracle (seqid, start) so it diffs directly against gold/chunkN.gold.tsv.
//! Usage: tirvish <fasta>

use std::fs;
use tirvish_rs::pipeline::run;

fn read_fasta(path: &str) -> Vec<(String, Vec<u8>)> {
    let data = fs::read_to_string(path).expect("read fasta");
    let mut out = Vec::new();
    let (mut name, mut seq) = (String::new(), Vec::new());
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
    let path = std::env::args().nth(1).expect("usage: tirvish <fasta>");
    let contigs = read_fasta(&path);
    let mut els = run(&contigs);
    // sort like parse_tirvish gold output (by seqid as gt emits, then start)
    els.sort_by(|a, b| a.seqid.cmp(&b.seqid).then(a.start.cmp(&b.start)));
    println!("seqid\tstart\tstop\ttir1\ttir2\ttsd1\ttsd2\tsim");
    for el in &els {
        println!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}",
            el.seqid, el.start, el.stop, el.tir1, el.tir2, el.tsd1, el.tsd2, el.sim
        );
    }
    eprintln!("{} elements", els.len());
}
