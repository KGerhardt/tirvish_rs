#!/usr/bin/env python3
"""
Translate a `gt tirvish` GFF3 into the exact tabular shape `tirvish_rs` emits, so
the two can be diffed directly.

The ingestion is a verbatim template of TIR-Learner's one_tirvish GFF parse: every
non-'#' line contributes a (start, stop) pair (GFF columns 4 and 5) to the current
record in order; seqid is column 1; any '#' line flushes a non-empty record. A gt
tirvish element is six consecutive feature lines, read BY POSITION (no type check):

  next_result[0] repeat_region                    -> full element (incl. TSDs)
  next_result[1] target_site_duplication          -> TSD1
  next_result[2] terminal_inverted_repeat_element -> body (no TSDs)
  next_result[3] terminal_inverted_repeat         -> TIR1  (gt emits the two TIRs
  next_result[4] terminal_inverted_repeat         -> TIR2   sorted by (start,end))
  next_result[5] target_site_duplication          -> TSD2

Coordinates are emitted verbatim (1-based inclusive) -- nothing is anchored or
recomputed; tirvish_rs computes the identical coordinates internally. `sim` is gt's
tir_similarity (on the *_element line), formatted to 2 decimals to match
tirvish_rs. No TA%/N%/WFA filtering: this is the raw prediction multiset.

Output columns (identical order to tirvish_rs):
  seqid full_start full_stop tsd1_start tsd1_stop body_start body_stop
  tir1_start tir1_stop tir2_start tir2_stop tsd2_start tsd2_stop sim

Usage: python3 gff_to_tsv.py CHUNK.gff [out.tsv]   (default out: stdout)
"""
import sys
import re

sim_regex = re.compile(r'tir_similarity=([\d.]+)')

HEADER = ("seqid\tfull_start\tfull_stop\ttsd1_start\ttsd1_stop\tbody_start\tbody_stop\t"
          "tir1_start\ttir1_stop\ttir2_start\ttir2_stop\ttsd2_start\ttsd2_stop\tsim")


def parse(path):
    rows = []
    bad = []
    next_result = []
    seqid = None
    sim = ""

    def flush():
        if not next_result:
            return
        if len(next_result) != 6:
            bad.append((seqid, len(next_result)))
            return
        full = next_result[0]; tsd1 = next_result[1]; body = next_result[2]
        tir1 = next_result[3]; tir2 = next_result[4]; tsd2 = next_result[5]
        sim_s = f"{float(sim):.2f}" if sim else "0.00"
        rows.append((seqid,
                     full[0], full[1], tsd1[0], tsd1[1], body[0], body[1],
                     tir1[0], tir1[1], tir2[0], tir2[1], tsd2[0], tsd2[1], sim_s))

    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                flush()
                next_result = []
                sim = ""
                continue
            segs = line.rstrip('\n').split('\t')
            if len(segs) < 5:
                continue
            seqid = segs[0]
            next_result.append((int(segs[3]), int(segs[4])))
            if len(segs) >= 9:
                mm = sim_regex.search(segs[8])
                if mm:
                    sim = mm.group(1)
    flush()
    return rows, bad


def main():
    gff = sys.argv[1]
    rows, bad = parse(gff)
    out = open(sys.argv[2], 'w') if len(sys.argv) > 2 else sys.stdout
    out.write(HEADER + "\n")
    for r in rows:
        out.write("\t".join(map(str, r)) + "\n")
    if out is not sys.stdout:
        out.close()
    if bad:
        sys.stderr.write(f"WARNING: {len(bad)} records with !=6 features skipped: "
                         f"{bad[:5]}{'...' if len(bad) > 5 else ''}\n")


if __name__ == "__main__":
    main()
