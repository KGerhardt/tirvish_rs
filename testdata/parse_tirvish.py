#!/usr/bin/env python3
"""
Parse a gt tirvish GFF3 file into the gold candidate set.

Field extraction mirrors TIR-Learner's tirvish_new.one_tirvish EXACTLY (same
per-record feature order), but applies NONE of its TA%/N%/WFA filtering: the
gold set is the RAW gt tirvish prediction multiset, since that is the spec our
Rust port must reproduce. TIR-Learner's downstream filters run identically on
top of either seed source and so are out of scope for the faithfulness test.

Per-record GFF feature order (gt tirvish, 6 lines between `#` separators):
  [0] repeat_region              -> full element (with TSDs)
  [1] target_site_duplication    -> TSD1
  [2] terminal_inverted_repeat_element -> element without TSDs
  [3] terminal_inverted_repeat   -> TIR1 (left)
  [4] terminal_inverted_repeat   -> TIR2 (right)
  [5] target_site_duplication    -> TSD2

Output TSV (one row per prediction):
  seqid  start  stop  tir1  tir2  tsd1  tsd2  sim
The first 7 are the comparison KEY (start/stop = full element incl. TSDs, 1-based
as gt emits; tir*/tsd* = lengths). `sim` = gt's tir_similarity attribute, a
diagnostic (NOT part of the key) -- the port computes this for the -similar gate,
so it's handy for debugging where similarity values diverge.

Usage: python3 parse_tirvish.py gold/chunkN.tirvish.gff [out.tsv]
"""
import sys, re

genome_split_regex = re.compile(r'(.+);;(\d+)')
sim_regex = re.compile(r'tir_similarity=([\d.]+)')


def parse(path):
    rows = []
    odd = 0  # records whose feature count != 6 (TSD-less or malformed)
    next_result = []
    sim = ""
    seqid = None

    def flush():
        nonlocal odd
        if not next_result:
            return
        if len(next_result) != 6:
            odd += 1
            return
        m = genome_split_regex.match(seqid)
        short_id = m.group(1) if m else seqid
        full_start, full_stop = min(next_result[0]), max(next_result[0])
        tsd1 = next_result[1][1] - next_result[1][0] + 1
        tsd2 = next_result[5][1] - next_result[5][0] + 1
        tir1 = next_result[3][1] - next_result[3][0] + 1
        tir2 = next_result[4][1] - next_result[4][0] + 1
        rows.append((short_id, full_start, full_stop, tir1, tir2, tsd1, tsd2, sim))

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
            start, end = int(segs[3]), int(segs[4])
            next_result.append((start, end))
            if len(segs) >= 9:  # tir_similarity lives on the *_element line
                mm = sim_regex.search(segs[8])
                if mm:
                    sim = mm.group(1)
    flush()  # final record (no trailing # may follow)
    return rows, odd


def main():
    gff = sys.argv[1]
    out = sys.argv[2] if len(sys.argv) > 2 else gff.replace('.tirvish.gff', '.gold.tsv')
    rows, odd = parse(gff)
    with open(out, 'w') as fh:
        fh.write("seqid\tstart\tstop\ttir1\ttir2\ttsd1\ttsd2\tsim\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    print(f"{gff}: {len(rows)} predictions -> {out}"
          + (f"  (WARNING: {odd} records with !=6 features, skipped)" if odd else ""))


if __name__ == "__main__":
    main()
