#!/usr/bin/env bash
# Run `gt tirvish` (original) and `tirvish_rs` side by side on the four shrimp test
# chunks, time each invocation with /usr/bin/time -v, and diff the predictions to
# demonstrate the outputs are identical.
#
# Requirements:
#   - genometools `gt` on PATH            (or set GT=/path/to/gt)
#   - a built tirvish_rs binary           (or set RS=/path/to/tirvish; else this
#                                           script cargo-builds it)
#   - python3, /usr/bin/time
#
# Usage:
#   ./run_compare.sh
#   GT=/opt/genometools/bin/gt RS=../target/release/tirvish ./run_compare.sh
#
# Exit status is non-zero if any chunk's predictions differ.

set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"

GT="${GT:-gt}"
RS="${RS:-$HERE/../target/release/tirvish}"
if [[ ! -x "$RS" ]]; then
  echo "building tirvish_rs (release)..."
  (cd "$HERE/.." && cargo build --release)
fi

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
tar -xzf chunks.tar.gz -C "$WORK"

# gt flags = exactly what TIR-Learner's tirvish_new.one_tirvish runs.
rc=0
for fa in "$WORK"/chunk*.fa; do
  base="$(basename "$fa" .fa)"
  echo "=================== $base ==================="
  idx="$WORK/$base.idx"

  # ---- original: gt suffixerator (mirrored ESA) + gt tirvish ----
  /usr/bin/time -v "$GT" suffixerator -db "$fa" -indexname "$idx" \
      -tis -suf -lcp -des -ssp -sds -dna -mirrored \
      >/dev/null 2> "$WORK/$base.gt_suf.time"
  /usr/bin/time -v "$GT" tirvish -index "$idx" \
      -seed 20 -mintirlen 10 -maxtirlen 1000 -mintirdist 10 -maxtirdist 5000 \
      -similar 80 -mintsd 2 -maxtsd 11 -vic 13 -seqids yes \
      > "$WORK/$base.gt.gff" 2> "$WORK/$base.gt_tirvish.time"

  # ---- tirvish_rs ----
  /usr/bin/time -v "$RS" "$fa" > "$WORK/$base.rs.tsv" 2> "$WORK/$base.rs.time"

  # ---- diff predictions (gt GFF -> TSV; compare on the element key, cols 1-7) ----
  python3 parse_tirvish.py "$WORK/$base.gt.gff" "$WORK/$base.gt.tsv" >/dev/null
  if diff -q \
       <(tail -n +2 "$WORK/$base.gt.tsv" | cut -f1-7 | sort) \
       <(tail -n +2 "$WORK/$base.rs.tsv" | cut -f1-7 | sort) >/dev/null; then
    n=$(( $(wc -l < "$WORK/$base.rs.tsv") - 1 ))
    echo "OUTPUTS IDENTICAL: $n predictions match"
  else
    echo "!!! MISMATCH !!!"
    diff <(tail -n +2 "$WORK/$base.gt.tsv" | cut -f1-7 | sort) \
         <(tail -n +2 "$WORK/$base.rs.tsv" | cut -f1-7 | sort) | head -20
    rc=1
  fi

  for stage in gt_suf gt_tirvish rs; do
    label=$([[ $stage == gt_suf ]] && echo "gt suffixerator" || { [[ $stage == gt_tirvish ]] && echo "gt tirvish" || echo "tirvish_rs"; })
    printf -- "--- %s ---\n" "$label"
    grep -E "Elapsed \(wall|Maximum resident" "$WORK/$base.$stage.time" || true
  done
done

echo
[[ $rc -eq 0 ]] && echo "ALL CHUNKS IDENTICAL." || echo "SOME CHUNKS DIFFERED."
exit $rc
