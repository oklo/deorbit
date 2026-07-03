#!/bin/bash
# Pierazzo-2008 45-deg oblique Al benchmark ladder (Table 1 lower block; oracle = pierazzo2008_oracle.md).
# Runs cheap-first so partial results land early; each rung logs a PZ45SUMMARY line + keeps its decay files.
# Launch DETACHED (survives session close + sleep):
#   cd state/pz45 && nohup caffeinate -i -m -s ../../run_pierazzo45_ladder.sh > ladder45.log 2>&1 &
set -u
cd "$(dirname "$0")"
BIN=$(pwd)/hydro_gpu
RCFL=${RCFL:-100}   # 100 = prater-validated free-surface treatment (voidzero-to-ambient + vacuum flux): ~200x fewer steps, and 20 km/s cppr>=16 NaNs without it. RCFL=0 env = legacy ambient-live (kept for cross-checks).
OUT=state/pz45_rcfl$RCFL
mkdir -p "$OUT"; cd "$OUT"
REACH=12; HALF=6; HUP=6; HDN=12
echo "=== pz45 ladder START $(date '+%F %T') (reach=$REACH half=$HALF hup=$HUP hdn=$HDN rcfl=$RCFL) ==="
for rung in "12 5000" "12 20000" "16 5000" "16 20000" "20 5000" "20 20000"; do
  set -- $rung; cppr=$1; U=$2
  tag="u${U}_c${cppr}"
  echo "--- $tag START $(date '+%F %T') ---"
  "$BIN" pierazzo "$cppr" $REACH $HALF "$U" -1 $RCFL 45 $HUP $HDN > "out_$tag.txt" 2> "err_$tag.txt"
  for f in down ray row; do mv -f "pierazzo45_$f.txt" "decay45_${f}_$tag.txt" 2>/dev/null; done
  grep -E "PZ45SUMMARY|GATE|paper" "out_$tag.txt"
  echo "--- $tag DONE $(date '+%F %T') ---"
done
echo "=== pz45 ladder DONE $(date '+%F %T') ==="
