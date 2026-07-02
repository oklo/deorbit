#!/bin/bash
# Pierazzo-2008 Al peak-pressure benchmark: 3D scale-up run (attacks the M5 caveat that the
# near-field n~1.3 is unconverged and the far-field ~2 is unseen because the old grid reached
# only r/a~7). Two parts: a RESOLUTION-convergence ladder at reach(r/a)=20 (cppr 8,10,12) and
# one DEEP-domain run (reach=30) to pin the asymptotic far-field slope.
# NOTE: pierazzo's CFL is limited by live low-density ambient cells (rcfl floor destabilises it),
# so 3D cost scales ~cppr^5-6; hence the modest top cppr. Runs FOREGROUND (normal QoS).
# Usage: caffeinate -i -s ./run_pierazzo_ladder.sh
set -u
cd "$(dirname "$0")"
BIN="$PWD/hydro_gpu"
OUT="$PWD/state/pierazzo"; mkdir -p "$OUT"
LOG="$OUT/ladder.log"
U=10000                 # Pierazzo Al benchmark impact velocity (m/s)
DEADLINE=$((16*3600))   # do not START a new rung after 16 h elapsed
START=$(date +%s)

# label  cppr  reach  half
RUNS=(
  "conv_cppr8   8  20 10"
  "conv_cppr10 10  20 10"
  "conv_cppr12 12  20 10"
  "deep_r30    10  30 12"
)

echo "=== Pierazzo 3D scale-up START $(date '+%F %T') | U=$U deadline=16h ===" | tee -a "$LOG"
for spec in "${RUNS[@]}"; do
  set -- $spec; label=$1; c=$2; reach=$3; half=$4
  elapsed=$(( $(date +%s) - START ))
  if [ "$elapsed" -gt "$DEADLINE" ]; then
    echo "!!! deadline reached (${elapsed}s) -- skipping $label and the rest" | tee -a "$LOG"; break
  fi
  t0=$(date +%s)
  echo "--- $label (cppr=$c reach=$reach half=$half) START $(date '+%F %T') elapsed=$((elapsed/60))min ---" | tee -a "$LOG"
  ( cd "$OUT" && "$BIN" pierazzo "$c" "$reach" "$half" "$U" 2>"rung_${label}.stderr" ) > "$OUT/rung_${label}.out"
  [ -f "$OUT/pierazzo_decay.txt" ] && mv "$OUT/pierazzo_decay.txt" "$OUT/decay_${label}.txt"
  grep -hE "PZSUMMARY|Pierazzo Al-sphere|GATE" "$OUT/rung_${label}.out" | sed "s/^/  [$label] /" | tee -a "$LOG"
  t1=$(date +%s); echo "--- $label DONE $(date '+%F %T') wall=$(( (t1-t0)/60 )) min ---" | tee -a "$LOG"
done
echo "=== Pierazzo 3D scale-up DONE $(date '+%F %T') total=$(( ($(date +%s)-START)/60 )) min ===" | tee -a "$LOG"
