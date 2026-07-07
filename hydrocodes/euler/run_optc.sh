#!/bin/bash
# Option C calibration runner (PICKUP-optionC.md 3.5): one implicit-Bingham crater run.
#   usage: ./run_optc.sh <a_m> <eta_Pa_s> <Y_B_Pa> [TDEC_s] [tend_s] [trace]
# Fixed: U=12 km/s, Mars g, cppr5, Rfac30/Zfac22, ROCK on, Y_d0=5 MPa, VOID_AMB=0.27, VISC_IMP=1.
# TDEC default = 2x the crater gravitational timescale 2pi*sqrt(D/g) at the anchor
# (a=300->550s, 1000->900s, 2000->1200s, 3000->1400s); tend default = 4x TDEC so decay
# outlasts collapse (the TDEC-insensitivity acceptance needs arrest by stress balance).
# Detach long runs with:  nohup caffeinate -i -m -s ./run_optc.sh ... &   (2026-07 gotcha)
set -e
cd "$(dirname "$0")"
a=$1; eta=$2; yb=$3
case $a in 300) td=550;; 1000) td=900;; 2000) td=1200;; 3000) td=1400;; *) td=900;; esac
tdec=${4:-$td}
tend=${5:-$(( 4 * tdec ))}
tag="a${a}_e${eta}_yb${yb}_td${tdec}"
mkdir -p state/optc
if [ "${6:-}" = "trace" ]; then
  mkdir -p "state/optc/trace_${tag}"
  export CRATER_TRACE=1 CRATER_SNAP_DT=200 CRATER_SNAP_DIR="state/optc/trace_${tag}"
fi
echo "[optc] $tag  TDEC=$tdec tend=$tend  trace=${6:-no}"
exec ./hydro_cpu crater "$a" 12000 "$tdec" "$eta" 3.71 5 "$tend" 30 22 \
  "state/optc/${tag}_prof.txt" 1 5e6 0.27 "$yb" 1 \
  > "state/optc/${tag}.out" 2> "state/optc/${tag}.log"
