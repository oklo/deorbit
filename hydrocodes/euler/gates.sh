#!/bin/bash
# euler validation gates -- ordered build (CPU oracle, Metal lib, GPU host) + run every gate.
# Usage:  ./gates.sh           # build + all gates (CPU sedov is the slow 64^3 FP64 oracle, ~5 min)
#         ./gates.sh quick     # skip the slow CPU sedov (GPU sedov still runs, fast)
# Each gate prints PASS/CHECK; unknown modes are fatal (exit 2). Run from the euler/ directory.
set -e
cd "$(dirname "$0")"
QUICK=${1:-}

echo "== build =="
clang++ -std=c++17 -O2 hydro_cpu.cpp -o hydro_cpu
xcrun -sdk macosx metal -O3 -ffast-math -fmodules-cache-path=.clang-module-cache -c hydro.metal -o hydro.air
xcrun -sdk macosx metallib hydro.air -o hydro.metallib
clang++ -std=c++17 -O2 -I../common/metal-cpp hydro_gpu.cpp \
  -framework Metal -framework Foundation -framework QuartzCore -o hydro_gpu
echo "   ok"

# CPU oracle gates (FP64). sedov is 64^3 and slow -> skipped in 'quick'.
# visc_diff/visc_agree/visc_energy/bingham_slope = Option C implicit-Bingham gates (PICKUP-optionC.md 3.4 b-e;
# 3.4a feature-off bit-identity = all the pre-existing gates below still passing with VISC_IMP default-off).
CPU="sod shear yield vacuum tracer tensile freefall atmos alimpact pierazzo2d af_activate sedov_axi lame af_visc vib_advect friction substrate visc_diff visc_agree visc_energy bingham_slope"
[ "$QUICK" = "quick" ] || CPU="sod sedov shear yield vacuum tracer tensile freefall atmos alimpact pierazzo2d af_activate sedov_axi lame af_visc vib_advect friction substrate visc_diff visc_agree visc_energy bingham_slope"
echo "== CPU gates =="
for m in $CPU; do ./hydro_cpu "$m" 2>/dev/null | grep -i 'gate' | sed "s/^/  [$m] /"; done

# GPU gates (FP32). Should match the CPU oracle.
echo "== GPU gates =="
for m in sod sedov surface shear yield tensile freefall atmos pierazzo pierazzo2d vacuum tracer af_activate sedov_axi lame af_visc vib_advect friction substrate; do
  ./hydro_gpu "$m" 2>/dev/null | grep -i 'gate' | sed "s/^/  [$m] /"
done
# 45-deg oblique Pierazzo smoke gate (low-res, ~1 s; full oracle head-to-head = run_pierazzo45_ladder.sh at cppr>=16)
./hydro_gpu pierazzo 6 8 4 5000 -1 0 45 2>/dev/null | grep -i 'gate' | sed "s/^/  [pierazzo45] /"
# Prater-1970 Al-alloy strength validation (E23; oracle = prater1970_oracle.md). quick = cppr6 smoke (~20s+50s);
# full adds the cppr10 experiment-envelope gate on both codes (~1+2 min).
./hydro_cpu prater 6 2>/dev/null | grep -i 'gate' | sed "s/^/  [prater cpu6] /"
./hydro_gpu prater 6 2>/dev/null | grep -i 'gate' | sed "s/^/  [prater gpu6] /"
if [ "$QUICK" != "quick" ]; then
  ./hydro_cpu prater 10 2>/dev/null | grep -i 'gate' | sed "s/^/  [prater cpu10] /"
  ./hydro_gpu prater 10 2>/dev/null | grep -i 'gate' | sed "s/^/  [prater gpu10] /"
fi

echo "== done. (collapse = M6 AF-calibration demo, not a pass/fail gate yet) =="
