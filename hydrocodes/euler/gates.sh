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
CPU="sod shear yield vacuum tracer tensile freefall atmos alimpact af_activate substrate"
[ "$QUICK" = "quick" ] || CPU="sod sedov shear yield vacuum tracer tensile freefall atmos alimpact af_activate substrate"
echo "== CPU gates =="
for m in $CPU; do ./hydro_cpu "$m" 2>/dev/null | grep -i 'gate' | sed "s/^/  [$m] /"; done

# GPU gates (FP32). Should match the CPU oracle.
echo "== GPU gates =="
for m in sod sedov surface shear yield tensile freefall atmos pierazzo vacuum tracer af_activate substrate; do
  ./hydro_gpu "$m" 2>/dev/null | grep -i 'gate' | sed "s/^/  [$m] /"
done

echo "== done. (collapse = M6 AF-calibration demo, not a pass/fail gate yet) =="
