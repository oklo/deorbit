#!/bin/bash
# Overnight 1km / 41M Orcus pipeline: relax substrate -> build impact IC -> impact.
# Same physics as the 3km run (40km, 10 km/s, 5deg); only the resolution changes.
set -e
cd /Users/greglaughlin/Projects/deorbit/gpu
LOG=orcus_1km.log
{
  echo "=== $(date) START 1km/41M pipeline ==="
  echo "[1/3 relax] slab -> t=150 s, damp=0.99, g=3.71 (~2.3 h)"
  rm -f gpu_snap_*.bin
  ./gpu_ic orcus_slab_1km.bin 1000 150 14400 --frac --g 3.71 --damp 0.99 --snap 50
  LAST=$(ls gpu_snap_*.bin | tail -1); cp "$LAST" /tmp/orcus_relaxed_1km_snap.bin
  echo "[2/3 build] relaxed impact IC from $LAST"
  cd ../sph
  ../.venv/bin/python make_orcus_relaxed_ic.py --snap /tmp/orcus_relaxed_1km_snap.bin --dx 1000 --out orcus_relaxed_1km_ic.bin
  cd ../gpu
  echo "[3/3 impact] relaxed IC -> t=200 s, snap 5 (~11.5 h)"
  rm -f gpu_snap_*.bin
  ./gpu_ic orcus_relaxed_1km_ic.bin 1000 200 54000 --frac --g 3.71 --snap 5
  echo "=== $(date) DONE 1km/41M pipeline ==="
} >> $LOG 2>&1
