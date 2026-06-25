#!/bin/bash
# Self-chaining driver for the fracture run: waits for the currently-running
# 12 h cap to finish, then keeps restarting from the lossless checkpoint
# (frac_chk.bin -> preserves damage D, stress S, flaws) until t_end is reached.
# Fully detached (nohup) so it survives the controlling session ending.
# Needs only: laptop awake (lid open + plugged); each cap is caffeinate-wrapped.
cd /Users/greglaughlin/Projects/deorbit/sph || exit 1
TEND=1.2          # target sim time (raise to go longer)
CAP=43200         # walltime per cap (s) = 12 h
MAXR=6            # safety cap on restart count (0.6->1.2 needs ~1-2)
log(){ echo "[$(date '+%F %T')] $*" >> frac_driver.log; }

log "driver up; t_end=$TEND, cap=${CAP}s. waiting for the current cap to finish..."
while pgrep -f "run_ic2 --frac cayambe_ic.bin" >/dev/null; do sleep 120; done
log "current cap finished; entering restart chain"

for i in $(seq 1 $MAXR); do
    reason=$(grep "done (" frac.log | tail -1)
    if echo "$reason" | grep -q "t_end"; then log "REACHED t_end ($reason) -- stopping"; break; fi
    if [ ! -f frac_chk.bin ]; then log "no frac_chk.bin -- aborting"; break; fi
    log "restart #$i from frac_chk.bin -> t_end=$TEND  (previous: $reason)"
    caffeinate -dims ./run_ic2 --frac --restart frac_chk.bin "$TEND" "$CAP" >> frac.log 2>&1
    log "restart #$i exited"
done
log "driver done"
