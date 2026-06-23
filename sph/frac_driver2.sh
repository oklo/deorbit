#!/bin/bash
# Seamless EXTENSION driver: keep the frac sim running to TARGET past t=1.2 with
# no idle CPU gap. The current run has t_end=1.2 baked in and will exit there;
# this driver tight-polls (15 s) so it restarts from the lossless checkpoint
# within ~15 s of that exit, and chains across 12 h walltime caps until TARGET.
# Holds its own blanket sleep-prevention. Detached (survives the session).
cd /Users/greglaughlin/Projects/deorbit/sph || exit 1
TARGET=2.0        # extend to here (blob exits the box ~1.7; curtain develops longer)
CAP=43200         # 12 h per cap
log(){ echo "[$(date '+%F %T')] $*" >> frac_driver2.log; }

nohup caffeinate -dims -w $$ >/dev/null 2>&1 &      # blanket: no sleep for this driver's life
log "driver2 up; TARGET=$TARGET (seamless extension past t=1.2)"

while true; do
    t=$(grep -oE 't=[0-9]+\.[0-9]+' frac.log | tail -1 | cut -d= -f2); [ -z "$t" ] && t=0
    if awk "BEGIN{exit !($t+0 >= $TARGET)}"; then log "reached TARGET (t=$t) -- stopping"; break; fi
    if pgrep -f 'run_ic2 --frac' >/dev/null; then sleep 15; continue; fi   # a run is active -> wait
    [ -f frac_chk.bin ] || { log "no frac_chk.bin -- abort"; break; }
    log "no run active at t=$t (< $TARGET) -> restart from checkpoint to TARGET"
    caffeinate -dims ./run_ic2 --frac --restart frac_chk.bin "$TARGET" "$CAP" >> frac.log 2>&1
    log "run segment exited"
done
log "driver2 done"
