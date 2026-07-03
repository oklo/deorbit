#!/usr/bin/env python3
"""Phase 3 AF calibration sweep (free local-compute grind).

Runs the `crater` driver over impactor sizes x dimensionless AF parameters,
dumping each final surface profile for ROBUST offline depth-diameter measurement
and appending a summary row to state/results.csv. Idempotent/resumable: a run
already in results.csv is skipped, so re-invoking continues where it left off.

The AF block model is parameterised in the literature (Wuennemann & Ivanov 2003)
in dimensionless, impactor-scaled form so ONE pair (Tfrac, Efrac) applies to all
sizes -- that is what lets a single AF calibration reproduce the whole
depth-diameter curve. We map:
    TDEC   = Tfrac * a / c_s          (decay time in impactor acoustic units)
    ETA_AF = Efrac * rho * c_s * a    (viscosity, impactor-scaled)
with basalt c_s~3000 m/s, rho=2700. So the daemon sweeps (Tfrac, Efrac, a) and
computes the dimensional driver args.

Usage:
    sweep_crater.py --increment N     # run up to N pending jobs, then exit (cloud routine)
    sweep_crater.py --daemon          # loop forever, grinding pending jobs (local supervisor)
    sweep_crater.py --list            # show pending/done counts
"""
import csv, os, subprocess, sys, time, itertools

HERE   = os.path.dirname(os.path.abspath(__file__))
STATE  = os.path.join(HERE, "state")
PROF   = os.path.join(STATE, "profiles")
RESCSV = os.path.join(STATE, "results.csv")
# Binary: CPU oracle by default; --gpu (or --bin PATH) selects the GPU crater driver.
# NOTE: the GPU is a single shared device -> use --jobs 1 with --gpu (concurrent GPU procs contend).
BIN    = os.path.join(HERE, "hydro_gpu" if "--gpu" in sys.argv else "hydro_cpu")
if "--bin" in sys.argv:
    BIN = sys.argv[sys.argv.index("--bin") + 1]

# --- physical constants / fixed run params ---
CS, RHO = 3000.0, 2700.0     # basalt sound speed, density (AF scaling)
U, G, CPPR = 12000.0, 3.71, 5.0   # impact speed (m/s), Mars gravity, cells per impactor radius (5 = first-pass curve; converge cppr 8/12 on key sizes later)
VAMB = 0.27   # void reset density (crater arg14): AMBIENT, not the lithostatic reference. The legacy reference-reset REFILLED the collapsing floor with rock (2026-07-03 A/B: a=3000 d/D 0.018->0.165) -- results before that date were tuned on top of the artifact (archived state/results_legacy_vamb0.csv).

# --- sweep grid ---
SIZES = [300, 500, 750, 1000, 1500, 2000, 3000]   # impactor radii (m) -> craters spanning the simple/complex transition
# (label, Tfrac, Efrac, Yd0_Pa). First pass: vary the damaged-cohesion Y_d0 (the creep knob) with AF off (the
# strength baseline) and with a fixed reasonable AF-on setting (Tfrac=30, Efrac=0.01). ROCK friction always on.
SETTINGS = []
for yd in (1.0e6, 5.0e6, 1.0e7):
    tag = f"{int(yd/1e6)}M"
    SETTINGS.append((f"off_yd{tag}", 0.0,  0.0,  yd))    # AF off, friction only
    SETTINGS.append((f"af_yd{tag}",  30.0, 0.01, yd))    # AF on (1x viscosity), + friction

FIELDS = ["af_label", "a_m", "Tfrac", "Efrac", "Yd0_Pa", "TDEC_s", "ETA_Pas",
          "D_app_m", "depth_m", "transient_m", "dD", "profile", "vamb"]

def jobs():
    for (label, Tfrac, Efrac, Yd0), a in itertools.product(SETTINGS, SIZES):
        TDEC = Tfrac * a / CS if Tfrac > 0 else 0.0
        ETA  = Efrac * RHO * CS * a if Efrac > 0 else 0.0
        yield label, a, Tfrac, Efrac, Yd0, TDEC, ETA

def done_keys():
    if not os.path.exists(RESCSV):
        return set()
    with open(RESCSV) as f:
        return {(r["af_label"], r["a_m"]) for r in csv.DictReader(f)}

def append_row(row):
    new = not os.path.exists(RESCSV)
    with open(RESCSV, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        if new:
            w.writeheader()
        w.writerow(row)

def job_key(j):
    return (j[0], str(int(j[1])))   # (af_label, a) -- unique per run

def launch(job):
    label, a, Tfrac, Efrac, Yd0, TDEC, ETA = job
    os.makedirs(PROF, exist_ok=True)
    prof = os.path.join(PROF, f"{label}_a{int(a)}.txt")
    # crater a U TDEC ETA g cppr tend Rfac Zfac profile rock Yd0 vamb  (tend auto via -1; rock=1)
    args = [BIN, "crater", str(a), str(U), f"{TDEC:.6g}", f"{ETA:.6g}",
            str(G), str(CPPR), "-1", "30", "22", prof, "1", f"{Yd0:.6g}", f"{VAMB:g}"]
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    return dict(p=p, job=job, prof=prof, t0=time.time())

def reap(h):
    label, a, Tfrac, Efrac, Yd0, TDEC, ETA = h["job"]
    out, _ = h["p"].communicate()
    res = next((l for l in out.splitlines() if l.startswith("RESULT")), None)
    if not res:
        print(f"  !! no RESULT for {label} a={a} (run failed?)", flush=True); return
    _, _a, D, depth, trans, dD = res.split()
    append_row(dict(af_label=label, a_m=int(a), Tfrac=Tfrac, Efrac=Efrac, Yd0_Pa=f"{Yd0:.4g}",
                    TDEC_s=f"{TDEC:.4g}", ETA_Pas=f"{ETA:.4g}", D_app_m=D,
                    depth_m=depth, transient_m=trans, dD=dD, profile=os.path.relpath(h["prof"], HERE), vamb=VAMB))
    print(f"  done {label} a={int(a)}m: depth={float(depth)/1e3:.2f}km transient={float(trans)/1e3:.2f}km ({time.time()-h['t0']:.0f}s)", flush=True)

def pending():
    dk = done_keys()
    return [j for j in jobs() if job_key(j) not in dk]

def main():
    os.makedirs(STATE, exist_ok=True)
    if "--list" in sys.argv:
        p = pending(); tot = len(SETTINGS) * len(SIZES)
        print(f"jobs: {tot} total, {tot-len(p)} done, {len(p)} pending"); return
    K = int(sys.argv[sys.argv.index("--jobs") + 1]) if "--jobs" in sys.argv else 1   # max concurrent runs
    inc = int(sys.argv[sys.argv.index("--increment") + 1]) if "--increment" in sys.argv else None
    daemon = "--daemon" in sys.argv
    running, inflight, done = [], set(), 0          # handles, in-flight keys, completed count
    while True:
        # fill the pool up to K (skip done + in-flight); respect the increment budget
        dk = done_keys() | inflight
        for j in jobs():
            if len(running) >= K: break
            if inc is not None and done + len(running) >= inc: break
            if job_key(j) in dk: continue
            running.append(launch(j)); inflight.add(job_key(j)); dk.add(job_key(j))
        if not running:
            if daemon:
                print("all jobs done; sleeping 1h", flush=True); time.sleep(3600); continue
            print("all jobs done."); break
        time.sleep(2)
        still = []
        for h in running:
            if h["p"].poll() is None: still.append(h); continue
            reap(h); inflight.discard(job_key(h["job"])); done += 1
        running = still
        if inc is not None and done >= inc:
            print(f"increment {inc} complete ({len(pending())} pending)", flush=True); break

if __name__ == "__main__":
    main()
