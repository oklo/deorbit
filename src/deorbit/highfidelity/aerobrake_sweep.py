#!/usr/bin/env python3
"""Survival-map sweep: for a grid of capture ICs (periapsis -> first apogee,
Moon phase / epoch), run the full-fidelity aerobraking sequence and classify
the outcome. The key result is the boundary in FIRST APOGEE above which the
Moon destabilizes the capture (steep impact or ejection) vs. below which the
body gently aerobrakes to an ULTRA-OBLIQUE GRAZING impact.

Regimes (see aerobrake.classify):
  grazing_impact  -- gentle aerobrake to |fpa| < 5 deg  (the target scenario)
  steep_impact    -- Moon pumps periapsis into a steep (|fpa| >= 5 deg) impact
  ejected         -- lost beyond Earth's Hill sphere / unbound by perturbations
  lunar_encounter -- close approach to the Moon
  aerobraking     -- unresolved within the pass/walltime budget

Local daemon (free, heavy):  ../.venv/bin/python aerobrake_sweep.py --daemon
Bounded step:                ../.venv/bin/python aerobrake_sweep.py --increment 1200
"""
import os, sys, json, time, math, argparse
from . import aerobrake

# IC grid (periapsis altitude in km controls the first apogee)
PERI_KM = [6, 8, 10, 12, 14, 16, 18, 20, 24, 28]
EPOCH_D = [0.0, 6.75, 13.5, 20.25]          # ~quarter-phases over a lunar month
VINF_KMS = 0.2
INC_DEG = 5.0
PASS_CAP = 220
WALLTIME_PER_IC = 45.0


def grid():
    # epoch outer, periapsis inner: sweeps a full first-apogee transition curve
    # at one Moon phase first, then repeats at other phases.
    return [(p, e) for e in EPOCH_D for p in PERI_KM]


def run():
    ap = argparse.ArgumentParser()
    ap.add_argument("--daemon", action="store_true")
    ap.add_argument("--increment", type=float, default=None, help="seconds")
    ap.add_argument("--statedir", default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "state_sweep"))
    a = ap.parse_args()
    sd = a.statedir
    os.makedirs(sd, exist_ok=True)
    resf = os.path.join(sd, "survival.json")
    statusf = os.path.join(sd, "status.txt")
    statusj = os.path.join(sd, "status.json")
    pidf = os.path.join(sd, "daemon.pid")
    scratch = os.path.join(sd, "scratch")

    st = {"cell": 0, "results": [], "started": time.time()}
    if os.path.exists(os.path.join(sd, "state.json")):
        try:
            st.update(json.load(open(os.path.join(sd, "state.json"))))
        except Exception:
            pass
    with open(pidf, "w") as fh:
        fh.write(str(os.getpid()))

    cells = grid()

    def summarize():
        rs = st["results"]
        from collections import Counter
        c = Counter(r["regime"] for r in rs)
        graz = [r for r in rs if r["regime"] == "grazing_impact"]
        bad = [r for r in rs if r["regime"] in ("steep_impact", "ejected")]
        max_graz_apo = max((r["first_apo_km"] for r in graz), default=0.0)
        min_bad_apo = min((r["first_apo_km"] for r in bad), default=0.0)
        return c, max_graz_apo, min_bad_apo

    def save():
        json.dump(st, open(os.path.join(sd, "state.json") + ".tmp", "w"))
        os.replace(os.path.join(sd, "state.json") + ".tmp", os.path.join(sd, "state.json"))
        json.dump({"results": st["results"]}, open(resf + ".tmp", "w"), indent=1)
        os.replace(resf + ".tmp", resf)
        c, mga, mba = summarize()
        n = len(st["results"])
        with open(statusf, "w") as fh:
            fh.write("aerobrake survival map (full-fidelity REBOUND+ASSIST+drag)\n"
                     f"v_inf={VINF_KMS} km/s, inc={INC_DEG} deg; grid {len(cells)} cells\n"
                     f"done: {n}/{len(cells)}\n"
                     f"regimes: {dict(c)}\n"
                     f"max first-apogee that still GRAZES: {mga:,.0f} km\n"
                     f"min first-apogee that is LOST (steep/eject): {mba:,.0f} km\n"
                     f"(Moon orbit = 384,400 km)\n"
                     f"updated {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        c2, mga, mba = summarize()
        json.dump({"name": "deorbit — lunar survival map",
                   "goal": "first-apogee boundary: grazing aerobrake vs Moon-destabilized",
                   "metric": "max_grazing_first_apogee_km", "value": round(mga),
                   "detail": f"{dict(c2)}; lost above ~{mba:,.0f} km",
                   "iterations": n, "solved": False, "updated": time.time()},
                  open(statusj, "w"), indent=1)

    deadline = None if a.daemon else time.time() + (a.increment or 1800)
    while True:
        if st["cell"] >= len(cells):
            save()
            if not a.daemon:
                return
            time.sleep(30)
            continue
        peri, epoch = cells[st["cell"]]
        try:
            outcome, rows, meta = aerobrake.run(
                VINF_KMS, peri, INC_DEG, PASS_CAP, scratch,
                epoch_days=epoch, max_walltime=WALLTIME_PER_IC)
            regime, diag = aerobrake.classify(outcome, rows)
            full = [r for r in rows if "apo_alt_km" in r]
            first_apo = full[0]["apo_alt_km"] if full else 0.0
            rec = {"peri_km": peri, "epoch_d": epoch, "first_apo_km": first_apo,
                   "regime": regime, "n_passes": len(full),
                   "min_moon_km": min((r["moon_dist_km"] for r in full), default=None)}
            if diag and diag.get("event") == "impact":
                rec.update(fpa_deg=diag["fpa_deg"], speed_kms=diag["speed_kms"],
                           lat_deg=diag["lat_deg"], t_days=diag["t_days"])
            st["results"].append(rec)
        except Exception as e:
            st["results"].append({"peri_km": peri, "epoch_d": epoch,
                                   "regime": "error", "err": str(e)})
        st["cell"] += 1
        save()
        if deadline is not None and time.time() >= deadline:
            return


if __name__ == "__main__":
    run()
