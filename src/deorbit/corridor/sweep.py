#!/usr/bin/env python3
"""deorbit capture-corridor daemon (pure stdlib).

Sweeps the (v_inf, periapsis-altitude) plane for a fixed body (default: a 1 km
monolithic iron sphere) and classifies each grazing passage as impact / capture
/ skip via physics.simulate_passage. The deliverable is the CAPTURE CORRIDOR:
for which encounter speeds, and over how narrow a periapsis window, does the
body get bound rather than skip out or hit directly.

Operating model (see CLAUDE.md):
    python3 deorbit.py --daemon                 # free local 24/7, refines grid
    python3 deorbit.py --increment 1500         # bounded step (cloud routine), seconds

State (state/): checkpoints corridor.json + status.txt/status.json for the board.
Resolution grows over time (coarse -> fine) so longer runtime sharpens the corridor.
"""
import sys, os, json, time, argparse
from . import physics

V_INF_LO, V_INF_HI = 0.2e3, 5.0e3          # m/s
PERI_LO, PERI_HI = 0.0, 120e3              # m
MAX_LEVEL = 8


def _grid(level):
    n = 12 + 6 * level
    vs = [V_INF_LO + (V_INF_HI - V_INF_LO) * i / (n - 1) for i in range(n)]
    ps = [PERI_LO + (PERI_HI - PERI_LO) * j / (n - 1) for j in range(n)]
    return vs, ps


def _summary(results):
    cap = [r for r in results if r["outcome"] == "capture"]
    summ = {
        "n_cells": len(results),
        "n_capture": len(cap),
        "n_impact": sum(1 for r in results if r["outcome"] == "impact"),
        "n_skip": sum(1 for r in results if r["outcome"] == "skip"),
        "v_inf_max_capture_kms": (max(r["v_inf"] for r in cap) / 1e3) if cap else 0.0,
    }
    if cap:
        # widest periapsis capture window at any single v_inf
        by_v = {}
        for r in cap:
            by_v.setdefault(round(r["v_inf"], 1), []).append(r["peri_alt_km"])
        widths = {v: (max(p) - min(p)) for v, p in by_v.items()}
        best_v = max(widths, key=widths.get)
        summ["widest_window_km"] = widths[best_v]
        summ["widest_window_at_vinf_kms"] = best_v / 1e3
    else:
        summ["widest_window_km"] = 0.0
        summ["widest_window_at_vinf_kms"] = 0.0
    return summ


def run():
    ap = argparse.ArgumentParser()
    ap.add_argument("--daemon", action="store_true")
    ap.add_argument("--increment", type=float, default=None, help="seconds")
    ap.add_argument("--diameter", type=float, default=1000.0)
    ap.add_argument("--cd", type=float, default=1.0)
    ap.add_argument("--statedir", default=os.path.join(
        os.path.dirname(os.path.abspath(sys.argv[0])), "state"))
    a = ap.parse_args()
    sd = a.statedir
    os.makedirs(sd, exist_ok=True)
    corridorf = os.path.join(sd, "corridor.json")
    statef = os.path.join(sd, "state.json")
    statusf = os.path.join(sd, "status.txt")
    statusj = os.path.join(sd, "status.json")
    pidf = os.path.join(sd, "daemon.pid")
    logf = os.path.join(sd, "study.log")

    body = physics.Body(diameter_m=a.diameter, cd=a.cd)
    name = "deorbit — capture corridor"
    goal = ("km iron monolith: does low-v_inf grazing capture it (eps<0) or "
            "skip/impact? Map the corridor.")
    metric = "v_inf_max_capture_kms"

    st = {"level": 0, "cell": 0, "results": [], "started": time.time(),
          "diameter": a.diameter, "cd": a.cd, "beta": body.beta}
    if os.path.exists(statef):
        try:
            st.update(json.load(open(statef)))
        except Exception:
            pass
    with open(pidf, "w") as fh:
        fh.write(str(os.getpid()))

    def log(m):
        with open(logf, "a") as fh:
            fh.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {m}\n")

    def save():
        json.dump(st, open(statef + ".tmp", "w"))
        os.replace(statef + ".tmp", statef)
        summ = _summary(st["results"])
        json.dump({"summary": summ, "results": st["results"],
                   "body": {"diameter_m": a.diameter, "cd": a.cd, "beta": body.beta}},
                  open(corridorf + ".tmp", "w"), indent=1)
        os.replace(corridorf + ".tmp", corridorf)
        with open(statusf, "w") as fh:
            fh.write(f"{name}\ngoal: {goal}\n"
                     f"beta: {body.beta:.3e} kg/m^2  (D={a.diameter:.0f} m iron, Cd={a.cd})\n"
                     f"resolution level: {st['level']} (grid {12 + 6*st['level']}^2)\n"
                     f"cells done: {len(st['results'])}\n"
                     f"capture: {summ['n_capture']}  impact: {summ['n_impact']}  skip: {summ['n_skip']}\n"
                     f"max v_inf with capture: {summ['v_inf_max_capture_kms']:.2f} km/s\n"
                     f"widest peri window: {summ['widest_window_km']:.1f} km "
                     f"@ v_inf={summ['widest_window_at_vinf_kms']:.2f} km/s\n"
                     f"updated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        json.dump({"name": name, "goal": goal, "metric": metric,
                   "value": round(summ["v_inf_max_capture_kms"], 3),
                   "detail": (f"corridor: {summ['n_capture']} capture / "
                              f"{summ['n_impact']} impact / {summ['n_skip']} skip; "
                              f"window {summ['widest_window_km']:.0f} km"),
                   "iterations": len(st["results"]), "solved": False,
                   "updated": time.time()}, open(statusj, "w"), indent=1)

    deadline = None if a.daemon else time.time() + (a.increment or 1800)
    log(f"start pid={os.getpid()} level={st['level']} cell={st['cell']} "
        f"beta={body.beta:.3e} mode={'daemon' if a.daemon else f'inc {a.increment}s'}")

    vs, ps = _grid(st["level"])
    n = len(vs)
    while True:
        if st["cell"] >= n * n:
            # level complete -> refine
            log(f"level {st['level']} complete: {_summary(st['results'])}")
            if st["level"] >= MAX_LEVEL:
                save()
                log("max level reached; idling" if a.daemon else "max level; done")
                if not a.daemon:
                    return
                time.sleep(30)
                continue
            st["level"] += 1
            st["cell"] = 0
            st["results"] = []
            vs, ps = _grid(st["level"])
            n = len(vs)
        i, j = divmod(st["cell"], n)
        try:
            res = physics.simulate_passage(vs[i], ps[j], body)
            st["results"].append(res)
        except Exception as e:
            log(f"cell ({i},{j}) v={vs[i]:.0f} p={ps[j]:.0f} crash: {e}")
        st["cell"] += 1
        if st["cell"] % 25 == 0:
            save()
        if deadline is not None and time.time() >= deadline:
            save()
            log(f"increment done: level={st['level']} cells={len(st['results'])}")
            return


if __name__ == "__main__":
    run()
