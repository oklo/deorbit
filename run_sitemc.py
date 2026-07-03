#!/usr/bin/env python3
"""Site-selection Monte Carlo driver (pure stdlib; docs/site_mc_plan.md).

Each seed = one fully-recorded prior draw (deorbit.sitemc.sample.draw) flown
end-to-end (capture pass -> cascade -> terminal plunge -> DEM intercept).
Rows append to state/sitemc/catalog.csv; idempotent by seed (re-invoking
continues). Priors are applied post-hoc by reweighting the catalog.

Usage:
  python3 run_sitemc.py --n 20000 --jobs 7          # grind seeds 0..19999
  python3 run_sitemc.py --n 50 --jobs 2 --smoke     # quick smoke run
Launch detached:  nohup caffeinate -i -m -s python3 run_sitemc.py ... &
"""
import csv
import math
import multiprocessing as mp
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "src"))

STATE_DIR = os.path.join(HERE, "state", "sitemc")
CATALOG = os.path.join(STATE_DIR, "catalog.csv")

FIELDS = ["seed", "outcome", "v_inf", "rp_alt_km", "b_km", "i0_deg", "node0_deg",
          "asym_dec_deg", "cd", "f_rho", "f_H", "epoch_jd", "n_pass", "n_low_pass",
          "apo1_km", "t_days", "lat_deg", "lon_e_deg", "elev_m", "v_imp", "fpa_deg",
          "peak_name", "peak_dist_km", "min_alt_km", "steps"]

_worker = {}


def _init_worker():
    from deorbit.corridor.physics import RHO_IRON
    from deorbit.sitemc.terrain import Terrain
    _worker["terrain"] = Terrain()


def fly_seed(seed):
    import math as m
    from deorbit.corridor.physics import Body
    from deorbit.sitemc import sample
    from deorbit.sitemc.geodesy import A_EQ, B_POL
    from deorbit.sitemc.propagate import Config, fly

    d = sample.draw(seed)
    body = Body(cd=d["cd"])
    s0, geo = sample.entry_state(d["v_inf"], d["u_hat"], d["theta_B"],
                                 A_EQ + d["rp_alt_m"])
    cfg = Config(body.k, jd0=d["epoch_jd"], f_rho=d["f_rho"], f_H=d["f_H"],
                 terrain=_worker.get("terrain"))
    r = fly(s0, cfg, max_days=200.0)

    row = {
        "seed": seed, "outcome": r["outcome"], "v_inf": round(d["v_inf"], 1),
        "rp_alt_km": round(d["rp_alt_m"] / 1e3, 3), "b_km": round(geo["b_m"] / 1e3, 1),
        "i0_deg": round(geo["i_deg"], 3), "node0_deg": round(geo["node_deg"], 2),
        "asym_dec_deg": round(d["asym_dec_deg"], 2), "cd": round(d["cd"], 3),
        "f_rho": round(d["f_rho"], 3), "f_H": round(d["f_H"], 3),
        "epoch_jd": round(d["epoch_jd"], 3), "n_pass": r["n_pass"],
        "n_low_pass": r["n_low_pass"],
        "apo1_km": round(r["apo1_km"], 0) if r["apo1_km"] else "",
        "t_days": round(r["t_days"], 3), "lat_deg": "", "lon_e_deg": "",
        "elev_m": "", "v_imp": "", "fpa_deg": "", "peak_name": "",
        "peak_dist_km": "", "min_alt_km": round(r["min_alt_km"], 2),
        "steps": r["steps"],
    }
    if r["outcome"] == "impact":
        latc = m.radians(r["latc_deg"])
        lat_d = m.degrees(m.atan2(m.sin(latc), (B_POL / A_EQ) ** 2 * m.cos(latc)))
        lon = ((r["lon_e_deg"] + 180.0) % 360.0) - 180.0
        terr = _worker.get("terrain")
        elev = terr.elevation(lat_d, lon) if terr else 0.0
        pk_name, pk_dist = "", ""
        if terr:
            best = None
            for name, plat, plon, pelev, pcos in terr.peaks:
                dlon = abs((lon - plon + 180.0) % 360.0 - 180.0)
                dd = 111.32 * m.hypot(lat_d - plat, dlon * pcos)
                if best is None or dd < best[1]:
                    best = (name, dd)
            if best and best[1] < 300.0:
                pk_name, pk_dist = best[0], round(best[1], 1)
        row.update({"lat_deg": round(lat_d, 3), "lon_e_deg": round(lon, 3),
                    "elev_m": round(max(elev, 0.0), 0), "v_imp": round(r["v_end"], 0),
                    "fpa_deg": round(r["fpa_deg"], 2), "peak_name": pk_name,
                    "peak_dist_km": pk_dist})
    return row


def done_seeds():
    if not os.path.exists(CATALOG):
        return set()
    with open(CATALOG) as f:
        return {int(r["seed"]) for r in csv.DictReader(f)}


def main():
    n = int(sys.argv[sys.argv.index("--n") + 1]) if "--n" in sys.argv else 1000
    jobs = int(sys.argv[sys.argv.index("--jobs") + 1]) if "--jobs" in sys.argv else 1
    os.makedirs(STATE_DIR, exist_ok=True)
    todo = [s for s in range(n) if s not in done_seeds()]
    print(f"sitemc: {n} requested, {n - len(todo)} done, {len(todo)} to fly, jobs={jobs}",
          flush=True)
    if not todo:
        return
    new_file = not os.path.exists(CATALOG)
    fout = open(CATALOG, "a", newline="")
    w = csv.DictWriter(fout, fieldnames=FIELDS)
    if new_file:
        w.writeheader()
        fout.flush()
    n_done = 0
    with mp.get_context("fork").Pool(jobs, initializer=_init_worker) as pool:
        for row in pool.imap_unordered(fly_seed, todo, chunksize=4):
            w.writerow(row)
            n_done += 1
            if n_done % 25 == 0:
                fout.flush()
                print(f"  {n_done}/{len(todo)} flown (last: seed {row['seed']} "
                      f"-> {row['outcome']})", flush=True)
    fout.flush()
    fout.close()
    print("sitemc: increment complete.", flush=True)


if __name__ == "__main__":
    _init_worker()          # for --jobs 1 / direct calls
    main()
