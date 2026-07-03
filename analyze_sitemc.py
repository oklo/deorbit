#!/usr/bin/env python3
"""Site-MC catalog analysis (Phase D; docs/site_mc_plan.md).

Headline statistic: eta(elevation class) = P_impact(class) / P_area(class),
where the null is the cos(lat)-area-weighted DEM elevation distribution
REWEIGHTED TO THE IMPACT LATITUDE DISTRIBUTION (2-deg bins) — so eta measures
elevation steering AT FIXED latitude, separating the summit question from the
(real, expected) equatorward latitude steering. Bootstrap CIs over impacts.

Usage: python3 analyze_sitemc.py [state/sitemc/catalog.csv] [--json results/sitemc_summary.json]
"""
import csv
import json
import math
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "src"))

EDGES = [0.0, 1000.0, 3000.0, 5000.0]        # elevation classes (m): <=0 ocean, then bands
CLASS_NAMES = ["ocean", "0-1km", "1-3km", "3-5km", ">5km"]
LATBIN = 2.0


def elev_class(e):
    i = 0
    for k, edge in enumerate(EDGES):
        if e > edge:
            i = k + 1
    return i


def latitude_controlled_null(terr, lats, stride=6):
    """Area elevation distribution weighted to the impact-latitude histogram."""
    from deorbit.sitemc.terrain import NC, NR, CELL
    import struct
    hist = {}
    for la in lats:
        hist[int(la // LATBIN)] = hist.get(int(la // LATBIN), 0) + 1
    w = [0.0] * len(CLASS_NAMES)
    u = struct.unpack_from
    mm = terr.mm
    for key, cnt in hist.items():
        lat_lo, lat_hi = key * LATBIN, (key + 1) * LATBIN
        r_lo = max(0, int((90.0 - lat_hi) / CELL))
        r_hi = min(NR - 1, int((90.0 - lat_lo) / CELL))
        band = [0.0] * len(CLASS_NAMES)
        tot = 0.0
        for r in range(r_lo, r_hi + 1, stride):
            for c in range(0, NC - 1, stride):
                v = u("<h", mm, 2 * (r * NC + c))[0]
                band[elev_class(v)] += 1.0
                tot += 1.0
        if tot:
            for i in range(len(band)):
                w[i] += cnt * band[i] / tot
    s = sum(w)
    return [x / s for x in w] if s else w


def main():
    cat = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("--") \
        else os.path.join(HERE, "state", "sitemc", "catalog.csv")
    rows = list(csv.DictReader(open(cat)))
    n = len(rows)
    out = {"n_seeds": n}
    by = {}
    for r in rows:
        k = r["outcome"]
        if k == "impact":
            k = "impact_direct" if int(r["n_pass"]) < 2 else "impact_deorbit"
        by[k] = by.get(k, 0) + 1
    out["outcomes"] = by
    print(f"seeds: {n}   outcomes: {by}")

    deorb = [r for r in rows if r["outcome"] == "impact" and int(r["n_pass"]) >= 2]
    out["n_deorbit"] = len(deorb)
    if not deorb:
        print("no deorbit impacts yet.")
        _dump(out)
        return

    lats = [float(r["lat_deg"]) for r in deorb]
    abs_lats = sorted(abs(l) for l in lats)
    out["abs_lat_quartiles"] = [round(abs_lats[int(q * (len(abs_lats) - 1))], 2)
                                for q in (0.25, 0.5, 0.75, 0.9)]
    elevs = [float(r["elev_m"]) if r["elev_m"] else 0.0 for r in deorb]
    imp_frac = [0.0] * len(CLASS_NAMES)
    for e in elevs:
        imp_frac[elev_class(e if e > 0 else -1.0)] += 1.0
    imp_frac = [x / len(elevs) for x in imp_frac]

    from deorbit.sitemc.terrain import Terrain
    terr = Terrain()
    null = latitude_controlled_null(terr, lats)
    out["impact_frac"] = [round(x, 4) for x in imp_frac]
    out["null_frac_latctrl"] = [round(x, 4) for x in null]

    # bootstrap CI on eta
    rng = random.Random(1)
    B = 400
    etas = {i: [] for i in range(len(CLASS_NAMES))}
    for _ in range(B):
        samp = [elevs[rng.randrange(len(elevs))] for _ in elevs]
        f = [0.0] * len(CLASS_NAMES)
        for e in samp:
            f[elev_class(e if e > 0 else -1.0)] += 1.0
        for i in range(len(CLASS_NAMES)):
            if null[i] > 0:
                etas[i].append(f[i] / len(samp) / null[i])
    print(f"\ndeorbit impacts: {len(deorb)}   |lat| quartiles(25/50/75/90): "
          f"{out['abs_lat_quartiles']}")
    print(f"{'class':8s} {'P_impact':>9s} {'P_area|lat':>11s} {'eta':>7s}   95% CI")
    out["eta"] = {}
    for i, name in enumerate(CLASS_NAMES):
        if null[i] > 0:
            es = sorted(etas[i])
            eta = imp_frac[i] / null[i]
            lo, hi = es[int(0.025 * B)], es[int(0.975 * B)]
            out["eta"][name] = [round(eta, 3), round(lo, 3), round(hi, 3)]
            print(f"{name:8s} {imp_frac[i]:9.4f} {null[i]:11.4f} {eta:7.2f}   "
                  f"[{lo:.2f}, {hi:.2f}]")

    # named peaks
    near = {}
    for r in deorb:
        if r["peak_name"] and r["peak_dist_km"] and float(r["peak_dist_km"]) < 150.0:
            near[r["peak_name"]] = near.get(r["peak_name"], 0) + 1
    out["impacts_within_150km_of_peak"] = near
    print(f"\nimpacts within 150 km of a named summit: {near or 'none'}")

    # grazing + Moon-filter checks
    graz = sum(1 for r in deorb if r["fpa_deg"] and abs(float(r["fpa_deg"])) < 5.0)
    gentle = sum(1 for r in deorb if r["apo1_km"] and float(r["apo1_km"]) < 130000.0)
    out["grazing_frac"] = round(graz / len(deorb), 3)
    out["apo1_below_130k_frac"] = round(gentle / len(deorb), 3)
    land = sum(1 for e in elevs if e > 0)
    out["land_frac"] = round(land / len(elevs), 3)
    print(f"grazing(|fpa|<5deg): {graz}/{len(deorb)}   apo1<130k km: {gentle}/{len(deorb)}"
          f"   land fraction: {out['land_frac']}")
    _dump(out)


def _dump(out):
    if "--json" in sys.argv:
        p = sys.argv[sys.argv.index("--json") + 1]
        json.dump(out, open(p, "w"), indent=1)
        print(f"wrote {p}")


if __name__ == "__main__":
    main()
