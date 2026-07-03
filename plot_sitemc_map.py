#!/usr/bin/env python3
"""Flashlight diagram: global impact-location map from the site-MC catalog.

Each dot = one deorbit impact (capture -> aerobrake cascade -> terminal
plunge). Default = the raw isotropic-prior catalog. With --weights FILE
(csv: seed,weight), dots are REJECTION-SAMPLED (accept w/w_max) so point
density is an honest unweighted draw from the weighted distribution —
the right way to visualize a reweighted MC (weights stay reserved for
quantitative estimators in analyze_sitemc.py).

Output: figures/sitemc_flashlight.png/.pdf
Usage: uv run python plot_sitemc_map.py [--weights w.csv] [--out name]
"""
import csv
import math
import os
import random
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "src"))
CATALOG = os.path.join(HERE, "state", "sitemc", "catalog.csv")

INK = "#e8e4d8"
BG = "#101418"


def load_impacts():
    rows = list(csv.DictReader(open(CATALOG)))
    de = [r for r in rows if r["outcome"] == "impact" and int(r["n_pass"]) >= 2]
    return rows, de


def rejection(de, wfile):
    w = {}
    for r in csv.DictReader(open(wfile)):
        w[int(r["seed"])] = float(r["weight"])
    wmax = max(w.values())
    rng = random.Random(42)
    keep = [r for r in de if rng.random() < w.get(int(r["seed"]), 0.0) / wmax]
    print(f"rejection: kept {len(keep)}/{len(de)} (acceptance {len(keep)/len(de):.2%})")
    return keep


def main():
    _, de = load_impacts()
    tag = "isotropic arrival prior"
    if "--weights" in sys.argv:
        de = rejection(de, sys.argv[sys.argv.index("--weights") + 1])
        tag = "reweighted arrival prior (rejection-sampled)"
    graz = [r for r in de if r["fpa_deg"] and abs(float(r["fpa_deg"])) < 5.0]
    steep = [r for r in de if r not in graz]
    print(f"deorbit impacts: {len(de)} (grazing {len(graz)}, steep {len(steep)})")

    # DEM background (our own coastlines; stride-12 ~ 0.2 deg)
    from deorbit.sitemc.terrain import NC, NR
    dem = np.memmap(os.path.join(HERE, "data", "dem", "etopo1_ice_g_i2.bin"),
                    dtype="<i2", mode="r", shape=(NR, NC))
    st = 12
    sub = np.asarray(dem[::st, ::st], dtype=np.float32)
    lats = np.radians(90.0 - st * np.arange(sub.shape[0]) / 60.0)
    lons = np.radians(-180.0 + st * np.arange(sub.shape[1]) / 60.0)
    LO, LA = np.meshgrid(lons, lats)

    fig = plt.figure(figsize=(13.5, 6.4), facecolor=BG)
    ax = fig.add_axes([0.02, 0.05, 0.76, 0.88], projection="hammer", facecolor=BG)
    ax.set_xticklabels([]); ax.set_yticklabels([])
    ax.grid(alpha=0.12, color=INK, lw=0.4)
    land = np.ma.masked_less(sub, 0.0)
    ax.pcolormesh(LO, LA, land, cmap="Greys_r", vmin=-4000, vmax=9000,
                  rasterized=True, alpha=0.55, shading="auto")
    ax.contour(LO, LA, sub, levels=[0.0], colors="#4a5560", linewidths=0.35)

    # dots colored by impact obliquity |fpa| (log scale: 0.3-70 deg spans the
    # grazing cluster and the Moon-pumped steep tail); steepest drawn on top
    from matplotlib.colors import LogNorm
    lo = np.radians([float(r["lon_e_deg"]) for r in de])
    la = np.radians([float(r["lat_deg"]) for r in de])
    fpa = np.abs([float(r["fpa_deg"]) for r in de])
    fpa = np.clip(fpa, 0.3, 70.0)
    order = np.argsort(fpa)
    sc = ax.scatter(lo[order], la[order], s=9, c=fpa[order],
                    cmap="viridis", norm=LogNorm(vmin=0.3, vmax=70.0),
                    alpha=0.85, lw=0, zorder=5)
    cax = fig.add_axes([0.055, 0.115, 0.20, 0.022])
    cb = fig.colorbar(sc, cax=cax, orientation="horizontal")
    cb.set_label("impact flight-path angle |fpa| (deg)", color=INK, fontsize=8)
    cb.set_ticks([0.3, 1, 3, 10, 30, 70])
    cb.set_ticklabels(["0.3", "1", "3", "10", "30", "70"])
    cb.ax.tick_params(colors=INK, labelsize=7)
    cb.outline.set_edgecolor("#3a4148")

    from deorbit.site_selection import PEAKS
    plo = np.radians([p[2] for p in PEAKS]); pla = np.radians([p[1] for p in PEAKS])
    ax.scatter(plo, pla, marker="^", s=26, c="#5ad7ff", lw=0.4,
               edgecolors=BG, label="named summits", zorder=6)

    ax.legend(loc="lower right", fontsize=8, facecolor=BG, edgecolor="none",
              labelcolor=INK, framealpha=0.4)
    ax.set_title(f"km iron monolith — deorbit impact sites ({tag})",
                 color=INK, fontsize=12, pad=14)

    # latitude marginal
    axh = fig.add_axes([0.82, 0.13, 0.15, 0.72], facecolor=BG)
    la_all = [float(r["lat_deg"]) for r in de]
    bins = np.arange(-90, 91, 5)
    axh.hist(la_all, bins=bins, orientation="horizontal", color="#ffd76a", alpha=0.85)
    null = np.cos(np.radians((bins[:-1] + bins[1:]) / 2.0))
    null *= len(la_all) * 5.0 / np.degrees(2.0)   # cos-lat area null, same norm
    axh.plot(null, (bins[:-1] + bins[1:]) / 2.0, color="#5ad7ff", lw=1.0,
             label="area null")
    axh.set_ylim(-90, 90)
    axh.set_yticks([-60, -30, 0, 30, 60])
    axh.tick_params(colors=INK, labelsize=7)
    for s in axh.spines.values():
        s.set_color("#3a4148")
    axh.set_title("impact latitude", color=INK, fontsize=9)
    axh.legend(fontsize=7, facecolor=BG, edgecolor="none", labelcolor=INK)

    os.makedirs(os.path.join(HERE, "figures"), exist_ok=True)
    out = sys.argv[sys.argv.index("--out") + 1] if "--out" in sys.argv else "sitemc_flashlight"
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(HERE, "figures", f"{out}.{ext}"), dpi=180,
                    facecolor=BG, bbox_inches="tight")
    print(f"wrote figures/{out}.png/.pdf")


if __name__ == "__main__":
    main()
