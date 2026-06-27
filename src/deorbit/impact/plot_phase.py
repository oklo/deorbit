#!/usr/bin/env python3
"""Render the bounce-off (theta, v/v_esc) phase diagram from results/bounce_off_phase_diagram.json
(produced by `run_bounce.py --report`). FIGURE tool -- uses matplotlib (not the stdlib daemon core).

    python3 -m deorbit.impact.plot_phase        # -> figures/bounce_off_phase_diagram.png/.pdf

Primary visual = the 69 ground-truth SPH runs colored by regime (C/I/D/E). Overlaid: the
bisection-located boundary crossings (black markers) and simple guide curves through them, plus
the two anchor cases (Cayambe, Orcus). First-draft layout; refine before publication.
"""
import os
import json

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

_REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

# regime -> (color, marker, label)
STYLE = {"C": ("#6e4b2a", "s", "crater (retained)"),
         "I": ("#1f6fb2", "o", "intact ricochet"),
         "D": ("#e8862e", "^", "dispersed ricochet"),
         "E": ("#c0392b", "*", "escape")}


def _curve(points):
    """Sort (x,y) points by x for a guide polyline."""
    p = sorted(points)
    return [q[0] for q in p], [q[1] for q in p]


def main():
    src = os.path.join(_REPO, "results", "bounce_off_phase_diagram.json")
    d = json.load(open(src))
    runs, bnds = d["runs"], d["boundaries"]

    fig, ax = plt.subplots(figsize=(8.5, 6.2))

    # --- ground-truth SPH runs, by regime ---
    for reg, (c, m, lab) in STYLE.items():
        xs = [r["theta"] for r in runs if r["regime"] == reg]
        ys = [r["vr"] for r in runs if r["regime"] == reg]
        ax.scatter(xs, ys, c=c, marker=m, s=70, edgecolors="k", linewidths=0.4,
                   zorder=3, label=f"{lab}  (n={len(xs)})")

    # --- bisection-located boundary crossings + guide curves ---
    id_pts, de_pts, crater_th = [], [], []
    for b in bnds:
        # each boundary is a transition point in the (theta, v/v_esc) plane
        x = b["value"] if b["vary"] == "theta" else b["fixed"]
        y = b["fixed"] if b["vary"] == "theta" else b["value"]
        trans = {b["from"], b["to"]}
        if trans == {"I", "D"}:
            id_pts.append((x, y))
        elif trans == {"D", "E"}:
            de_pts.append((x, y))
        elif trans == {"D", "C"}:
            crater_th.append(x)
        ax.scatter([x], [y], facecolors="none", edgecolors="k", marker="D",
                   s=120, linewidths=1.3, zorder=4)

    if id_pts:
        cx, cy = _curve(id_pts)
        ax.plot(cx, cy, "--", color="#1f6fb2", lw=1.8, zorder=2, label="intact↔dispersed (melt onset)")
    if de_pts:
        cx, cy = _curve(de_pts)
        ax.plot(cx, cy, "--", color="#c0392b", lw=1.8, zorder=2, label="dispersed↔escape")
    if crater_th:
        thc = float(np.mean(crater_th))
        ax.axvline(thc, ls="-.", color="#6e4b2a", lw=2.0, zorder=2,
                   label=f"crater onset θ≈{thc:.0f}°")

    # --- anchor cases ---
    anchors = [("Cayambe\n(Earth, 7.4 km/s, 2°)", 2.0, 7400.0 / d["v_esc_earth"], (14, 12)),
               ("Orcus-like\n(Mars, 10 km/s, 30°)", 30.0, 10000.0 / 5027.0, (-10, -28))]
    for name, th, vr, off in anchors:
        ax.scatter([th], [vr], marker="P", s=240, c="gold", edgecolors="k", linewidths=1.2, zorder=6)
        ax.annotate(name, (th, vr), textcoords="offset points", xytext=off, fontsize=9,
                    ha="center", fontweight="bold",
                    arrowprops=dict(arrowstyle="->", lw=1.0))

    ax.axhline(1.0, color="grey", lw=0.8, ls=":", zorder=1)
    ax.text(44.5, 1.02, "v = v_esc", ha="right", fontsize=8, color="grey")
    ax.set_xlabel("impact angle θ from horizontal (deg)", fontsize=12)
    ax.set_ylabel("v / v_esc", fontsize=12)
    ax.set_title("Grazing-impact bounce-off phase diagram — km iron on rock\n"
                 f"({d['n_runs']} GPU-SPH runs, CPPR=8; {d['converged']}/{d['n_brackets']} boundaries bisected)",
                 fontsize=12)
    ax.set_xlim(0, 45); ax.set_ylim(0.3, 2.6)
    ax.legend(loc="upper right", fontsize=8, framealpha=0.95, ncol=1)
    ax.grid(alpha=0.25)
    fig.tight_layout()

    outdir = os.path.join(_REPO, "figures")
    os.makedirs(outdir, exist_ok=True)
    png = os.path.join(outdir, "bounce_off_phase_diagram.png")
    fig.savefig(png, dpi=150)
    fig.savefig(os.path.join(outdir, "bounce_off_phase_diagram.pdf"))
    print(f"wrote {png}\n      {png[:-4]}.pdf")


if __name__ == "__main__":
    main()
