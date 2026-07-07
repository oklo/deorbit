#!/usr/bin/env python3
"""Collins et al. (2020) Fig. 2-style crater-development frames from the
crater mode's tracer + snapshot output (CRATER_TRACE=1, CRATER_SNAP_DT>0).

Each panel is a full cross-section (the axisymmetric solution mirrored about
r=0): the pre-impact target carried by a Lagrangian tracer lattice drawn as
deformed mesh lines; the uppermost 'sediment analog' band in sandy brown;
tracers shocked above the melt threshold in red; a mid band of tracers
shaded by recorded peak pressure (white->blue). The free surface comes from
the density field.

Usage:
  uv run python plot_crater_frames.py state/trace/a3000 --times 100 800 2400 4600 \
      [--melt-gpa 60] [--sed-km 3] [--out figures/frames_a3000]
"""
import argparse
import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize

INK = "#1a1d21"


def load_snap(path):
    hdr = open(path).readline()
    t = float(hdr.split("t=")[1].split()[0])
    nx = int(hdr.split("nx=")[1].split()[0])
    nz = int(hdr.split("nz=")[1].split()[0])
    dx = float(hdr.split("dx=")[1].split()[0])
    d = np.loadtxt(path)
    F = {n: d[:, 2 + j].reshape(nx, nz) for j, n in
         enumerate(["rho", "P", "vib", "af", "D", "Pmax"])}
    return t, nx, nz, dx, F


def load_tracers(path):
    hdr = open(path).readline()
    ksurf = int(hdr.split("ksurf=")[1].split()[0])
    d = np.loadtxt(path)
    return ksurf, d[:, 1].astype(int), d[:, 2].astype(int), d[:, 3], d[:, 4], d[:, 5]


def nearest_file(dirpath, stem, t):
    files = glob.glob(os.path.join(dirpath, f"{stem}_*.txt"))
    ts = np.array([int(os.path.basename(f).split("_")[1].split(".")[0]) for f in files])
    return files[int(np.argmin(np.abs(ts - t)))]


def panel(ax, dirpath, t_req, melt_pa, sed_m, mesh_every, zsurf_m):
    sf = nearest_file(dirpath, "snap", t_req)
    tf = nearest_file(dirpath, "tracers", t_req)
    t, nx, nz, dx, F = load_snap(sf)
    ksurf, i0, k0, r, z, pm = load_tracers(tf)

    km = 1e-3
    zk = (z - zsurf_m) * km                     # depth-referenced z (0 = pre-impact surface)
    # classes by ORIGINAL position
    depth0 = zsurf_m - (k0 + 0.5) * dx
    sed = depth0 <= sed_m
    melt = pm >= melt_pa
    shock_band = (~sed) & (~melt)

    for sign in (1.0, -1.0):
        rr = sign * r * km
        # deformed mesh: rows (constant k0) and columns (constant i0), every Nth
        for kk in range(0, k0.max() + 1, mesh_every):
            m = (k0 == kk) & ~sed
            if m.sum() > 3:
                o = np.argsort(i0[m])
                ax.plot(rr[m][o], zk[m][o], color="0.25", lw=0.3, zorder=3)
        for ii in range(0, i0.max() + 1, mesh_every):
            m = (i0 == ii) & ~sed
            if m.sum() > 3:
                o = np.argsort(k0[m])
                ax.plot(rr[m][o], zk[m][o], color="0.25", lw=0.3, zorder=3)
        # peak-shock shading (white->blue) for moderately shocked tracers
        hot = shock_band & (pm > 1e9)
        ax.scatter(rr[hot], zk[hot], s=0.9, c=pm[hot] / 1e9, cmap="Blues",
                   norm=Normalize(0, melt_pa / 1e9), lw=0, zorder=4)
        # sediment-analog band and melt
        ax.scatter(rr[sed & ~melt], zk[sed & ~melt], s=0.9, color="#d8a35c", lw=0, zorder=5)
        ax.scatter(rr[melt], zk[melt], s=1.1, color="#c22a1e", lw=0, zorder=6)

    # target silhouette from the density field (grey under the free surface)
    rho = F["rho"]
    ii = (np.arange(nx) + 0.5) * dx * km
    surf_k = np.array([np.max(np.where(rho[i] > 1350)[0]) if (rho[i] > 1350).any()
                       else 0 for i in range(nx)])
    surf_z = ((surf_k + 1.0) * dx - zsurf_m) * km
    for sign in (1.0, -1.0):
        ax.fill_between(sign * ii, -zsurf_m * km, surf_z, color="0.88", zorder=1)
        ax.plot(sign * ii, surf_z, color="0.4", lw=0.5, zorder=2)
    ax.set_title(f"T = {t:.0f} s", fontsize=9, loc="right", color=INK)
    return dx


def main():
    p = argparse.ArgumentParser()
    p.add_argument("dirpath")
    p.add_argument("--times", type=float, nargs="+", required=True)
    p.add_argument("--melt-gpa", type=float, default=60.0)
    p.add_argument("--sed-km", type=float, default=3.0)
    p.add_argument("--mesh-every", type=int, default=5)
    p.add_argument("--rmax-km", type=float, default=0.0)
    p.add_argument("--zmin-km", type=float, default=0.0)
    p.add_argument("--out", default="")
    a = p.parse_args()

    # zsurf from the run geometry: Zdom - 5*a_imp; recover from any snapshot header
    sf0 = sorted(glob.glob(os.path.join(a.dirpath, "snap_*.txt")))[0]
    _, nx, nz, dx, _ = load_snap(sf0)
    ksurf, *_ = load_tracers(sorted(glob.glob(os.path.join(a.dirpath, "tracers_*.txt")))[0])
    zsurf_m = ksurf * dx

    n = len(a.times)
    fig, axes = plt.subplots(n, 1, figsize=(9.5, 2.6 * n), facecolor="white",
                             sharex=True, sharey=True)
    axes = np.atleast_1d(axes)
    for ax, tt in zip(axes, a.times):
        panel(ax, a.dirpath, tt, a.melt_gpa * 1e9, a.sed_km * 1e3, a.mesh_every, zsurf_m)
        rmax = a.rmax_km or 0.85 * nx * dx * 1e-3
        zmin = -(a.zmin_km or zsurf_m * 1e-3)
        ax.set_xlim(-rmax, rmax)
        ax.set_ylim(zmin, 0.35 * (-zmin))
        ax.set_ylabel("z [km]", fontsize=8)
        ax.tick_params(labelsize=7, colors=INK)
    axes[-1].set_xlabel("r [km]", fontsize=8)
    sm = plt.cm.ScalarMappable(norm=Normalize(0, a.melt_gpa), cmap="Blues")
    cb = fig.colorbar(sm, ax=axes, shrink=0.5, pad=0.01)
    cb.set_label("tracer peak pressure [GPa]", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    out = a.out or os.path.join("figures", "frames_" + os.path.basename(a.dirpath.rstrip("/")))
    for ext in ("png", "pdf"):
        fig.savefig(f"{out}.{ext}", dpi=150, facecolor="white", bbox_inches="tight")
    print(f"wrote {out}.png/.pdf")


if __name__ == "__main__":
    main()
