#!/usr/bin/env python3
"""3D portrait of one capture->aerobrake->impact trajectory on a wireframe
globe. Re-flies a catalog seed deterministically with full-path tracking.

Default seed 402938: v_inf=1.5 km/s, i0=48 deg, apo1=68,392 km, 7 aerobraking
passes / 2.06 days, grazing impact (fpa -0.84 deg, 7.8 km/s) 78 km from
Cayambe — the Cayambe scenario realized spontaneously in the MC.

Two panels: (a) the full cascade (shrinking ellipses, J2 precession),
(b) zoom to <2.5 R_E: final passes + plunge over the coastline globe.
Coastlines drawn at the IMPACT epoch (inertial frame; earlier ground
tracks sweep under the rotating Earth).

Usage: uv run python plot_sitemc_traj.py [--seed N] [--out name]
"""
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d.art3d import Line3DCollection

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "src"))

from deorbit.corridor.physics import Body
from deorbit.sitemc import ephem, sample
from deorbit.sitemc.geodesy import A_EQ
from deorbit.sitemc.propagate import Config, fly
from deorbit.sitemc.terrain import Terrain

INK = "#1a1d21"


def refly(seed):
    d = sample.draw(seed)
    body = Body(cd=d["cd"])
    s0, geo = sample.entry_state(d["v_inf"], d["u_hat"], d["theta_B"],
                                 A_EQ + d["rp_alt_m"])
    cfg = Config(body.k, jd0=d["epoch_jd"], f_rho=d["f_rho"], f_H=d["f_H"],
                 terrain=Terrain())
    res = fly(s0, cfg, max_days=200.0, track_every=1)
    print(f"seed {seed}: {res['outcome']} passes={res['n_pass']} "
          f"apo1={res['apo1_km']:.0f} km t={res['t_days']:.2f} d "
          f"lat={res['latc_deg']:.2f} v={res['v_end']/1e3:.1f} km/s "
          f"fpa={res['fpa_deg']:.2f} deg  ({len(res['path'])} path points)")
    return d, geo, res


def coast_segments(gmst_deg):
    """Coastline polylines (from our DEM) on the unit sphere, rotated to the
    inertial frame by GMST."""
    from deorbit.sitemc.terrain import NC, NR
    dem = np.memmap(os.path.join(HERE, "data", "dem", "etopo1_ice_g_i2.bin"),
                    dtype="<i2", mode="r", shape=(NR, NC))
    st = 10
    sub = np.asarray(dem[::st, ::st], dtype=np.float32)
    lat = 90.0 - st * np.arange(sub.shape[0]) / 60.0
    lon = -180.0 + st * np.arange(sub.shape[1]) / 60.0
    f = plt.figure()
    cs = plt.contour(lon, lat, sub, levels=[0.0])
    segs = [np.array(s) for s in cs.allsegs[0] if len(s) > 20]
    plt.close(f)
    out = []
    for s in segs:
        lo = np.radians(s[:, 0] + gmst_deg)
        la = np.radians(s[:, 1])
        out.append(np.column_stack([np.cos(la) * np.cos(lo),
                                    np.cos(la) * np.sin(lo),
                                    np.sin(la)]) * 1.002)
    return out


def view_dir(elev_deg, azim_deg):
    """Unit vector from origin toward the (orthographic) camera, matching
    mplot3d's view_init convention."""
    el, az = math.radians(elev_deg), math.radians(azim_deg)
    return np.array([math.cos(el) * math.cos(az),
                     math.cos(el) * math.sin(az), math.sin(el)])


def visible_mask(P, vhat, R=1.0):
    """True where an exterior point is NOT occluded by the sphere of radius R
    (orthographic camera at +inf along vhat): occluded iff inside the
    silhouette cylinder AND behind the center plane."""
    d = P @ vhat
    perp2 = (P * P).sum(axis=1) - d * d
    return ~((d < 0.0) & (perp2 < R * R))


def _plot_runs(ax, poly, mask, **kw):
    p = poly.copy()
    p[~mask] = np.nan
    ax.plot(p[:, 0], p[:, 1], p[:, 2], **kw)


def draw_globe(ax, coasts, vhat):
    """Hidden-line globe: near-side wireframe + coastlines, faint far side,
    crisp limb circle. All culling done in data space (ortho camera)."""
    # wireframe as explicit polylines so we can cull per point
    th = np.linspace(0, 2 * np.pi, 181)
    lines = []
    for lam in np.radians(np.arange(0, 180, 15)):          # meridians
        lines.append(np.column_stack([np.cos(th) * np.cos(lam),
                                      np.cos(th) * np.sin(lam), np.sin(th)]))
    for phi in np.radians(np.arange(-75, 76, 15)):         # parallels
        lines.append(np.column_stack([np.cos(phi) * np.cos(th),
                                      np.cos(phi) * np.sin(th),
                                      np.full_like(th, np.sin(phi))]))
    for ln in lines:
        m = visible_mask(ln, vhat)
        _plot_runs(ax, ln, ~m, color="#b9c0c7", lw=0.25, alpha=0.18, zorder=1)
        _plot_runs(ax, ln, m, color="#b9c0c7", lw=0.3, alpha=0.65, zorder=2)
    for c in coasts:
        _plot_runs(ax, c, visible_mask(c, vhat), color="#4a5560", lw=0.6, zorder=3)
    # limb (silhouette) circle: plane through origin perpendicular to vhat
    e1 = np.cross(vhat, [0.0, 0.0, 1.0])
    if np.linalg.norm(e1) < 1e-6:
        e1 = np.cross(vhat, [1.0, 0.0, 0.0])
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(vhat, e1)
    limb = np.outer(np.cos(th), e1) + np.outer(np.sin(th), e2)
    ax.plot(limb[:, 0], limb[:, 1], limb[:, 2], color="#8a939c", lw=0.8, zorder=3)


def traj_collection(P, tdays, norm, vhat, decimate=1):
    """Time-colored trajectory with sphere-occluded segments removed."""
    pts = P[::decimate]
    td = tdays[::decimate]
    vis = visible_mask(pts, vhat)
    keep = vis[:-1] & vis[1:]                  # segment visible if both ends are
    segs = np.stack([pts[:-1], pts[1:]], axis=1)[keep]
    lc = Line3DCollection(segs, cmap="managua", norm=norm, lw=0.9, zorder=4)
    lc.set_array((0.5 * (td[:-1] + td[1:]))[keep])
    return lc


def main():
    seed = int(sys.argv[sys.argv.index("--seed") + 1]) if "--seed" in sys.argv else 402938
    d, geo, res = refly(seed)
    path = np.array(res["path"])
    t = path[:, 0]
    P = path[:, 1:4] / A_EQ
    tdays = t / 86400.0
    jd_end = d["epoch_jd"] + res["t_end"] / 86400.0
    gmst_deg = math.degrees(ephem.gmst_rad(jd_end)) % 360.0
    coasts = coast_segments(gmst_deg)
    imp = P[-1]
    norm = Normalize(vmin=0.0, vmax=tdays[-1])

    fig = plt.figure(figsize=(12.5, 6.6), facecolor="white")
    fig.suptitle(f"capture → aerobrake → grazing impact  (seed {seed}: "
                 f"v∞={d['v_inf']/1e3:.1f} km/s, i₀={geo['i_deg']:.0f}°, "
                 f"{res['n_pass']} passes, apo₁={res['apo1_km']:.0f} km, "
                 f"{res['t_days']:.2f} d, impact {res['v_end']/1e3:.1f} km/s "
                 f"@ {res['fpa_deg']:.2f}°)", color=INK, fontsize=11)

    az_imp = math.degrees(math.atan2(imp[1], imp[0]))
    for k, (dec, ttl) in enumerate(((2, "full cascade"), (1, "endgame (zoom)"))):
        ax = fig.add_subplot(1, 2, k + 1, projection="3d", facecolor="white")
        ax.set_proj_type("ortho")          # culling math assumes ortho camera
        ax.computed_zorder = False         # we own the paint order
        elev, azim = 18.0, az_imp + (25.0 if k == 0 else 5.0)
        vhat = view_dir(elev, azim)
        draw_globe(ax, coasts, vhat)
        if k == 0:
            sel = np.ones(len(P), bool)
            lo_ = np.minimum(P.min(axis=0), -1.0)
            hi_ = np.maximum(P.max(axis=0), 1.0)
            ctr = 0.5 * (lo_ + hi_)
            half = 0.54 * (hi_ - lo_).max()
        else:
            sel = np.linalg.norm(P, axis=1) < 3.0
            ctr = np.zeros(3)
            half = 2.1
        lc = traj_collection(P[sel], tdays[sel], norm, vhat, decimate=dec)
        lc.set_linewidth(1.1 if k == 0 else 0.9)
        ax.add_collection3d(lc)
        if visible_mask(imp[None, :], vhat)[0]:
            ax.scatter(*imp, color="#d62728", marker="*", s=90, zorder=5)
        if k == 0:
            ax.scatter(*P[0], color=INK, marker="o", s=12, zorder=5)
            ax.text(P[0][0], P[0][1] - 0.4 * half, P[0][2] - 0.12 * half,
                    "arrival", color=INK, fontsize=8, ha="center", zorder=5)
        else:
            ax.text(*imp, "  impact", color="#d62728", fontsize=9, zorder=5)
        ax.set_box_aspect((1, 1, 1))
        ax.set_xlim(ctr[0] - half, ctr[0] + half)
        ax.set_ylim(ctr[1] - half, ctr[1] + half)
        ax.set_zlim(ctr[2] - half, ctr[2] + half)
        ax.view_init(elev=elev, azim=azim)
        ax.set_axis_off()
        ax.set_title(ttl, color=INK, fontsize=10, pad=0, y=0.93)
    fig.subplots_adjust(left=-0.04, right=1.04, top=0.99, bottom=0.0, wspace=-0.18)

    cax = fig.add_axes([0.41, 0.10, 0.18, 0.02])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap="managua"),
                      cax=cax, orientation="horizontal")
    cb.set_label("days since arrival", color=INK, fontsize=8)
    cb.ax.tick_params(colors=INK, labelsize=7)
    cb.outline.set_edgecolor("#c3c9cf")

    out = sys.argv[sys.argv.index("--out") + 1] if "--out" in sys.argv else "sitemc_traj"
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(HERE, "figures", f"{out}.{ext}"), dpi=180,
                    facecolor="white", bbox_inches="tight")
    print(f"wrote figures/{out}.png/.pdf")


if __name__ == "__main__":
    main()
