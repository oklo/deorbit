#!/usr/bin/env python3
"""3D portrait of one capture->aerobrake->impact trajectory: a TRUE
perspective render (pinhole camera, per-point line-of-sight occlusion).

mplot3d cannot do this honestly (parallel projection + whole-artist depth
sort), so this is a tiny software renderer:
  * finite eye point E, look-at target, pinhole projection X=x_c/z_c;
  * a trajectory point P is occluded iff the ray E->P intersects the
    Earth sphere strictly between eye and point (quadratic per point);
  * points ON the sphere (coastlines/graticule) are visible iff on the
    near cap, P.E > R^2; the far side of an opaque planet is not drawn;
  * the horizon (tangent-cone circle, P.Ehat = R^2/|E|) replaces the
    orthographic limb;
  * visible trajectory segments are depth-sorted far->near and drawn
    with mildly depth-scaled linewidths (perspective depth cue).

Default seed 402938: v_inf=1.5 km/s, i0=48 deg, apo1=68,392 km, 7 passes
/ 2.06 days, grazing impact (fpa -0.84 deg, 7.8 km/s) 78 km from Cayambe.
Coastlines drawn at the impact-epoch GMST (inertial frame).

Usage: uv run python plot_sitemc_traj.py [--seed N] [--out name]
"""
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize

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
                                    np.sin(la)]))
    return out


class Camera:
    """Pinhole camera: eye E, look-at target, world-z up."""

    def __init__(self, eye, target, up=(0.0, 0.0, 1.0)):
        self.E = np.asarray(eye, float)
        f = np.asarray(target, float) - self.E
        self.f = f / np.linalg.norm(f)
        r = np.cross(self.f, up)
        self.r = r / np.linalg.norm(r)
        self.u = np.cross(self.r, self.f)

    def project(self, P):
        """(n,3) world -> screen X,Y and depth z (along the optical axis)."""
        D = np.atleast_2d(P) - self.E
        z = D @ self.f
        X = (D @ self.r) / z
        Y = (D @ self.u) / z
        return X, Y, z

    def occluded(self, P, R=1.0):
        """Ray eye->P hits the sphere strictly between them (free points)."""
        D = np.atleast_2d(P) - self.E
        a = (D * D).sum(axis=1)
        b = 2.0 * (D @ self.E)
        c = self.E @ self.E - R * R
        disc = b * b - 4.0 * a * c
        hit = disc > 0.0
        s0 = np.where(hit, (-b - np.sqrt(np.maximum(disc, 0.0))) / (2.0 * a), np.inf)
        return hit & (s0 > 1e-9) & (s0 < 1.0 - 1e-9)

    def near_cap(self, P, R=1.0):
        """Visibility for points lying ON the sphere of radius ~R."""
        return (np.atleast_2d(P) @ self.E) > R * R

    def horizon(self, R=1.0, n=241):
        """Tangent-cone circle: the true perspective silhouette."""
        dE = np.linalg.norm(self.E)
        eh = self.E / dE
        c = eh * (R * R / dE)
        rh = R * math.sqrt(max(0.0, 1.0 - (R / dE) ** 2))
        e1 = np.cross(eh, [0.0, 0.0, 1.0])
        if np.linalg.norm(e1) < 1e-6:
            e1 = np.cross(eh, [1.0, 0.0, 0.0])
        e1 /= np.linalg.norm(e1)
        e2 = np.cross(eh, e1)
        th = np.linspace(0, 2 * np.pi, n)
        return c + rh * (np.outer(np.cos(th), e1) + np.outer(np.sin(th), e2))


def draw_sphere_lines(ax, cam, poly, R=1.0, **kw):
    X, Y, z = cam.project(poly)
    vis = cam.near_cap(poly, R) & (z > 1e-3)
    X = np.where(vis, X, np.nan)
    Y = np.where(vis, Y, np.nan)
    ax.plot(X, Y, **kw)
    return X, Y


def main():
    seed = int(sys.argv[sys.argv.index("--seed") + 1]) if "--seed" in sys.argv else 402938
    d, geo, res = refly(seed)
    path = np.array(res["path"])
    P = path[:, 1:4] / A_EQ
    tdays = path[:, 0] / 86400.0
    jd_end = d["epoch_jd"] + res["t_end"] / 86400.0
    gmst_deg = math.degrees(ephem.gmst_rad(jd_end)) % 360.0
    coasts = coast_segments(gmst_deg)
    imp = P[-1]
    norm = Normalize(vmin=0.0, vmax=tdays[-1])
    az_imp = math.atan2(imp[1], imp[0])

    fig = plt.figure(figsize=(12.5, 6.6), facecolor="white")
    fig.suptitle(f"capture → aerobrake → grazing impact  (seed {seed}: "
                 f"v∞={d['v_inf']/1e3:.1f} km/s, i₀={geo['i_deg']:.0f}°, "
                 f"{res['n_pass']} passes, apo₁={res['apo1_km']:.0f} km, "
                 f"{res['t_days']:.2f} d, impact {res['v_end']/1e3:.1f} km/s "
                 f"@ {res['fpa_deg']:.2f}°)", color=INK, fontsize=11)

    panels = []
    # (title, camera, point-selection, decimation)
    # endgame: eye 3.2 R_E out, 25 deg around from the impact meridian
    n_hat = np.array([math.cos(0.31) * math.cos(az_imp + 0.09),
                      math.cos(0.31) * math.sin(az_imp + 0.09), math.sin(0.31)])
    panels.append(("endgame (zoom)", Camera(3.2 * n_hat, (0, 0, 0)),
                   np.linalg.norm(P, axis=1) < 3.4, 1))
    # full cascade: look at the orbit-cloud centroid from ~2.4x its extent
    ctr = 0.5 * (P.min(axis=0) + P.max(axis=0))
    half = 0.5 * (P.max(axis=0) - P.min(axis=0)).max()
    n2 = np.array([math.cos(0.31) * math.cos(az_imp + 0.44),
                   math.cos(0.31) * math.sin(az_imp + 0.44), math.sin(0.31)])
    panels.insert(0, ("full cascade", Camera(ctr + 2.4 * half * n2, ctr),
                      np.ones(len(P), bool), 2))

    for k, (ttl, cam, sel, dec) in enumerate(panels):
        ax = fig.add_subplot(1, 2, k + 1, facecolor="white")
        ax.set_aspect("equal")
        ax.axis("off")
        # opaque planet: near-side graticule + coastlines + true horizon
        th = np.linspace(0, 2 * np.pi, 181)
        for lam in np.radians(np.arange(0, 180, 15)):
            ln = np.column_stack([np.cos(th) * np.cos(lam),
                                  np.cos(th) * np.sin(lam), np.sin(th)])
            draw_sphere_lines(ax, cam, ln, color="#c3c9cf", lw=0.3, zorder=1)
        for phi in np.radians(np.arange(-75, 76, 15)):
            ln = np.column_stack([np.cos(phi) * np.cos(th),
                                  np.cos(phi) * np.sin(th),
                                  np.full_like(th, np.sin(phi))])
            draw_sphere_lines(ax, cam, ln, color="#c3c9cf", lw=0.3, zorder=1)
        for c in coasts:
            draw_sphere_lines(ax, cam, c, color="#4a5560", lw=0.6, zorder=2)
        hx, hy, _ = cam.project(cam.horizon())
        ax.plot(hx, hy, color="#8a939c", lw=0.8, zorder=2)

        # trajectory: cull by line of sight, project, depth-sort far->near
        pts = P[sel][::dec]
        td = tdays[sel][::dec]
        occ = cam.occluded(pts)
        X, Y, z = cam.project(pts)
        good = ~occ[:-1] & ~occ[1:] & (z[:-1] > 0.05) & (z[1:] > 0.05)
        seg = np.stack([np.column_stack([X[:-1], Y[:-1]]),
                        np.column_stack([X[1:], Y[1:]])], axis=1)[good]
        zmid = (0.5 * (z[:-1] + z[1:]))[good]
        tmid = (0.5 * (td[:-1] + td[1:]))[good]
        order = np.argsort(-zmid)                      # paint far first
        zref = np.median(zmid)
        lw = np.clip(1.0 * zref / zmid[order], 0.45, 1.9)
        lc = LineCollection(seg[order], cmap="managua", norm=norm,
                            linewidths=lw, zorder=3)
        lc.set_array(tmid[order])
        ax.add_collection(lc)

        # markers (impact is on the near cap by camera construction)
        ix, iy, _ = cam.project(imp)
        ax.plot(ix, iy, marker="*", color="#d62728", ms=11, zorder=4)
        if k == 0:
            ax0x, ax0y, _ = cam.project(P[0])
            ax.plot(ax0x, ax0y, marker="o", color=INK, ms=3, zorder=4)
            ax.annotate("arrival", (ax0x[0], ax0y[0]), textcoords="offset points",
                        xytext=(6, -10), color=INK, fontsize=8)
        else:
            ax.annotate("impact", (ix[0], iy[0]), textcoords="offset points",
                        xytext=(8, 2), color="#d62728", fontsize=9)
        # frame the panel on what was drawn
        allx = np.concatenate([seg[:, :, 0].ravel(), hx[~np.isnan(hx)]])
        ally = np.concatenate([seg[:, :, 1].ravel(), hy[~np.isnan(hy)]])
        mx = 0.06 * max(allx.max() - allx.min(), ally.max() - ally.min())
        ax.set_xlim(allx.min() - mx, allx.max() + mx)
        ax.set_ylim(ally.min() - mx, ally.max() + mx)
        ax.set_title(ttl, color=INK, fontsize=10)

    fig.subplots_adjust(left=0.02, right=0.98, top=0.90, bottom=0.10, wspace=0.06)
    cax = fig.add_axes([0.41, 0.06, 0.18, 0.02])
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
