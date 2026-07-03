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


def draw_sphere_lines(ax, cam, poly, R=1.0, extra_mask=None, **kw):
    X, Y, z = cam.project(poly)
    vis = cam.near_cap(poly, R) & (z > 1e-3)
    if extra_mask is not None:
        vis = vis & extra_mask
    X = np.where(vis, X, np.nan)
    Y = np.where(vis, Y, np.nan)
    ax.plot(X, Y, **kw)
    return X, Y


def _circular_arc(pts, mask):
    """Extract the (single) contiguous True run from a cyclic mask, returned
    in traversal order."""
    n = len(mask)
    if mask.all():
        return pts
    start = next(i for i in range(n) if mask[i] and not mask[i - 1])
    idx = [(start + k) % n for k in range(n) if mask[(start + k) % n]]
    return pts[idx]


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
    plane_view(P, tdays, norm, geo, coasts, out, d)


def plane_view(P, tdays, norm, geo, coasts, out, draw):
    """Vantage tipped -45 deg from the capture-plane normal (perigee side
    toward the viewer). All translucent surfaces (orbital-plane sheet,
    equatorial rotation arrow, orbit-following motion arrow) are true 3D
    planar objects tessellated into cells; each cell is line-of-sight
    culled against the globe, so occlusion is exact everywhere with no
    cover polygons. Unlabeled (caption carries the annotations); the
    apsidal table prints to stdout."""
    r = np.linalg.norm(P, axis=1)
    iapo = [i for i in range(1, len(r) - 1)
            if r[i] > 2.0 and r[i] >= r[i - 1] and r[i] > r[i + 1]]
    iper = [i for i in range(1, len(r) - 1)
            if r[i] < 1.1 and r[i] <= r[i - 1] and r[i] < r[i + 1]]
    hhat = np.asarray(geo["h_hat"])
    e1 = P[iapo[0]] - (P[iapo[0]] @ hhat) * hhat
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(hhat, e1)
    print("apogee-direction angle in the initial orbit plane (deg), per pass:")
    angs = []
    for k, i in enumerate(iapo):
        a = math.degrees(math.atan2(P[i] @ e2, P[i] @ e1))
        hk = np.cross(P[max(0, i - 40)], P[min(len(P) - 1, i + 40)])
        hk /= np.linalg.norm(hk)
        tilt = math.degrees(math.acos(abs(float(np.clip(hk @ hhat, -1, 1)))))
        angs.append(a)
        print(f"  pass {k+1}: apo_dir {a:+7.2f} deg   plane tilt vs initial "
              f"{tilt:5.2f} deg   r_apo {r[i]:.2f} R_E")
    print(f"apogee-direction spread {max(angs)-min(angs):.2f} deg; "
          f"cascade {tdays[-1]:.2f} d = {tdays[-1]*86400/86164.1:.2f} sidereal rev")

    # tallest (first) loop in-plane extents: sheet = 1.1x its semi-axes,
    # centered on the ELLIPSE (focus offset ~5 R_E along the apsidal axis)
    i0l, i1l = iper[0], (iper[1] if len(iper) > 1 else len(P) - 1)
    uu = P[i0l:i1l] @ e1
    vv = P[i0l:i1l] @ e2
    uc, vc = 0.5 * (uu.min() + uu.max()), 0.5 * (vv.min() + vv.max())
    hu = 1.1 * 0.5 * (uu.max() - uu.min())
    hv = 1.1 * 0.5 * (vv.max() - vv.min())
    print(f"sheet: center ({uc:+.2f},{vc:+.2f}) R_E in-plane, half-extents "
          f"{hu:.2f} x {hv:.2f} R_E (1.1x loop-1 semi-axes)")

    TILT = math.radians(45.0)
    EYE_DIST = 7.5
    n_hat = math.cos(TILT) * hhat - math.sin(TILT) * e1
    cam = Camera(EYE_DIST * n_hat, (0.0, 0.0, 0.0), up=e1)

    fig = plt.figure(figsize=(6.9, 13.8), facecolor="white")
    ax = fig.add_subplot(111, facecolor="white")
    ax.set_aspect("equal")
    ax.axis("off")
    th = np.linspace(0, 2 * np.pi, 361)

    # ---- helpers: tessellated translucent surfaces, cell-level occlusion ----
    def cells_visible(centers):
        _, _, zc = cam.project(centers)
        return (~cam.occluded(centers)) & (zc > 0.05)

    def fill_quads(quads, color, alpha, zorder):
        from matplotlib.collections import PolyCollection
        cen = quads.mean(axis=1)
        ok = cells_visible(cen)
        q = quads[ok]
        if not len(q):
            return
        flat = q.reshape(-1, 3)
        X, Y, _ = cam.project(flat)
        scr = np.stack([X.reshape(-1, 4).T, Y.reshape(-1, 4).T]).T
        pc = PolyCollection(list(scr), facecolor=color, alpha=alpha,
                            edgecolor=color, linewidths=0.3, zorder=zorder)
        ax.add_collection(pc)

    def fill_strip(outer, inner, color, alpha, zorder):
        """Planar strip: merge contiguous visible spans into single fills
        (seamless at low alpha)."""
        cen = 0.5 * (outer + inner)
        ok = cells_visible(cen)
        k0 = None
        for k in range(len(ok) + 1):
            on = k < len(ok) and ok[k]
            if on and k0 is None:
                k0 = k
            if not on and k0 is not None:
                if k - k0 > 1:
                    ox, oy, _ = cam.project(outer[k0:k])
                    ixn, iyn, _ = cam.project(inner[k0:k][::-1])
                    ax.fill(np.concatenate([ox, ixn]), np.concatenate([oy, iyn]),
                            color=color, alpha=alpha, lw=0, zorder=zorder)
                k0 = None

    def fill_head(tip, base_c, side_hat, half_w, color, alpha, zorder):
        tri = np.array([tip, base_c + half_w * side_hat, base_c - half_w * side_hat])
        if cells_visible(tri.mean(axis=0)[None, :])[0]:
            tx, ty, _ = cam.project(tri)
            ax.fill(tx, ty, color=color, alpha=alpha, lw=0, zorder=zorder)

    # ---- globe: near-side graticule + coastlines + horizon ----
    for lam in np.radians(np.arange(0, 180, 15)):
        ln = np.column_stack([np.cos(th) * np.cos(lam),
                              np.cos(th) * np.sin(lam), np.sin(th)])
        draw_sphere_lines(ax, cam, ln, color="#c3c9cf", lw=0.3, zorder=4)
    for phi in np.radians(np.arange(-75, 76, 15)):
        ln = np.column_stack([np.cos(phi) * np.cos(th),
                              np.cos(phi) * np.sin(th),
                              np.full_like(th, np.sin(phi))])
        draw_sphere_lines(ax, cam, ln, color="#c3c9cf", lw=0.3, zorder=4)
    for c in coasts:
        draw_sphere_lines(ax, cam, c, color="#4a5560", lw=0.55, zorder=4)
    H = cam.horizon(n=361)
    hx, hy, _ = cam.project(H)
    ax.plot(hx, hy, color="#8a939c", lw=0.8, zorder=4)
    # in-plane unit circle (plane-planet intersection), visible arc
    circ3 = np.outer(np.cos(th), e1) + np.outer(np.sin(th), e2)
    draw_sphere_lines(ax, cam, circ3, color="#aab4bd", lw=0.5, zorder=4)

    # ---- orbital-plane sheet: rectangle 1.1x the tallest loop's semi-axes,
    # centered on the ellipse (cell-tessellated, cell ~0.075 R_E) ----
    nu = max(2, int(round(2 * hu / 0.075)))
    nv = max(2, int(round(2 * hv / 0.075)))
    gu = np.linspace(uc - hu, uc + hu, nu + 1)
    gv = np.linspace(vc - hv, vc + hv, nv + 1)
    GA, GB = np.meshgrid(gu[:-1], gv[:-1], indexing="ij")
    dua, dvb = gu[1] - gu[0], gv[1] - gv[0]
    corners = []
    for du, dv in ((0, 0), (dua, 0), (dua, dvb), (0, dvb)):
        A = (GA + du).ravel()[:, None]
        B = (GB + dv).ravel()[:, None]
        corners.append(A * e1 + B * e2)
    sheet_quads = np.stack(corners, axis=1)
    fill_quads(sheet_quads, "#7d94aa", 0.10, 2)
    # crisp rectangle outline, LOS-culled
    tu = np.linspace(uc - hu, uc + hu, 121)
    tv = np.linspace(vc - hv, vc + hv, 61)
    per = np.concatenate([
        np.column_stack([tu, np.full_like(tu, vc - hv)]),
        np.column_stack([np.full_like(tv, uc + hu), tv]),
        np.column_stack([tu[::-1], np.full_like(tu, vc + hv)]),
        np.column_stack([np.full_like(tv, uc - hu), tv[::-1]])])
    sq3 = per[:, :1] * e1 + per[:, 1:2] * e2
    occ_sq = cam.occluded(sq3)
    sqx, sqy, sqz = cam.project(sq3)
    bad = occ_sq | (sqz <= 0.05)
    ax.plot(np.where(~bad, sqx, np.nan), np.where(~bad, sqy, np.nan),
            color="#93a7bb", lw=0.8, zorder=2)

    # ---- trajectory, with the incoming hyperbola reconstructed ----
    from deorbit.corridor.physics import MU
    s0, _ = sample.entry_state(draw["v_inf"], draw["u_hat"], draw["theta_B"],
                               A_EQ + draw["rp_alt_m"])
    rv, vv = np.array(s0[:3]), np.array(s0[3:])
    hv = np.cross(rv, vv)
    ev = np.cross(vv, hv) / MU - rv / np.linalg.norm(rv)
    ecc = np.linalg.norm(ev)
    pfoc = (hv @ hv) / MU
    p_hat = ev / ecc
    q_hat = np.cross(hv / np.linalg.norm(hv), p_hat)
    r0n = np.linalg.norm(rv)
    nu0 = -math.acos(max(-1.0, min(1.0, (pfoc / r0n - 1.0) / ecc)))
    r_far = 16.0 * A_EQ
    nu_far = -math.acos(max(-1.0, min(1.0, (pfoc / r_far - 1.0) / ecc)))
    nus = np.linspace(nu_far, nu0, 250)
    rr = pfoc / (1.0 + ecc * np.cos(nus))
    hyp = (rr[:, None] * (np.cos(nus)[:, None] * p_hat
                          + np.sin(nus)[:, None] * q_hat)) / A_EQ
    Pd = np.vstack([hyp, P])
    tdd = np.concatenate([np.full(len(hyp), -1e-4), tdays])

    occ = cam.occluded(Pd)
    X, Y, z = cam.project(Pd)
    good = ~occ[:-1] & ~occ[1:] & (z[:-1] > 0.05) & (z[1:] > 0.05)
    seg = np.stack([np.column_stack([X[:-1], Y[:-1]]),
                    np.column_stack([X[1:], Y[1:]])], axis=1)[good]
    zmid = (0.5 * (z[:-1] + z[1:]))[good]
    tmid = (0.5 * (tdd[:-1] + tdd[1:]))[good]
    order = np.argsort(-zmid)
    zref = np.median(zmid)
    lw = np.clip(0.95 * zref / zmid[order], 0.5, 1.8)
    lc = LineCollection(seg[order], cmap="managua", norm=norm,
                        linewidths=lw, zorder=5)
    lc.set_array(tmid[order])
    ax.add_collection(lc)

    # ---- apsidal lines (perigee through apogee), impact star ----
    cmap = plt.get_cmap("managua")

    def plane_line(dirv, color, t0, t1):
        ray = np.linspace(t0, t1, 500)[:, None] * dirv[None, :]
        rocc = cam.occluded(ray)
        rin = np.linalg.norm(ray, axis=1) < 0.999
        rx, ry, rz = cam.project(ray)
        rbad = rocc | rin | (rz <= 0.05)
        ax.plot(np.where(~rbad, rx, np.nan), np.where(~rbad, ry, np.nan),
                ls="--", lw=0.7, color=color, alpha=0.85, zorder=5)

    EXTRUDE = 0.8                              # beyond the sheet edges
    for i in (iapo[0], iapo[-1]):              # initial + final apsidal lines
        plane_line(P[i] / r[i], cmap(norm(tdays[i])),
                   uc - hu - EXTRUDE, uc + hu + EXTRUDE)
    plane_line(e2, "#9aa2ab", vc - hv - EXTRUDE, vc + hv + EXTRUDE)
    if cam.near_cap(P[-1])[0]:
        ix, iy, _ = cam.project(P[-1])
        ax.plot(ix, iy, marker="*", color="#d62728", ms=12, zorder=6)

    # ---- Earth-rotation arrow: annular ribbon IN the equatorial plane ----
    eq_dir = cam.E - cam.E[2] * np.array([0.0, 0.0, 1.0])
    phi0 = math.atan2(eq_dir[1], eq_dir[0])
    phe = phi0 + np.radians(np.linspace(-28.0, 2.0, 40))
    rin_, rout_ = 1.05, 1.21
    ring = np.column_stack([np.cos(phe), np.sin(phe), np.zeros_like(phe)])
    fill_strip(rout_ * ring, rin_ * ring, "#c08a3e", 0.30, 4.5)
    pm = phe[-1]
    rmid = 0.5 * (rin_ + rout_)
    t_hat = np.array([-math.sin(pm), math.cos(pm), 0.0])
    r_hat = np.array([math.cos(pm), math.sin(pm), 0.0])
    fill_head(rmid * r_hat + 0.26 * t_hat, rmid * r_hat, r_hat, 0.145,
              "#c08a3e", 0.30, 4.5)

    # ---- incoming-motion arrow: ON the hyperbolic approach leg ----
    rhyp = np.linalg.norm(hyp, axis=1)
    hsel = (rhyp > 2.7) & (rhyp < 3.6)          # just outside the sheet, left
    Q = hyp[hsel][::2]
    T = np.gradient(Q, axis=0)
    T /= np.linalg.norm(T, axis=1)[:, None]
    M = np.cross(T, hhat)
    M /= np.linalg.norm(M, axis=1)[:, None]
    hw = 0.055
    fill_strip(Q + hw * M, Q - hw * M, "#3b6ea5", 0.35, 4.5)
    fill_head(Q[-1] + 0.20 * T[-1], Q[-1], M[-1], 0.10, "#3b6ea5", 0.35, 4.5)

    # ---- frame on the globe neighborhood ----
    # frame: comfortably enclose everything drawn (loops, sheet, globe)
    allx = np.concatenate([seg[:, :, 0].ravel(), sqx[~np.isnan(sqx)], hx])
    ally = np.concatenate([seg[:, :, 1].ravel(), sqy[~np.isnan(sqy)], hy])
    mx = 0.035 * max(allx.max() - allx.min(), ally.max() - ally.min())
    x0, x1 = allx.min() - mx, allx.max() + mx
    y0, y1 = ally.min() - mx, ally.max() + mx
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    asp = (y1 - y0) / (x1 - x0)
    W = 12.5 if asp <= 1.0 else 12.5 / asp
    fig.set_size_inches(W, W * asp)
    fig.subplots_adjust(left=0.005, right=0.995, top=0.995, bottom=0.005)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(HERE, "figures", f"{out}_plane.{ext}"), dpi=180,
                    facecolor="white", bbox_inches="tight")
    print(f"wrote figures/{out}_plane.png/.pdf")


if __name__ == "__main__":
    main()
