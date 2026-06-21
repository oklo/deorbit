#!/usr/bin/env python3
"""Kernel-faithful visualization of an SPH snapshot.

SPH fields are continuous: A(r) = sum_j m_j (A_j/rho_j) W(|r-r_j|, h_j). Instead
of drawing particles as points (which fakes a granularity that isn't there), we
RECONSTRUCT the field on a grid at the simulation's own resolution h and render
that:
  * deposit mass / mass*u / iron-mass onto a voxel grid,
  * blur by ~h (the cubic spline ~ a Gaussian of width h) -> smooth fields
    rho(r), u(r), iron_fraction(r),
  * photoreal: marching-cubes ISOSURFACE at rho=rho0/2 (the free surface) shaded
    by material, with shock-heated regions glowing,
  * tufte: a midplane SLICE of the internal-energy field with contours.

  ./recon.py snap.csv --dx 30 --downsample 2 --time 0.4 --out recon
"""
import argparse, csv
import numpy as np
from scipy.ndimage import gaussian_filter

RHO0 = {0: 2700.0, 1: 7800.0}      # basalt, iron
ETA = 1.3


def load(path):
    if path.endswith(".bin"):     # full-res binary: int64 n, then n*[x y z vx vy vz u rho h mat]
        with open(path, "rb") as f:
            n = int(np.fromfile(f, dtype="<i8", count=1)[0])
            a = np.fromfile(f, dtype="<f8", count=n * 10).reshape(n, 10)
        return dict(x=a[:, 0], y=a[:, 1], z=a[:, 2], u=a[:, 6], rho=a[:, 7],
                    h=a[:, 8], mat=a[:, 9].astype(int))
    cols = {k: [] for k in ("x", "y", "z", "u", "mat", "rho", "h")}   # CSV (baseline)
    rdr = csv.DictReader(open(path)); has = rdr.fieldnames
    for r in rdr:
        for k in ("x", "y", "z", "u"): cols[k].append(float(r[k]))
        cols["mat"].append(int(r["mat"]))
        cols["rho"].append(float(r["rho"]) if "rho" in has else np.nan)
        cols["h"].append(float(r["h"]) if "h" in has else np.nan)
    return {k: np.array(v) for k, v in cols.items()}


def reconstruct(d, dx, downsample, dz, pad=4):
    mat = d["mat"]; rho0 = np.array([RHO0[m] for m in mat])
    mass = rho0 * dx**3 * downsample            # scaled for the kept fraction
    h = d["h"].copy()
    if np.isnan(h).any(): h[:] = ETA * dx       # approximate if absent
    P = np.column_stack([d["x"], d["y"], d["z"]])
    lo = P.min(0) - pad * dz; hi = P.max(0) + pad * dz
    dims = np.ceil((hi - lo) / dz).astype(int)
    nx, ny, nz = dims
    idx = np.clip(((P - lo) / dz).astype(int), 0, dims - 1)
    flat = (idx[:, 0] * ny + idx[:, 1]) * nz + idx[:, 2]
    def grid(w):
        g = np.zeros(nx * ny * nz); np.add.at(g, flat, w); return g.reshape(nx, ny, nz)
    M = grid(mass); EU = grid(mass * d["u"]); MI = grid(mass * (mat == 1))
    sig = (ETA * dx) / dz * 0.6                  # blur ~ kernel width
    M = gaussian_filter(M, sig); EU = gaussian_filter(EU, sig); MI = gaussian_filter(MI, sig)
    rho = M / dz**3
    with np.errstate(divide="ignore", invalid="ignore"):
        u = np.where(M > 0, EU / M, 0.0)
        ironf = np.where(M > 0, MI / M, 0.0)
    return dict(rho=rho, u=u, ironf=ironf, lo=lo, dz=dz, dims=dims)


def photo(R, t, out):
    import pyvista as pv
    pv.OFF_SCREEN = True
    nx, ny, nz = R["dims"]
    g = pv.ImageData(dimensions=(nx, ny, nz), spacing=(R["dz"],) * 3, origin=tuple(R["lo"]))
    g.point_data["rho"] = R["rho"].flatten(order="F")
    g.point_data["u"] = R["u"].flatten(order="F")
    g.point_data["iron"] = R["ironf"].flatten(order="F")
    surf = g.contour([RHO0[0] / 2.0], scalars="rho")      # free surface
    surf = surf.smooth(n_iter=30)
    p = pv.Plotter(off_screen=True, window_size=(1700, 950))
    p.set_background("#0a0d12")
    # base material color: rock (warm grey) -> iron (steel), then glow by u
    irf = surf["iron"]; uu = surf["u"]; zpts = surf.points[:, 2]
    snow = np.clip((zpts - 5000.0) / 600.0, 0, 1)[:, None]   # glacier only on the high peak
    rocklo = np.array([0.34, 0.28, 0.22]); snowc = np.array([0.90, 0.92, 0.97])
    steel = np.array([0.40, 0.41, 0.43])        # neutral mid-grey iron (distinct from snow)
    terrain = rocklo[None] * (1 - snow) + snowc[None] * snow
    base = terrain * (1 - irf[:, None]) + steel[None] * irf[:, None]
    glow = np.clip(uu / 3e6, 0, 1)[:, None]                  # incandescent melt (lower clim -> pops)
    hot = np.array([1.0, 0.45, 0.10])
    rgb = np.clip(base * (1 - glow) + hot[None] * glow + glow * 0.8, 0, 1)
    surf["rgb"] = (rgb * 255).astype(np.uint8)
    p.add_mesh(surf, scalars="rgb", rgb=True, ambient=0.22, diffuse=0.9,
               specular=0.2, specular_power=12, smooth_shading=True)
    p.add_light(pv.Light(position=(-3, -4, 5), intensity=1.0))         # sun, windward
    p.add_light(pv.Light(position=(3, -2, 2), intensity=0.3, color="#9fb6ff"))  # sky
    b = surf.bounds
    sx = 0.5 * (b[0] + b[1])
    p.camera.focal_point = (sx * 0.3, 0, 0.7 * b[5])         # toward the summit/impact
    p.camera.position = (b[0] - 3000, b[2] - 8000, b[5] + 3500)   # windward, front, above
    p.camera.up = (0, 0, 1)
    p.camera.zoom(1.2)
    p.screenshot(out, scale=1); print("wrote", out)


def tufte(R, t, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    nx, ny, nz = R["dims"]
    j = ny // 2
    U = R["u"][:, j, :].T / 1e6                            # (nz,nx) MJ/kg
    rho = R["rho"][:, j, :].T
    x = (R["lo"][0] + np.arange(nx) * R["dz"]) / 1e3       # km
    z = (R["lo"][2] + np.arange(nz) * R["dz"]) / 1e3
    fig, ax = plt.subplots(figsize=(13, 5))
    im = ax.pcolormesh(x, z, U, cmap="magma", shading="auto", vmin=0, vmax=6)
    ax.contour(x, z, rho, levels=[RHO0[0] / 2], colors="white", linewidths=0.8)  # surface
    cb = fig.colorbar(im, ax=ax, pad=0.01); cb.set_label("specific internal energy  (MJ/kg)")
    ax.set_xlabel("downrange  (km)"); ax.set_ylabel("z  (km)")
    ax.set_title(f"reconstructed field (midplane)   t = {t:.2f} s")
    ax.set_aspect("equal"); ax.set_facecolor("#f7f7f7")
    fig.tight_layout(); fig.savefig(out, dpi=130); print("wrote", out)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("snap"); ap.add_argument("--dx", type=float, default=30.0)
    ap.add_argument("--downsample", type=int, default=1)   # 1 for full-res .bin snapshots
    ap.add_argument("--dz", type=float, default=18.0); ap.add_argument("--time", type=float, default=0.0)
    ap.add_argument("--out", default="recon"); ap.add_argument("--style", default="both")
    a = ap.parse_args()
    d = load(a.snap)
    print(f"{a.snap}: {len(d['x']):,} pts; reconstructing on dz={a.dz} m grid")
    R = reconstruct(d, a.dx, a.downsample, a.dz)
    print(f"grid {tuple(R['dims'])}  rho max {R['rho'].max():.0f}  u max {R['u'].max():.1e}")
    if a.style in ("photo", "both"): photo(R, a.time, a.out + "_photo.png")
    if a.style in ("tufte", "both"): tufte(R, a.time, a.out + "_tufte.png")
