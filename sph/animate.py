#!/usr/bin/env python3
"""Stitch a Tufte animation from a series of SPH snapshots (kernel-reconstructed
field, fixed axes/colour scale across frames). Output: an animated GIF (PIL).

  ./animate.py "cayambe_snap_0*.csv" --dx 30 --downsample 8 --dt 0.1 --out baseline_anim.gif
"""
import argparse, glob, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image
import recon

ap = argparse.ArgumentParser()
ap.add_argument("pattern")
ap.add_argument("--dx", type=float, default=30.0)
ap.add_argument("--downsample", type=int, default=8)
ap.add_argument("--dz", type=float, default=18.0)
ap.add_argument("--dt", type=float, default=0.1, help="sim time between frames (s)")
ap.add_argument("--xlim", default="-1.5,8.5"); ap.add_argument("--zlim", default="-1.6,1.0")
ap.add_argument("--clim", type=float, default=6.0)
ap.add_argument("--ms", type=int, default=500, help="GIF frame duration (ms)")
ap.add_argument("--out", default="baseline_anim.gif")
a = ap.parse_args()

snaps = sorted(glob.glob(a.pattern))
xlim = tuple(float(v) for v in a.xlim.split(","))
zlim = tuple(float(v) for v in a.zlim.split(","))
os.makedirs("frames", exist_ok=True)
paths = []
for k, s in enumerate(snaps):
    d = recon.load(s)
    R = recon.reconstruct(d, a.dx, a.downsample, a.dz)
    nx, ny, nz = R["dims"]; j = ny // 2
    U = R["u"][:, j, :].T / 1e6; rho = R["rho"][:, j, :].T
    x = (R["lo"][0] + np.arange(nx) * R["dz"]) / 1e3
    z = (R["lo"][2] + np.arange(nz) * R["dz"]) / 1e3
    t = a.dt * (k + 1)
    fig, ax = plt.subplots(figsize=(13, 4.2))
    im = ax.pcolormesh(x, z, U, cmap="magma", vmin=0, vmax=a.clim, shading="auto")
    ax.contour(x, z, rho, levels=[recon.RHO0[0] / 2], colors="white", linewidths=0.7)
    ax.set_xlim(xlim); ax.set_ylim(zlim); ax.set_aspect("equal")
    ax.set_facecolor("black")
    cb = fig.colorbar(im, ax=ax, pad=0.01); cb.set_label("specific internal energy  (MJ/kg)")
    ax.set_xlabel("downrange  (km)"); ax.set_ylabel("z  (km)")
    ax.set_title(f"grazing iron impact into basalt  -  t = {t:.2f} s")
    fig.tight_layout()
    fp = f"frames/anim_{k:03d}.png"; fig.savefig(fp, dpi=110); plt.close(fig)
    paths.append(fp); print("frame", k, t)

imgs = [Image.open(p).convert("P", palette=Image.ADAPTIVE) for p in paths]
imgs[0].save(a.out, save_all=True, append_images=imgs[1:], duration=a.ms, loop=0)
print("wrote", a.out, f"({len(imgs)} frames)")
