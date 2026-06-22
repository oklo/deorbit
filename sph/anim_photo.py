#!/usr/bin/env python3
"""Photoreal animation: render each full-res snapshot via the kernel
reconstruction (isosurface, shaded) with a FIXED camera, stitch a GIF.

  ./anim_photo.py "ic_snap_0*.bin" --dx 30 --dz 50 --out hero_photo_anim.gif
"""
import argparse, glob, os
import recon
from PIL import Image

ap = argparse.ArgumentParser()
ap.add_argument("pattern")
ap.add_argument("--dx", type=float, default=30.0)
ap.add_argument("--dz", type=float, default=50.0)
ap.add_argument("--dt", type=float, default=0.05)
ap.add_argument("--ms", type=int, default=400)
ap.add_argument("--out", default="hero_photo_anim.gif")
a = ap.parse_args()

snaps = sorted(glob.glob(a.pattern))
# fixed camera (windward, front, above; framed for the cone + rising plume)
cam = ((300.0, 0.0, 5300.0), (-7000.0, -13000.0, 10500.0))
os.makedirs("frames", exist_ok=True)
paths = []
for k, s in enumerate(snaps):
    d = recon.load(s)
    R = recon.reconstruct(d, a.dx, 1, a.dz)
    out = f"frames/photo_{k:03d}.png"
    recon.photo(R, a.dt * (k + 1), out, cam=cam)
    paths.append(out)
    print("frame", k, s, flush=True)

imgs = [Image.open(p).convert("P", palette=Image.ADAPTIVE) for p in paths]
imgs[0].save(a.out, save_all=True, append_images=imgs[1:], duration=a.ms, loop=0)
print("wrote", a.out, f"({len(imgs)} frames)")
