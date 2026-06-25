#!/usr/bin/env python3
"""Photoreal animation: render each full-res snapshot via the kernel
reconstruction (isosurface, shaded) with a FIXED camera, stitch a GIF.

  ./anim_photo.py "ic_snap_0*.bin" --dx 30 --dz 50 --out hero_photo_anim.gif
"""
import argparse, glob, os
import recon
from PIL import Image

def vec(s): return tuple(float(v) for v in s.split(","))

ap = argparse.ArgumentParser()
ap.add_argument("pattern")
ap.add_argument("--dx", type=float, default=30.0)
ap.add_argument("--dz", type=float, default=50.0)
ap.add_argument("--dt", type=float, default=0.05)
ap.add_argument("--ms", type=int, default=400)
ap.add_argument("--tag", default="photo", help="frame-cache prefix (use a new tag for a new camera)")
ap.add_argument("--cam-focal", type=vec, default=(300.0, 0.0, 5300.0))
ap.add_argument("--cam-pos", type=vec, default=(-7000.0, -13000.0, 10500.0))
ap.add_argument("--zoom", type=float, default=1.2)
ap.add_argument("--out", default="hero_photo_anim.gif")
a = ap.parse_args()

snaps = sorted(glob.glob(a.pattern))
cam = (a.cam_focal, a.cam_pos)                    # fixed camera (focal, position)
os.makedirs("frames", exist_ok=True)
paths = []
for k, s in enumerate(snaps):
    out = f"frames/{a.tag}_{k:03d}.png"
    if os.path.exists(out):                       # cache: don't re-render existing frames
        paths.append(out); continue
    d = recon.load(s)
    R = recon.reconstruct(d, a.dx, 1, a.dz)
    recon.photo(R, a.dt * (k + 1), out, cam=cam, zoom=a.zoom)
    paths.append(out)
    print("rendered frame", k, s, flush=True)

imgs = [Image.open(p).convert("P", palette=Image.ADAPTIVE) for p in paths]
imgs[0].save(a.out, save_all=True, append_images=imgs[1:], duration=a.ms, loop=0)
print("wrote", a.out, f"({len(imgs)} frames)")
