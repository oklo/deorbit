#!/usr/bin/env python3
"""Iron-only development animation: reconstruct the iron mass-density field
(kernel-smoothed by h -> space-filling), isosurface it, follow the iron centroid
with a camera that dollies back as the blob stretches. Output GIF + MP4.

  ./anim_iron.py "frac_snap_0*.bin" --dx 30 --dz 40 --out iron_dev
"""
import argparse, glob, os, subprocess
import recon
from PIL import Image
import imageio_ffmpeg

ap = argparse.ArgumentParser()
ap.add_argument("pattern"); ap.add_argument("--dx", type=float, default=30.0)
ap.add_argument("--dz", type=float, default=40.0); ap.add_argument("--dt", type=float, default=0.05)
ap.add_argument("--fps", type=int, default=6); ap.add_argument("--out", default="iron_dev")
a = ap.parse_args()

snaps = sorted(glob.glob(a.pattern)); os.makedirs("frames", exist_ok=True)
paths = []
for k, s in enumerate(snaps):
    out = f"frames/iron_{k:03d}.png"
    if not os.path.exists(out):
        d = recon.load(s)
        R = recon.reconstruct_iron(d, a.dx, a.dz)
        recon.photo_iron(R, a.dt * (k + 1), out)
        print("rendered", k, s, flush=True)
    paths.append(out)

imgs = [Image.open(p).convert("P", palette=Image.ADAPTIVE) for p in paths]
imgs[0].save(a.out + ".gif", save_all=True, append_images=imgs[1:], duration=int(1000/a.fps), loop=0)
ff = imageio_ffmpeg.get_ffmpeg_exe()
subprocess.run([ff, "-y", "-framerate", str(a.fps), "-i", "frames/iron_%03d.png",
    "-vf", "tpad=stop_mode=clone:stop_duration=1.5,scale=trunc(iw/2)*2:trunc(ih/2)*2",
    "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "20", "-movflags", "+faststart",
    a.out + ".mp4"], check=True)
print(f"wrote {a.out}.gif + {a.out}.mp4 ({len(paths)} frames)")
