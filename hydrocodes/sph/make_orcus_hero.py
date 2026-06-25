#!/usr/bin/env python3
"""Cayambe-style 3D photoreal hero animation for the Orcus impact.
Reuses recon.reconstruct; Orcus-specific photo: vertical exaggeration (shallow
basin), hot-ejecta incandescent glow, NO 'snow' term, oblique-aerial camera.
  test one frame:  ../.venv/bin/python make_orcus_hero.py --test ../gpu/gpu_snap_010.bin
  full animation:  ../.venv/bin/python make_orcus_hero.py
"""
import numpy as np, glob, os, subprocess, sys
import recon
from PIL import Image
import imageio_ffmpeg

EXAG = 1.0          # NO vertical exaggeration -- accuracy over drama
DZ   = 1500.0       # finer grid to resolve the true-scale (shallow) crater
SUB  = 2            # 1-in-2 (the crater-region crop already reduces the count)
EXTRA_BLUR = 1.5
# camera: 100 km altitude, 150 km off the midline, 300 km uprange of impact
# (crater at x ~ -165 km), looking downrange (+x). [alt -50; horiz pulled back 1.5x]
POS   = (-465000.0, -150000.0, 100000.0)
FOCAL = (-50000.0, 0.0, -2000.0)
ZOOM  = 1.0

def photo_orcus(R, out):
    import pyvista as pv; pv.OFF_SCREEN = True
    nx,ny,nz = R["dims"]
    g = pv.ImageData(dimensions=(nx,ny,nz), spacing=(R["dz"],)*3, origin=tuple(R["lo"]))
    g.point_data["rho"]=R["rho"].flatten(order="F"); g.point_data["u"]=R["u"].flatten(order="F")
    surf = g.contour([recon.RHO0[0]/2.0], scalars="rho").smooth(n_iter=20)
    pts = surf.points.copy(); zmean = np.median(pts[:,2])     # exaggerate RELIEF (not density)
    pts[:,2] = zmean + (pts[:,2]-zmean)*EXAG; surf.points = pts
    uu = surf["u"]
    rock = np.array([0.34,0.28,0.22])
    glow = np.clip(uu/3e6,0,1)[:,None]; hot = np.array([1.0,0.45,0.10])
    rgb = np.clip(rock[None]*(1-glow) + hot[None]*glow + glow*0.8, 0, 1)
    surf["rgb"] = (rgb*255).astype(np.uint8)
    p = pv.Plotter(off_screen=True, window_size=(1700,950)); p.set_background("#0a0d12")
    p.add_mesh(surf, scalars="rgb", rgb=True, ambient=0.25, diffuse=0.9, specular=0.2,
               specular_power=12, smooth_shading=True)
    p.add_light(pv.Light(position=(-3,-4,5), intensity=1.0))
    p.add_light(pv.Light(position=(3,-2,2), intensity=0.3, color="#9fb6ff"))
    p.camera.focal_point=FOCAL; p.camera.position=POS; p.camera.up=(0,0,1); p.camera.zoom(ZOOM)
    p.screenshot(out, scale=1); print("wrote", out, flush=True)

def render(snap, out):
    d = recon.load(snap)
    m = (d["x"]>-480e3)&(d["x"]<200e3)&(np.abs(d["y"])<180e3)   # crop to crater region (slab is huge)
    for k in d: d[k] = d[k][m]
    n = len(d["x"]); keep = np.arange(0, n, SUB)
    for k in d: d[k] = d[k][keep]
    d["mat"][d["mat"]==2] = 0          # impactor is basalt
    R = recon.reconstruct(d, 1000.0, SUB, DZ)          # true z -> correct density (~rho0)
    from scipy.ndimage import gaussian_filter
    R["rho"] = gaussian_filter(R["rho"], EXTRA_BLUR)   # smooth the surface (subsampled + jittery)
    photo_orcus(R, out)

if __name__ == "__main__":
    if len(sys.argv)>2 and sys.argv[1]=="--test":
        render(sys.argv[2], "orcus_hero_test.png"); sys.exit()
    snaps = sorted(glob.glob("../gpu/gpu_snap_*.bin"))
    snaps = [s for s in snaps if os.path.getsize(s) > 3.0e9][:-1]   # complete frames only
    os.makedirs("hframes", exist_ok=True); paths=[]
    for k,s in enumerate(snaps):
        out=f"hframes/h_{k:03d}.png"
        if not os.path.exists(out): render(s, out)
        paths.append(out)
    imgs=[Image.open(p).convert("P",palette=Image.ADAPTIVE) for p in paths]
    imgs[0].save("orcus_hero.gif", save_all=True, append_images=imgs[1:], duration=200, loop=0)
    ff=imageio_ffmpeg.get_ffmpeg_exe()
    subprocess.run([ff,"-y","-framerate","6","-i","hframes/h_%03d.png","-vf",
        "tpad=stop_mode=clone:stop_duration=1.5,scale=trunc(iw/2)*2:trunc(ih/2)*2",
        "-c:v","libx264","-pix_fmt","yuv420p","-crf","20","-movflags","+faststart","orcus_hero.mp4"],check=True)
    print(f"wrote orcus_hero.gif + orcus_hero.mp4 ({len(paths)} frames)")
