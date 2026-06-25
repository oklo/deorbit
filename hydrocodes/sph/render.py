#!/usr/bin/env python3
"""Render an SPH impact snapshot two ways:
  --style tufte : maximal information / minimal ink. Clean scene, particles
                  colored by specific internal energy (the shock-heating /melt
                  field) on a perceptually-uniform map, labeled colorbar + scale.
  --style photo : naturalistic. Dark scene, shaded basalt, metallic iron, and an
                  incandescent black-body glow where the material is shock-heated.
Both read the snapshot CSV (x,y,z,vx,vy,vz,mat,u) written by run_ic/cayambe.

  ./render.py snap.csv --style both --time 0.40 --out frame
"""
import argparse, csv, math
import numpy as np
import pyvista as pv

pv.OFF_SCREEN = True


def load(path):
    if path.endswith('.bin'):     # int64 n, then n*[x y z vx vy vz u rho h mat]
        with open(path, 'rb') as f:
            n = int(np.fromfile(f, dtype='<i8', count=1)[0])
            a = np.fromfile(f, dtype='<f8', count=n * 10).reshape(n, 10)
        return a[:, :3], a[:, 6], a[:, 9].astype(int)
    x, y, z, u, mat = [], [], [], [], []
    for r in csv.DictReader(open(path)):
        x.append(float(r['x'])); y.append(float(r['y'])); z.append(float(r['z']))
        u.append(float(r['u'])); mat.append(int(r['mat']))
    return np.column_stack([x, y, z]), np.array(u), np.array(mat)


def frame_camera(p, P, mat):
    iron = P[mat == 1]
    c = iron.mean(0) if len(iron) else P.mean(0)
    span = max(2500.0, 3.0 * (np.ptp(iron, axis=0).max() if len(iron) else 2000))
    p.camera.focal_point = (c[0], 0, c[2])
    p.camera.position = (c[0], -span * 2.2, c[2] + span * 0.45)
    p.camera.up = (0, 0, 1)


HOT = 3.0e5   # J/kg: above this the material is "shock-heated" (the information)


def tufte(P, u, mat, t, out, slab=80.0):
    # thin midplane cross-section: a clean side view, no through-lattice moire
    m = np.abs(P[:, 1]) < slab
    Ps, us = P[m], u[m]
    cold = us < HOT; hot = ~cold
    p = pv.Plotter(off_screen=True, window_size=(1700, 950))
    p.set_background('white')
    if cold.any():                       # context: faint grey, low ink
        p.add_points(pv.PolyData(Ps[cold]), color='#c9c9c9', point_size=3.0,
                     render_points_as_spheres=False)
    if hot.any():                        # information: shock-heated material
        hc = pv.PolyData(Ps[hot]); hc['u_MJ'] = us[hot] / 1e6
        p.add_points(hc, scalars='u_MJ', cmap='plasma', clim=[HOT / 1e6, 8.0],
                     point_size=6, render_points_as_spheres=True,
                     scalar_bar_args=dict(title='internal energy (MJ/kg)', n_labels=5,
                                          color='black', vertical=True, position_x=0.9,
                                          title_font_size=20, label_font_size=16))
    xr = [P[:, 0].min(), P[:, 0].max()]
    p.add_mesh(pv.Line((xr[0], 0, 0), (xr[1], 0, 0)), color='black', line_width=1.5)
    p.view_xz()                          # look along -y: x right, z up
    p.camera.zoom(1.3)
    p.add_text(f"grazing iron impact into basalt  -  t = {t:.2f} s   (midplane slice)",
               color='black', font_size=14, position='upper_left')
    p.show_grid(color='black', xtitle='downrange x (m)', ztitle='z (m)', ytitle='')
    p.screenshot(out, scale=1)
    print("wrote", out)


def photo(P, u, mat, t, out):
    p = pv.Plotter(off_screen=True, window_size=(1700, 950))
    p.set_background('#0a0d12')
    rock, iron = mat == 0, mat == 1
    cold = u < HOT
    # cold basalt (shaded grey) and cold iron (steel), lit so relief shows
    rc = rock & cold; ic = iron & cold; hotm = ~cold
    if rc.any():
        p.add_points(pv.PolyData(P[rc]), color='#5b5550', point_size=4,
                     render_points_as_spheres=True, ambient=0.25, diffuse=0.8)
    if ic.any():
        p.add_points(pv.PolyData(P[ic]), color='#9aa0a8', point_size=5,
                     render_points_as_spheres=True, ambient=0.3, diffuse=0.7,
                     specular=0.6, specular_power=20)
    if hotm.any():                       # incandescent (black-body) glow
        hc = pv.PolyData(P[hotm]); hc['u'] = u[hotm]
        p.add_points(hc, scalars='u', cmap='afmhot', clim=[HOT, 1.5e7],
                     point_size=6, render_points_as_spheres=True, show_scalar_bar=False)
    p.add_light(pv.Light(position=(-2, -3, 4), intensity=1.0))
    p.add_light(pv.Light(position=(2, -2, 1), intensity=0.3, color='#9fb6ff'))
    frame_camera(p, P, mat)
    p.camera.position = (p.camera.position[0], p.camera.position[1] * 0.8,
                         p.camera.position[2] * 0.7)
    p.screenshot(out, scale=1)
    print("wrote", out)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("snap"); ap.add_argument("--style", default="both")
    ap.add_argument("--time", type=float, default=0.0); ap.add_argument("--out", default="frame")
    a = ap.parse_args()
    P, u, mat = load(a.snap)
    print(f"{a.snap}: {len(P):,} pts  iron={int((mat==1).sum()):,}  u range {u.min():.1e}..{u.max():.1e}")
    if a.style in ("tufte", "both"): tufte(P, u, mat, a.time, a.out + "_tufte.png")
    if a.style in ("photo", "both"): photo(P, u, mat, a.time, a.out + "_photo.png")
