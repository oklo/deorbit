#!/usr/bin/env python3
"""Build an SPH initial condition for the oblique Cayambe impact from real SRTM
topography. Output: a binary particle file consumed by run_ic.cpp.

Pipeline:
  1. Mosaic the two SRTM tiles (N00W079 windward + N00W078 summit), 30 m.
  2. Find the summit (max elevation); local ENU frame centred there
     (x = East = downrange, prograde body arrives from the west).
  3. Fill basalt particles below the real terrain surface (constant-thickness
     slab) over a cropped windward window; freeze the lateral + basal edges.
  4. Place the 1 km iron impactor up-azimuth at an altitude that grazes the
     summit, moving at the Phase-1 terminal conditions (7.4 km/s, ~2 deg).

Binary format (little-endian f8): xlo xhi ylo yhi zlo zhi n, then n*[x y z vx vy
vz mass u mat fixed].  mat: 0=basalt, 1=iron.
"""
import numpy as np
import argparse, sys

ap = argparse.ArgumentParser()
ap.add_argument("--dx", type=float, default=30.0)
ap.add_argument("--fill", type=float, default=2500.0, help="rock fill depth below surface (m)")
ap.add_argument("--xw", type=float, default=3500.0, help="upwind extent (m, west of summit)")
ap.add_argument("--xd", type=float, default=6500.0, help="downrange extent (m, east of summit)")
ap.add_argument("--yw", type=float, default=3000.0, help="half cross-range width (m)")
ap.add_argument("--U", type=float, default=7400.0)
ap.add_argument("--gamma", type=float, default=2.0, help="grazing angle (deg)")
ap.add_argument("--a", type=float, default=500.0, help="impactor radius (m)")
ap.add_argument("--shear", type=float, default=400.0,
                help="depth below summit the impactor CENTER passes (m); it shears the top ~shear+a")
ap.add_argument("--out", default="cayambe_ic.bin")
A = ap.parse_args()

RHO_B, RHO_I = 2700.0, 7800.0
dx = A.dx; h2 = 2.0 * 1.3 * dx

# --- mosaic the two tiles: lon -79..-77, lat 0..1, 30 m ---
def load(t):
    d = np.fromfile(f"dem/{t}.hgt", dtype=">i2").astype(np.float64).reshape(3601, 3601)
    d[d < -1000] = np.nan
    return d
W079, W078 = load("N00W079"), load("N00W078")
COMB = np.hstack([W079[:, :3600], W078])          # cols: lon = -79 + col/3600
NC = COMB.shape[1]                                # 7201
def lonlat_to_rc(lon, lat):
    return (1.0 - lat) * 3600.0, (lon + 79.0) * 3600.0   # row, col (floats)

# summit = global max over the Cayambe neighbourhood (lon -78.1..-77.9, lat -0.05..0.15)
r0, r1 = int((1 - 0.15) * 3600), int((1 - (-0.05)) * 3600)
c0, c1 = int((-78.1 + 79) * 3600), int((-77.9 + 79) * 3600)
sub = COMB[r0:r1, c0:c1]
ij = np.unravel_index(np.nanargmax(sub), sub.shape)
srow, scol = r0 + ij[0], c0 + ij[1]
lat_s = 1.0 - srow / 3600.0
lon_s = -79.0 + scol / 3600.0
z_s = COMB[srow, scol]
print(f"summit: lat={lat_s:.4f} lon={lon_s:.4f} elev={z_s:.0f} m")

MPD_LAT = 111320.0
MPD_LON = 111320.0 * np.cos(np.radians(lat_s))

def terrain(x, y):                                # bilinear elevation at local (x,y) m
    lat = lat_s + y / MPD_LAT
    lon = lon_s + x / MPD_LON
    fr, fc = lonlat_to_rc(lon, lat)
    r = np.clip(fr.astype(int), 0, 3599); c = np.clip(fc.astype(int), 0, NC - 2)
    tr = fr - r; tc = fc - c
    z = (COMB[r, c] * (1 - tr) * (1 - tc) + COMB[r, c + 1] * (1 - tr) * tc
         + COMB[r + 1, c] * tr * (1 - tc) + COMB[r + 1, c + 1] * tr * tc)
    return z

# --- basalt fill over the windward window (summit at local origin) ---
xs = np.arange(-A.xw + dx / 2, A.xd, dx)
ys = np.arange(-A.yw + dx / 2, A.yw, dx)
X, Y = np.meshgrid(xs, ys, indexing="ij")
ZT = terrain(X.ravel(), Y.ravel()).reshape(X.shape)
nlev = int(A.fill / dx)
zbase_min = np.nanmin(ZT) - A.fill
parts = []
for k in range(nlev):                             # k=0 just below surface, downward
    z = ZT - (k + 0.5) * dx
    fx = (np.abs(X - (-A.xw)) < h2) | (np.abs(X - A.xd) < h2)
    fy = np.abs(np.abs(Y) - A.yw) < h2
    fb = (k >= nlev - 2)
    fixed = (fx | fy | fb).astype(np.float64)
    n = X.size
    parts.append(np.column_stack([
        X.ravel(), Y.ravel(), z.ravel(),
        np.zeros(n), np.zeros(n), np.zeros(n),
        np.full(n, RHO_B * dx**3), np.zeros(n), np.zeros(n), fixed.ravel()]))
basalt = np.vstack(parts)
print(f"basalt particles: {basalt.shape[0]:,} ({len(xs)}x{len(ys)}x{nlev})")

# --- iron impactor: bottom grazes the summit at x=0 ---
g = np.radians(A.gamma); vh = A.U * np.cos(g); vn = A.U * np.sin(g)
xb0 = -A.xw + A.a + dx                              # start just inside the upwind edge
# aim the impactor CENTRE to pass A.shear below the summit at x=0, so it ploughs
# through the upper cone (shears the top ~shear+a) instead of kissing the apex
zb0 = (z_s - A.shear) + np.tan(g) * (-xb0)
na = int(np.ceil(A.a / dx)) + 1
ix, iy, iz = np.meshgrid(*(np.arange(-na, na + 1),) * 3, indexing="ij")
xp = ix * dx; yp = iy * dx; zp = zb0 + iz * dx
m = (xp**2 + yp**2 + (zp - zb0)**2) <= A.a**2
xp, yp, zp = xp[m] + xb0, yp[m], zp[m]
ni = xp.size
iron = np.column_stack([xp, yp, zp,
                        np.full(ni, vh), np.zeros(ni), np.full(ni, -vn),
                        np.full(ni, RHO_I * dx**3), np.zeros(ni),
                        np.ones(ni), np.zeros(ni)])
print(f"iron particles: {ni:,}  start=({xb0:.0f},0,{zb0:.0f}) v=({vh:.0f},0,{-vn:.0f})")

data = np.vstack([basalt, iron])
N = data.shape[0]
xlo, xhi = -A.xw - dx, A.xd + dx
ylo, yhi = -A.yw - dx, A.yw + dx
zlo, zhi = zbase_min - dx, zb0 + A.a + 600.0
hdr = np.array([xlo, xhi, ylo, yhi, zlo, zhi, float(N)], dtype="<f8")
with open(A.out, "wb") as f:
    hdr.tofile(f); data.astype("<f8").tofile(f)
print(f"wrote {A.out}: N={N:,}  domain x[{xlo:.0f},{xhi:.0f}] y[{ylo:.0f},{yhi:.0f}] z[{zlo:.0f},{zhi:.0f}]")
print(f"  (~{N * 80 / 1e6:.0f} MB; CPPR = {A.a/dx:.1f})")
