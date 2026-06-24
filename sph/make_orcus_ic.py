#!/usr/bin/env python3
"""SPH initial condition for an ultra-oblique asteroid impact on Mars
(Orcus Patera test). A flat basalt half-space (surface at z=0) pre-stressed to
lithostatic equilibrium via a thermal internal-energy profile u(z)=g|z|/(a+b)
[basalt Tillotson a+b=2.0], so it does NOT slump when gravity is applied at
runtime (gpu_ic --g). A basalt impactor grazes at a shallow angle, heading +x
(the +x axis = the NNE-SSW long axis / downrange direction).

Same binary format as make_cayambe_ic.py (little-endian f8):
  header [xlo xhi ylo yhi zlo zhi n], then n*[x y z vx vy vz mass u mat fixed].
  mat 0 = basalt.
Run e.g.:  ./gpu_ic orcus_ic.bin 3000 200 21600 --g 3.71 --snap 5
"""
import numpy as np, argparse

ap = argparse.ArgumentParser()
ap.add_argument("--dx", type=float, default=3000.0)
ap.add_argument("--D", type=float, default=40e3, help="impactor diameter (m)")
ap.add_argument("--U", type=float, default=10000.0, help="impact speed (m/s)")
ap.add_argument("--theta", type=float, default=5.0, help="impact angle below horizontal (deg)")
ap.add_argument("--g", type=float, default=3.71, help="Mars gravity (m/s^2; must match gpu_ic --g)")
ap.add_argument("--xlo", type=float, default=-400e3); ap.add_argument("--xhi", type=float, default=450e3)
ap.add_argument("--yhalf", type=float, default=200e3)
ap.add_argument("--depth", type=float, default=120e3, help="slab depth (m)")
ap.add_argument("--x0", type=float, default=-250e3, help="impactor start x center (m)")
ap.add_argument("--slab-only", action="store_true", help="omit impactor (for substrate relaxation)")
ap.add_argument("--out", default="orcus_ic.bin")
A = ap.parse_args()

RHO = 2700.0; AB = 2.0           # basalt Tillotson (a+b) -> lithostatic thermal coeff
dx = A.dx; h2 = 2.0*1.3*dx; a = A.D/2.0; m_part = RHO*dx**3

# --- flat basalt slab, surface at z=0, filled downward ---
xs = np.arange(A.xlo+dx/2, A.xhi, dx)
ys = np.arange(-A.yhalf+dx/2, A.yhalf, dx)
zs = np.arange(-dx/2, -A.depth, -dx)
X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij"); n = X.size
u_litho = (A.g/AB)*np.abs(Z)     # lithostatic equilibrium via thermal pressure
fx = (np.abs(X-A.xlo) < h2) | (np.abs(X-A.xhi) < h2)
fy = (np.abs(Y-(-A.yhalf)) < h2) | (np.abs(Y-A.yhalf) < h2)
fb = (Z < -A.depth + 2*dx)
fixed = (fx | fy | fb).astype(np.float64)
slab = np.column_stack([X.ravel(), Y.ravel(), Z.ravel(),
    np.zeros(n), np.zeros(n), np.zeros(n),
    np.full(n, m_part), u_litho.ravel(), np.zeros(n), fixed.ravel()])
print(f"slab: {n:,} ({len(xs)}x{len(ys)}x{len(zs)}), fixed {int(fixed.sum()):,}, "
      f"u_bottom={u_litho.max():.2e} J/kg (melt thresh 1e6)")

# --- basalt impactor sphere, grazing at theta, heading +x ---
if A.slab_only:
    data = slab; print("slab-only (no impactor) -- for substrate relaxation")
else:
    th = np.radians(A.theta); vx = A.U*np.cos(th); vz = -A.U*np.sin(th)
    z0 = a + 5e3                      # bottom starts 5 km above the surface
    na = int(np.ceil(a/dx)) + 1
    ix, iy, iz = np.meshgrid(*(np.arange(-na, na+1),)*3, indexing="ij")
    mask = (ix*dx)**2 + (iy*dx)**2 + (iz*dx)**2 <= a*a
    xp = ix[mask]*dx + A.x0; yp = iy[mask]*dx; zp = iz[mask]*dx + z0; ni = xp.size
    imp = np.column_stack([xp, yp, zp,
        np.full(ni, vx), np.zeros(ni), np.full(ni, vz),
        np.full(ni, m_part), np.zeros(ni), np.full(ni, 2.0), np.zeros(ni)])  # mat=2: basalt impactor (trackable)
    KE = 0.5*ni*m_part*A.U**2
    print(f"impactor: {ni:,}, D={A.D/1e3:.0f}km, v=({vx:.0f},0,{vz:.0f}), "
          f"center=({A.x0/1e3:.0f},0,{z0/1e3:.0f})km, CPPR={a/dx:.1f}, KE={KE:.2e} J ({KE/4.184e18:.2e} Gt)")
    data = np.vstack([slab, imp])
N = data.shape[0]
xlo, xhi = A.xlo-dx, A.xhi+dx
ylo, yhi = -A.yhalf-dx, A.yhalf+dx
zlo, zhi = -A.depth-dx, (80e3 if A.slab_only else z0+a+50e3)
hdr = np.array([xlo, xhi, ylo, yhi, zlo, zhi, float(N)], dtype="<f8")
with open(A.out, "wb") as f:
    hdr.tofile(f); data.astype("<f8").tofile(f)
print(f"wrote {A.out}: N={N:,} ({N*80/1e6:.0f} MB)  domain x[{xlo/1e3:.0f},{xhi/1e3:.0f}] "
      f"y[{ylo/1e3:.0f},{yhi/1e3:.0f}] z[{zlo/1e3:.0f},{zhi/1e3:.0f}] km")
