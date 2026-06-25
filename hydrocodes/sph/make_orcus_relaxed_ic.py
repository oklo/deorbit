#!/usr/bin/env python3
"""Build the Orcus impact IC from a RELAXED slab snapshot (gpu_snap from a damped
gravity-only run): keep the settled positions + internal energy, restore mass +
fixed flags, zero the residual velocity, and drop the grazing impactor above the
settled surface. Same impactor as v1 (D=40km, 10 km/s, 5deg) -> only the substrate
changed.  -> orcus_relaxed_ic.bin
"""
import numpy as np, argparse
ap = argparse.ArgumentParser()
ap.add_argument("--snap", default="gpu_snap_004.bin")
ap.add_argument("--dx", type=float, default=3000.0)
ap.add_argument("--D", type=float, default=40e3); ap.add_argument("--U", type=float, default=10000.0)
ap.add_argument("--theta", type=float, default=5.0)
ap.add_argument("--xlo", type=float, default=-400e3); ap.add_argument("--xhi", type=float, default=450e3)
ap.add_argument("--yhalf", type=float, default=200e3); ap.add_argument("--depth", type=float, default=120e3)
ap.add_argument("--x0", type=float, default=-250e3); ap.add_argument("--out", default="orcus_relaxed_ic.bin")
A = ap.parse_args()
RHO = 2700.0; dx = A.dx; a = A.D/2; m_part = RHO*dx**3; wall = 6.0*dx   # must match make_orcus_ic.py

def load(fn):
    with open(fn,"rb") as f:
        n = int(np.fromfile(f,dtype=np.int64,count=1)[0])
        return np.fromfile(f,dtype=np.float64,count=10*n).reshape(n,10)  # x y z vx vy vz u rho h mat

s = load(A.snap); x,y,z,u = s[:,0],s[:,1],s[:,2],s[:,6]; n = len(x)
fx = (np.abs(x-A.xlo)<wall)|(np.abs(x-A.xhi)<wall)
fy = (np.abs(y-(-A.yhalf))<wall)|(np.abs(y-A.yhalf)<wall)
fb = (z < -A.depth + wall)
fixed = (fx|fy|fb).astype(np.float64)
slab = np.column_stack([x,y,z, s[:,3],s[:,4],s[:,5],   # keep relaxed velocities (~0.3 m/s; zeroing re-perturbs)
    np.full(n,m_part), u, np.zeros(n), fixed])                         # mat=0
surf = np.median(np.sort(z)[-15000:])
print(f"relaxed slab: {n:,} particles, fixed {int(fixed.sum()):,}, settled surface z={surf:.0f} m")

th = np.radians(A.theta); vx = A.U*np.cos(th); vz = -A.U*np.sin(th); z0 = surf + a + 5e3
na = int(np.ceil(a/dx))+1
ix,iy,iz = np.meshgrid(*(np.arange(-na,na+1),)*3, indexing="ij")
mask = (ix*dx)**2+(iy*dx)**2+(iz*dx)**2 <= a*a
xp = ix[mask]*dx+A.x0; yp = iy[mask]*dx; zp = iz[mask]*dx+z0; ni = xp.size
imp = np.column_stack([xp,yp,zp, np.full(ni,vx),np.zeros(ni),np.full(ni,vz),
    np.full(ni,m_part), np.zeros(ni), np.full(ni,2.0), np.zeros(ni)])
print(f"impactor: {ni:,}, D={A.D/1e3:.0f}km, v=({vx:.0f},0,{vz:.0f}), z0={z0:.0f} m ({(z0-surf)/1e3:.1f} km above settled surface)")

data = np.vstack([slab,imp]); N = len(data)
hdr = np.array([A.xlo-dx,A.xhi+dx,-A.yhalf-dx,A.yhalf+dx,-A.depth-dx,z0+a+50e3,float(N)],dtype="<f8")
with open(A.out,"wb") as f: hdr.tofile(f); data.astype("<f8").tofile(f)
print(f"wrote {A.out}: N={N:,}")
