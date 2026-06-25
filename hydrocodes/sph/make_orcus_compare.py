#!/usr/bin/env python3
"""Side-by-side MOLA-style shaded relief, CONSISTENT km scale:
  left  = REAL Orcus Patera (MOLA), rotated so its long axis runs left->right
  right = OUR SPH crater (same orientation, same km box, same scale)
Shared horizontal colorbar below both panels.
"""
import numpy as np, glob
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
from scipy import ndimage

BOX_W, BOX_H = 480.0, 260.0          # km box (same for both panels) -- the "scale"
ROT = -1                              # sign of rotation (flip if it comes out vertical)

# ---------------- REAL Orcus (MOLA) ----------------
dem = np.flipud(np.load("../mars_data/real_orcus_dem.npy"))   # rows now lat 8->21 (up), cols lon 173->185
pxlat = 463.0; pxlon = 463.0*np.cos(np.radians(14.5))
rel0 = dem - np.median(dem)
yy,xx = np.where(rel0 < -400)                                  # depressed pixels = basin
X=(xx-xx.mean())*pxlon; Y=(yy-yy.mean())*pxlat
w,v = np.linalg.eigh(np.cov(np.vstack([X,Y])))
ax = v[:,np.argmax(w)]; ang = np.degrees(np.arctan2(ax[1],ax[0]))
print("basin long-axis angle = %.1f deg from horizontal"%ang)
demr = ndimage.rotate(dem, ROT*ang, reshape=True, cval=np.nan, order=1)
relr = demr - np.nanmedian(demr)
# crop the box centered on the depressed centroid of the rotated map
m = relr < -400; yc,xc = np.argwhere(m).mean(0)
hw,hh = int(BOX_W*1e3/2/pxlon), int(BOX_H*1e3/2/pxlat)
R = relr[int(yc)-hh:int(yc)+hh, int(xc)-hw:int(xc)+hw]
R = np.nan_to_num(R, nan=0.0)
print("REAL box: %dx%d px"%(R.shape[1],R.shape[0]))

# ---------------- OUR crater ----------------
def load(fn):
    with open(fn,"rb") as f:
        n=int(np.fromfile(f,dtype=np.int64,count=1)[0]); return np.fromfile(f,dtype=np.float64,count=10*n).reshape(n,10)
s=load(sorted(glob.glob("../gpu/gpu_snap_*.bin"))[-1]); mat=s[:,9]; x,y,z=s[:,0],s[:,1],s[:,2]
vmag=np.linalg.norm(s[:,3:6],axis=1)
sb=(mat==0)&(vmag<500.0); xs,ys,zs=x[sb],y[sb],z[sb]; bs=4000.0   # |v|<500: settled ground (mask airborne ejecta)
ix=np.floor((xs-xs.min())/bs).astype(int); iy=np.floor((ys-ys.min())/bs).astype(int)
nx,ny=ix.max()+1,iy.max()+1; g=np.full((nx,ny),-1e9); np.maximum.at(g,(ix,iy),zs)
surf=np.where(g>-1e8,g,np.nan); surf=surf-np.nanmedian(surf)
xc_km=xs.min()/1e3+np.arange(nx)*bs/1e3; yc_km=ys.min()/1e3+np.arange(ny)*bs/1e3
# center box on the deep crater centroid
dm=surf< -1500; cxk=xc_km[np.argwhere(dm)[:,0]].mean(); cyk=yc_km[np.argwhere(dm)[:,1]].mean()
def near(arr,c,half): return (arr>c-half)&(arr<c+half)
xi=near(xc_km,cxk,BOX_W/2); yi=near(yc_km,cyk,BOX_H/2)
O=np.nan_to_num(surf[np.ix_(xi,yi)].T, nan=0.0)
print("OURS box: %.0f x %.0f km"%(xi.sum()*bs/1e3,yi.sum()*bs/1e3))

# ---------------- measured basin dimensions (depression extent) ----------------
def dims(arr,pxx,pyy,thr):
    m=arr<thr
    if not m.any(): return 0.,0.,0.
    rr,cc=np.where(m); L=(cc.max()-cc.min()+1)*pxx/1e3; W=(rr.max()-rr.min()+1)*pyy/1e3
    return L,W,L/max(W,1e-9)
Ld,Wd,ed=dims(R,pxlon,pxlat,-400.); Lo,Wo,eo=dims(O,bs,bs,-600.)
print("REAL basin %.0fx%.0f km e=%.2f | OURS basin %.0fx%.0f km e=%.2f"%(Ld,Wd,ed,Lo,Wo,eo))

# ---------------- plot: same scale, colorbar below ----------------
ls=LightSource(azdeg=315,altdeg=45); cmap=plt.cm.terrain; vmin,vmax=-2.5,2.0
fig=plt.figure(figsize=(14,6.5))
a1=fig.add_axes([0.04,0.18,0.45,0.74]); a2=fig.add_axes([0.51,0.18,0.45,0.74])
rgb1=ls.shade(R,cmap=cmap,vert_exag=3,dx=pxlon,dy=pxlat,blend_mode="soft",vmin=vmin*1e3,vmax=vmax*1e3)
a1.imshow(rgb1,extent=[0,BOX_W,0,BOX_H],origin="lower",aspect="equal")
a1.set_title("REAL Orcus Patera (MOLA)"); a1.set_xlabel("km"); a1.set_ylabel("km")
a1.text(0.03,0.96,"%.0f x %.0f km   e=%.1f"%(Ld,Wd,ed),transform=a1.transAxes,va="top",
        fontsize=11,bbox=dict(fc="white",alpha=0.75,ec="none"))
rgb2=ls.shade(O,cmap=cmap,vert_exag=28,dx=bs,dy=bs,blend_mode="soft",vmin=vmin*1e3,vmax=vmax*1e3)
a2.imshow(rgb2,extent=[0,BOX_W,0,BOX_H],origin="lower",aspect="equal")
a2.set_title("OUR SPH impact (40 km, 10 km/s, 5deg; transient)"); a2.set_xlabel("km")
a2.text(0.03,0.96,"%.0f x %.0f km   e=%.1f"%(Lo,Wo,eo),transform=a2.transAxes,va="top",
        fontsize=11,bbox=dict(fc="white",alpha=0.75,ec="none"))
cax=fig.add_axes([0.30,0.07,0.40,0.035])
sm=plt.cm.ScalarMappable(cmap=cmap,norm=plt.Normalize(vmin,vmax))
fig.colorbar(sm,cax=cax,orientation="horizontal",label="elevation rel. to plain (km)")
fig.savefig("orcus_compare.png",dpi=150); print("wrote orcus_compare.png")
