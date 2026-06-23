#!/usr/bin/env python3
"""Two-panel animation of the Orcus oblique-impact run from gpu_snap_*.bin.
  top:  top-down slab surface relief (crater shape + rim), impactor = red dot
  side: x-z cut near y=0 (slab grey, impactor red) -> grazing skim + skip-off
Run from the dir holding the snapshots:  ../.venv/bin/python ../sph/make_orcus_anim.py
"""
import numpy as np, glob, subprocess, os
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import imageio_ffmpeg
from PIL import Image

def load(fn):
    with open(fn,"rb") as f:
        n=int(np.fromfile(f,dtype=np.int64,count=1)[0])
        return np.fromfile(f,dtype=np.float64,count=10*n).reshape(n,10)

snaps=sorted(glob.glob("gpu_snap_*.bin")); dt=5.0; bs=4000.0
XL,XR,YB,YT,ZB,ZT=-403,453,-203,203,-90,55   # km, fixed frame
os.makedirs("oframes",exist_ok=True); paths=[]
for k,s in enumerate(snaps):
    r=load(s); mat=r[:,9]; x,y,z=r[:,0],r[:,1],r[:,2]
    slab=(mat==0); imp=(mat==2)
    xs,ys,zs=x[slab],y[slab],z[slab]
    ix=np.clip(np.floor((xs-XL*1e3)/bs).astype(int),0,None); iy=np.clip(np.floor((ys-YB*1e3)/bs).astype(int),0,None)
    nx,ny=ix.max()+1,iy.max()+1; g=np.full((nx,ny),-1e9); np.maximum.at(g,(ix,iy),zs)
    surf=np.where(g>-1e8,g,np.nan); rel=(surf-np.nanmedian(surf))/1e3
    icx=x[imp].mean()/1e3 if imp.any() else None; icz=z[imp].mean()/1e3 if imp.any() else None
    fig,(a1,a2)=plt.subplots(2,1,figsize=(12,7),gridspec_kw={"height_ratios":[2,1]})
    a1.imshow(np.clip(rel.T,-6,2),origin="lower",extent=[XL,XR,YB,YT],aspect="equal",cmap="terrain",vmin=-6,vmax=2)
    if icx is not None and XL<icx<XR: a1.plot(icx,0,"o",ms=7,mfc="red",mec="k")  # impactor (rides y~0)
    a1.set_xlim(XL,XR); a1.set_ylim(YB,YT); a1.set_ylabel("y (km)")
    a1.set_title(f"Orcus oblique impact  D=40km  10 km/s  5deg   t={dt*(k+1):.0f} s   (top-down surface relief, km)")
    sel=np.abs(y)<6000
    a2.scatter(x[sel&slab]/1e3,z[sel&slab]/1e3,s=0.4,c="0.55",lw=0)
    if imp.any(): a2.scatter(x[imp]/1e3,z[imp]/1e3,s=3,c="red",lw=0)
    a2.set_xlim(XL,XR); a2.set_ylim(ZB,ZT); a2.set_xlabel("x downrange (km)"); a2.set_ylabel("z (km)")
    a2.axhline(0,color="k",lw=0.4,ls=":")
    fig.tight_layout(); out=f"oframes/f_{k:03d}.png"; fig.savefig(out,dpi=85); plt.close(fig); paths.append(out)
    print("frame",k,s,flush=True)

imgs=[Image.open(p).convert("P",palette=Image.ADAPTIVE) for p in paths]
imgs[0].save("orcus_anim.gif",save_all=True,append_images=imgs[1:],duration=180,loop=0)
ff=imageio_ffmpeg.get_ffmpeg_exe()
subprocess.run([ff,"-y","-framerate","6","-i","oframes/f_%03d.png",
  "-vf","tpad=stop_mode=clone:stop_duration=2,scale=trunc(iw/2)*2:trunc(ih/2)*2",
  "-c:v","libx264","-pix_fmt","yuv420p","-crf","20","-movflags","+faststart","orcus_anim.mp4"],check=True)
print(f"wrote orcus_anim.gif + orcus_anim.mp4 ({len(paths)} frames)")
