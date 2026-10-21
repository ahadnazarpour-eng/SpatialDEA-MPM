"""01_magnetic_derivatives.py -- RTP, analytic signal and tilt derivative of the aeromagnetic grid.

Part of SpatialDEA-MPM (Nazarpour, submitted to Mathematical Geosciences).
Run from the repository root (see README.md). Set SDEA_BASE to override the base directory.
"""
import os, sys
BASE = os.environ.get("SDEA_BASE", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import geopandas as gpd, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
g=gpd.read_file(BASE+'/data/Mag.shp')
x=np.round(g.geometry.x.values,4); y=np.round(g.geometry.y.values,4)
xs=np.unique(x); ys=np.unique(y)[::-1]
G=np.full((len(ys),len(xs)),np.nan)
G[np.searchsorted(-ys,-y), np.searchsorted(xs,x)]=g.GRID_CODE.values
dx=0.01*111320*np.cos(np.radians(33.75)); dy=0.01*110950
print(G.shape, np.isnan(G).sum(), dx, dy)
# pad with reflection, detrend
p=25; Gp=np.pad(G-G.mean(),p,mode='reflect')
ny,nx=Gp.shape; kx=2*np.pi*np.fft.fftfreq(nx,dx); ky=2*np.pi*np.fft.fftfreq(ny,dy)
KX,KY=np.meshgrid(kx,ky); K=np.hypot(KX,KY)
F=np.fft.fft2(Gp)
# RTP (assume TMI; I=50.3, D=4.6 for Shazand ~ IGRF2015) 
I,D=np.radians(50.3),np.radians(4.6)
l=np.cos(I)*np.sin(D); m=np.cos(I)*np.cos(D); n=np.sin(I)
with np.errstate(divide='ignore',invalid='ignore'):
    th=(n+1j*(l*KX+m*KY)/np.where(K==0,1,K))
    Frtp=F/(th**2); Frtp[K==0]=0
rtp=np.real(np.fft.ifft2(Frtp))
def derivs(Fx):
    dz=np.real(np.fft.ifft2(Fx*K)); ddx=np.real(np.fft.ifft2(Fx*1j*KX)); ddy=np.real(np.fft.ifft2(Fx*1j*KY))
    return dz,ddx,ddy
dz,ddx,ddy=derivs(Frtp)
AS=np.sqrt(dz**2+ddx**2+ddy**2); TDR=np.arctan2(dz,np.hypot(ddx,ddy))
cut=lambda a:a[p:-p,p:-p]
out=dict(lon=xs,lat=ys,tmi=G,rtp=cut(rtp)+G.mean(),AS=cut(AS),TDR=cut(TDR))
np.savez(BASE+'/work/mag_grids.npz',**out)
fig,ax=plt.subplots(1,4,figsize=(18,5))
for a,k in zip(ax,['tmi','rtp','AS','TDR']): im=a.imshow(out[k],extent=[xs[0]-.005,xs[-1]+.005,ys[-1]-.005,ys[0]+.005],cmap='turbo'); a.set_title(k); plt.colorbar(im,ax=a,shrink=.7)
fig.savefig('mag.png',dpi=70)
