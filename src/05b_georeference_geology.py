"""05b_georeference_geology.py -- Projective georeferencing of map-derived lithology to the 30 m grid.

Part of SpatialDEA-MPM (Nazarpour, submitted to Mathematical Geosciences).
Run from the repository root (see README.md). Set SDEA_BASE to override the base directory.
"""
import os, sys
BASE = os.environ.get("SDEA_BASE", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np, rasterio
from skimage.transform import ProjectiveTransform
from pyproj import Transformer
from scipy import ndimage
src=np.array([[927,64],[2142,61],[933,1518],[2151,1514]],float)   # (col,row)
dst=np.array([[49.0,34.0],[49.5,34.0],[49.0,33.5],[49.5,33.5]])
pt=ProjectiveTransform(); pt.estimate(dst,src)   # lonlat -> pixel
ref=rasterio.open(BASE+'/work/gee/shz_dem_nasadem.tif'); T=ref.transform; H,W=ref.shape
cols,rows=np.meshgrid(np.arange(W)+.5,np.arange(H)+.5)
X=T.c+cols*T.a; Y=T.f+rows*T.e
lon,lat=Transformer.from_crs(32639,4326,always_xy=True).transform(X,Y)
pix=pt(np.c_[lon.ravel(),lat.ravel()]); pc=np.round(pix[:,0]).astype(int).reshape(H,W); pr=np.round(pix[:,1]).astype(int).reshape(H,W)
M=np.load(BASE+'/work/geomap_masks.npz')
ok=(pr>=0)&(pr<M['gran'].shape[0])&(pc>=0)&(pc<M['gran'].shape[1])
out={}
for k in ['gran','horn']:
    a=np.zeros((H,W),bool); a[ok]=M[k][pr[ok],pc[ok]]; out[k]=a
# distance to granitoid contact (outside and inside)
g=out['gran']; edge=g ^ ndimage.binary_erosion(g)
out['dist_gran_contact']=(ndimage.distance_transform_edt(~edge)*30).astype('f4')
out['dist_gran']=(ndimage.distance_transform_edt(~g)*30).astype('f4')
np.savez_compressed(BASE+'/work/geology_px.npz',**out)
print('granite frac',g.mean().round(3),'hornfels',out['horn'].mean().round(3))
# check: deposits location vs granite distance
import geopandas as gpd
dep=gpd.read_file(BASE+'/data/Deposits.shp').to_crs(32639)
for _,r in dep.iterrows():
    rr,cc=ref.index(r.geometry.x,r.geometry.y); print(f"{r.NAME:20s}{r['type']:9s} dGranContact={out['dist_gran_contact'][rr,cc]:7.0f} inGran={g[rr,cc]}")
