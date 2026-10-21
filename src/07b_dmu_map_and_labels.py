"""07b_dmu_map_and_labels.py -- Smallest-catchment DMU map and deposit-catchment links.

Part of SpatialDEA-MPM (Nazarpour, submitted to Mathematical Geosciences).
Run from the repository root (see README.md). Set SDEA_BASE to override the base directory.
"""
import os, sys
BASE = os.environ.get("SDEA_BASE", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as _np
if not hasattr(_np,"in1d"): _np.in1d=_np.isin
import numpy as np, pandas as pd, geopandas as gpd, rasterio, json
from rasterio.features import geometry_mask
from pysheds.grid import Grid
grid=Grid.from_raster(BASE+'/work/gee/shz_dem_nasadem.tif'); dem=grid.read_raster(BASE+'/work/gee/shz_dem_nasadem.tif')
fdir=grid.flowdir(grid.resolve_flats(grid.fill_depressions(grid.fill_pits(dem))))
d=pd.read_csv(BASE+'/work/samples_snapped.csv')
dep=gpd.read_file(BASE+'/data/Deposits.shp').to_crs(32639)
ref=rasterio.open(BASE+'/work/gee/shz_dem_nasadem.tif')
bufm=[~geometry_mask([g.buffer(1000)],out_shape=ref.shape,transform=ref.transform) for g in dep.geometry]
hit=np.zeros((len(d),len(dep)),bool)
# also build "compact" DMU map: each pixel -> smallest upstream catchment (<=50 km2) containing it
best=np.full(ref.shape,1e9,np.float32); lab=np.zeros(ref.shape,np.int32)
for k in d.index:
    cm=np.asarray(grid.catchment(x=int(d.loc[k,'col']),y=int(d.loc[k,'row']),fdir=fdir,xytype='index'))>0
    a=d.loc[k,'area_up_km2']
    if a<=50:
        upd=cm&(a<best); best[upd]=a; lab[upd]=k+1
    for j,m in enumerate(bufm): hit[k,j]=(cm&m).any()
np.save(BASE+'/work/dep_hit.npy',hit)
p=ref.profile; p.update(dtype='int32',count=1,nodata=0,compress='deflate')
with rasterio.open(BASE+'/work/dmu_map.tif','w',**p) as o: o.write(lab,1)
area=d.area_up_km2.values
for j,r in dep.iterrows():
    c=np.where(hit[:,j]&(area<=50))[0]
    sel=c[np.argmin(area[c])] if len(c) else -1
    print(f"{r.NAME:20s} {r['type']:9s} n_up={len(c):2d} smallest={sel} ({area[sel] if sel>=0 else 0:.1f} km2)")
print('mapped fraction',(lab>0).mean())
