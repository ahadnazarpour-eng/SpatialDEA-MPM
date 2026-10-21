"""07a_catchments.py -- D8 routing, sample snapping and nested catchments.

Part of SpatialDEA-MPM (Nazarpour, submitted to Mathematical Geosciences).
Run from the repository root (see README.md). Set SDEA_BASE to override the base directory.
"""
import os, sys
BASE = os.environ.get("SDEA_BASE", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as _np
if not hasattr(_np,"in1d"): _np.in1d=_np.isin
import numpy as np, pandas as pd, rasterio
from pysheds.grid import Grid
grid = Grid.from_raster(BASE+'/work/gee/shz_dem_nasadem.tif')
dem = grid.read_raster(BASE+'/work/gee/shz_dem_nasadem.tif')
pit = grid.fill_pits(dem); fl = grid.fill_depressions(pit); inf = grid.resolve_flats(fl)
fdir = grid.flowdir(inf); acc = grid.accumulation(fdir)
d=pd.read_excel(BASE+'/data/Shazand_BLEG.xlsx'); d.columns=['no','id','x','y','au']
A=np.asarray(acc); T=grid.affine
snap=[]; 
for x,y in zip(d.x,d.y):
    c,r = ~T*(x,y); c,r=int(c),int(r); w=7  # 210 m window
    r0,c0=max(r-w,0),max(c-w,0); sub=A[r0:r+w+1, c0:c+w+1]
    RR,CC=np.mgrid[0:sub.shape[0],0:sub.shape[1]]; dist=np.hypot(RR+r0-r,CC+c0-c)
    ok=sub>=278
    if ok.any(): i=np.unravel_index(np.argmin(np.where(ok,dist,1e9)),sub.shape)
    else: i=np.unravel_index(np.argmax(sub),sub.shape)
    rr,cc=r0+i[0], c0+i[1]; snap.append((rr,cc,A[rr,cc]))
d['row'],d['col'],d['acc']=zip(*snap)
xs,ys = T*(d.col.values+0.5, d.row.values+0.5); d['xs']=xs; d['ys']=ys
d['snapdist']=np.hypot(d['xs']-d.x,d['ys']-d.y)
lab=np.zeros(A.shape,np.int32)
for k in d.sort_values('acc').index:
    cm=np.asarray(grid.catchment(x=int(d.loc[k,'col']),y=int(d.loc[k,'row']),fdir=fdir,xytype='index'))
    lab[(cm>0)&(lab==0)]=k+1
d['area_up_km2']=d.acc*900/1e6
d['area_excl_km2']=[ (lab==k+1).sum()*900/1e6 for k in d.index]
print(d[['acc','snapdist','area_up_km2','area_excl_km2']].describe())
prof=rasterio.open(BASE+'/work/gee/shz_dem_nasadem.tif').profile; prof.update(dtype='int32',count=1,nodata=0,compress='deflate')
with rasterio.open(BASE+'/work/catchments.tif','w',**prof) as o: o.write(lab,1)
np.save(BASE+'/work/acc.npy',A.astype(np.float32))
d.to_csv(BASE+'/work/samples_snapped.csv',index=False)
print('labelled fraction', (lab>0).mean())
