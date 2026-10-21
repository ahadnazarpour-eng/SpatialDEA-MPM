"""04_evidence_layers.py -- Hydroxyl / ferric composites and magnetic contacts on the 30 m grid.

Part of SpatialDEA-MPM (Nazarpour, submitted to Mathematical Geosciences).
Run from the repository root (see README.md). Set SDEA_BASE to override the base directory.
"""
import os, sys
BASE = os.environ.get("SDEA_BASE", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np, rasterio, geopandas as gpd
from scipy.interpolate import RegularGridInterpolator
from scipy import ndimage
from pyproj import Transformer
ref=rasterio.open(BASE+'/work/gee/shz_dem_nasadem.tif'); T=ref.transform; H,W=ref.shape
alt=np.load(BASE+'/work/alteration_z.npz')
OH=np.nanmean(np.stack([alt['A_phyllic'],alt['A_argillic'],alt['S_hydroxyl']]),0)
FE=np.nanmean(np.stack([alt['A_ferric'],alt['S_ferric']]),0)
def zs(v): return (v-np.nanmean(v))/np.nanstd(v)
OH,FE=zs(OH),zs(FE)
# pixel coordinates -> lon/lat
cols,rows=np.meshgrid(np.arange(W)+.5,np.arange(H)+.5)
X=T.c+cols*T.a; Y=T.f+rows*T.e
lon,lat=Transformer.from_crs(32639,4326,always_xy=True).transform(X,Y)
mg=np.load(BASE+'/work/mag_grids.npz')
interp=lambda G: RegularGridInterpolator((mg['lat'][::-1],mg['lon']),G[::-1],bounds_error=False,fill_value=None)((lat,lon))
AS=interp(mg['AS']); TDR=interp(mg['TDR'])
# magnetic contacts: TDR zero crossings
zc=(np.sign(TDR)!=np.sign(np.roll(TDR,1,0)))|(np.sign(TDR)!=np.sign(np.roll(TDR,1,1)))
zc[0,:]=zc[:,0]=False
dist_contact=ndimage.distance_transform_edt(~zc)*30.0
np.savez_compressed(BASE+'/work/evidence_px.npz',OH=OH.astype('f4'),FE=FE.astype('f4'),AS=AS.astype('f4'),TDR=TDR.astype('f4'),dist_mag_contact=dist_contact.astype('f4'),mask=alt['mask'])
print('OH anomaly frac',np.nanmean(OH>1.5),'FE',np.nanmean(FE>1.5),'dist contact median',np.median(dist_contact))
