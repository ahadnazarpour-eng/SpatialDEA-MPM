"""02_mosaic_gee_tiles.py -- Mosaic ASTER / Sentinel-2 GeoTIFF tiles exported from Earth Engine.

Part of SpatialDEA-MPM (Nazarpour, submitted to Mathematical Geosciences).
Run from the repository root (see README.md). Set SDEA_BASE to override the base directory.
"""
import os, sys
BASE = os.environ.get("SDEA_BASE", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import rasterio, glob, numpy as np
from rasterio.merge import merge
ref=rasterio.open(BASE+'/work/gee/shz_dem_nasadem.tif')
for name,bands in [('aster',['B01','B02','B3N','B04','B05','B06','B07','B08','B09']),('s2',['B2','B3','B4','B8','B11','B12'])]:
    fs=[rasterio.open(f) for f in sorted(glob.glob(f'{BASE}/work/gee/shz_{name}_tile*.tif'))]
    arr,tr=merge(fs,bounds=ref.bounds,res=30,nodata=-32768)
    print(name,arr.shape,tr[:6])
    p=ref.profile.copy(); p.update(count=arr.shape[0],dtype='int16',nodata=-32768,compress='deflate')
    with rasterio.open(f'{BASE}/work/{name}_refl.tif','w',**p) as o:
        o.write(arr); o.descriptions=tuple(bands)
