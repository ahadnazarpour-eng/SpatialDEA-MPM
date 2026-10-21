"""07c_deposit_table.py -- Deposit table with catchment links (used for validation labels).

Part of SpatialDEA-MPM (Nazarpour, submitted to Mathematical Geosciences).
"""
import os, sys
BASE = os.environ.get("SDEA_BASE", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np, pandas as pd, geopandas as gpd, rasterio
from rasterio.features import geometry_mask
dep=gpd.read_file(BASE+'/data/Deposits.shp').to_crs(32639)
dep['target']=dep['type'].isin(['Cu','Au_Vein','W'])
lab=rasterio.open(BASE+'/work/catchments.tif'); L=lab.read(1)
rows=[]
for i,r in dep.iterrows():
    buf=r.geometry.buffer(1000)
    m=~geometry_mask([buf],out_shape=L.shape,transform=lab.transform)
    ids=np.unique(L[m]); ids=ids[ids>0]-1
    v=L[lab.index(r.geometry.x,r.geometry.y)]-1
    rows.append({'name':r.NAME,'type':r['type'],'target':r.target,'x':r.geometry.x,'y':r.geometry.y,'basin_at_point':int(v),'basins_1km':list(map(int,ids))})
D=pd.DataFrame(rows); print(D.to_string())
# near-duplicates among targets
t=D[D.target]
from scipy.spatial.distance import pdist,squareform
dd=squareform(pdist(t[['x','y']].values)); np.fill_diagonal(dd,1e9); print('min pair dist targets', dd.min(1).round(0))
D.to_json(BASE+'/work/deposit_labels.json',orient='records')
