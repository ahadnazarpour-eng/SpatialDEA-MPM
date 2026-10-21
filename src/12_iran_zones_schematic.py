"""12_iran_zones_schematic.py -- Schematic structural zones for the location map.

Part of SpatialDEA-MPM (Nazarpour, submitted to Mathematical Geosciences).
Run from the repository root (see README.md). Set SDEA_BASE to override the base directory.
"""
import os, sys
BASE = os.environ.get("SDEA_BASE", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np, geopandas as gpd, json
from shapely.geometry import LineString
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
# GCPs: pixel in Deevsalar et al. (2018) Fig.1b (793x613 zoom) -> lon,lat
gcp_px=np.array([[45,60],[148,393],[473,548],[360,160]],float)
gcp_ll=np.array([[44.05,39.4],[48.0,30.0],[56.3,26.6],[54.0,36.9]])
A=np.c_[gcp_px,np.ones(4)]; M,res,*_=np.linalg.lstsq(A,gcp_ll,rcond=None)
f=lambda P: np.c_[np.asarray(P,float),np.ones(len(P))]@M
print('GCP residuals (deg):',np.round(f(gcp_px)-gcp_ll,2))
ZMT=[(72,205),(110,228),(150,262),(185,285),(215,325),(250,345),(290,370),(330,395),(355,420),(370,445),(400,470),(440,476),(470,500),(480,548)]
SaSZ_NE=[(60,128),(130,180),(200,240),(270,300),(330,360),(395,432)]
UDMA_NE=[(57,100),(130,150),(210,210),(290,280),(370,350),(430,420),(470,465)]
lines={'Main Zagros Thrust':ZMT,'SaSZ–UDMA boundary':SaSZ_NE,'UDMA NE boundary':UDMA_NE}
G=gpd.GeoDataFrame({'name':list(lines)},geometry=[LineString(f(v)) for v in lines.values()],crs=4326)
G.to_file(BASE+'/work/zones_schematic.geojson',driver='GeoJSON')
C=gpd.read_file(BASE+'/data/countries_LSIB.geojson')
fig,ax=plt.subplots(figsize=(7,6)); C.boundary.plot(ax=ax,color='k',lw=.6)
G.plot(ax=ax,column='name',legend=True,lw=2); ax.plot(49.25,33.75,'r*',ms=14); ax.set_xlim(43,64); ax.set_ylim(24,40.5)
fig.savefig(BASE+'/work/zones_check.png',dpi=80)
# calibrate: MZT must pass the Main Zagros Reverse Fault mapped in the SW corner of the Shazand sheet (~49.05E, 33.55N)
from shapely.affinity import translate
z=G[G.name=='Main Zagros Thrust'].geometry.iloc[0]
from shapely.geometry import LineString as LS
xs=np.array(z.coords)[:,0]; ys=np.array(z.coords)[:,1]; o=np.argsort(xs); lat_at=np.interp(49.05,xs[o],ys[o])
dy=33.55-lat_at; print('shift lat',dy)
G['geometry']=[translate(g,0,dy) for g in G.geometry]
G.to_file(BASE+'/work/zones_schematic.geojson',driver='GeoJSON')
fig,ax=plt.subplots(figsize=(7,6))
for g in C.geometry: 
    for p in getattr(g,'geoms',[g]):
        if p.geom_type=='LineString': x,y=p.xy; ax.plot(x,y,'k',lw=.5)
G.plot(ax=ax,column='name',legend=True,lw=2); ax.plot(49.25,33.75,'r*',ms=14); ax.set_xlim(43,64); ax.set_ylim(24,40.5)
fig.savefig(BASE+'/work/zones_check.png',dpi=80)
