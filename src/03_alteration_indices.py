"""03_alteration_indices.py -- ASTER and Sentinel-2 alteration indices and vegetation mask.

Part of SpatialDEA-MPM (Nazarpour, submitted to Mathematical Geosciences).
Run from the repository root (see README.md). Set SDEA_BASE to override the base directory.
"""
import os, sys
BASE = os.environ.get("SDEA_BASE", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import rasterio, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
A=rasterio.open(BASE+'/work/aster_refl.tif'); S=rasterio.open(BASE+'/work/s2_refl.tif')
a={b:A.read(i+1).astype('f4') for i,b in enumerate(['B1','B2','B3','B4','B5','B6','B7','B8','B9'])}
s={b:S.read(i+1).astype('f4') for i,b in enumerate(['B2','B3','B4','B8','B11','B12'])}
valid=np.all([v>0 for v in a.values()],0)&np.all([v>0 for v in s.values()],0)
for d in (a,s):
    for k in d: d[k]=np.where(valid,d[k]/10000,np.nan)
ndvi=(s['B8']-s['B4'])/(s['B8']+s['B4'])
mask=valid&(ndvi<0.30)
I={}
I['A_phyllic']=(a['B5']+a['B7'])/a['B6']          # Al-OH (sericite/muscovite)
I['A_argillic']=(a['B4']/a['B5'])*(a['B8']/a['B6'])  # KLI kaolinite (Ninomiya)
I['A_propylitic']=(a['B6']+a['B9'])/(a['B7']+a['B8']) # Mg-OH/carbonate/chlorite
I['A_ferric']=a['B2']/a['B1']                     # Fe3+
I['S_ferric']=s['B4']/s['B2']
I['S_hydroxyl']=s['B11']/s['B12']
I['S_ferrous']=s['B12']/s['B8']+s['B3']/s['B4']
out={}
for k,v in I.items():
    v=np.where(mask,v,np.nan); m=np.nanmean(v); sd=np.nanstd(v)
    z=(v-m)/sd; out[k]=z.astype('f4')
    print(f'{k:14s} mean={m:.3f} sd={sd:.3f} anom(z>1.65)={np.nanmean(z>1.65):.3f}')
print('masked fraction (veg/invalid):',1-mask.mean(), 'veg', np.nanmean(ndvi>=0.3))
np.savez_compressed(BASE+'/work/alteration_z.npz',ndvi=ndvi.astype('f4'),mask=mask,**out)
fig,ax=plt.subplots(2,4,figsize=(20,11))
for axx,(k,v) in zip(ax.flat,list(out.items())+[('NDVI',ndvi)]):
    im=axx.imshow(v,cmap='RdYlBu_r',vmin=-2.5 if k!='NDVI' else 0,vmax=2.5 if k!='NDVI' else .6); axx.set_title(k); axx.axis('off'); plt.colorbar(im,ax=axx,shrink=.6)
fig.tight_layout(); fig.savefig(BASE+'/work/alteration.png',dpi=55)
