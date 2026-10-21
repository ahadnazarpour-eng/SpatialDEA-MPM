"""09_dmu_evidence_table.py -- Catchment-mean outputs for every DMU.

Part of SpatialDEA-MPM (Nazarpour, submitted to Mathematical Geosciences).
Run from the repository root (see README.md). Set SDEA_BASE to override the base directory.
"""
import os, sys
BASE = os.environ.get("SDEA_BASE", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as _np
if not hasattr(_np,"in1d"): _np.in1d=_np.isin
import numpy as np, pandas as pd
from pysheds.grid import Grid
grid=Grid.from_raster(BASE+'/work/gee/shz_dem_nasadem.tif'); dem=grid.read_raster(BASE+'/work/gee/shz_dem_nasadem.tif')
fdir=grid.flowdir(grid.resolve_flats(grid.fill_depressions(grid.fill_pits(dem))))
d=pd.read_csv(BASE+'/work/samples_snapped.csv')
E=np.load(BASE+'/work/evidence_px.npz'); Gm=np.load(BASE+'/work/geology_px.npz'); Fz=np.load(BASE+'/work/faults_px.npz')
valid=E['mask']
H=[250,500,1000,2000,3000]
lay={}
for h in H:
    lay[f'MC_h{h}']=np.exp(-(E['dist_mag_contact']/h)**2)
    lay[f'GC_h{h}']=np.exp(-(Gm['dist_gran_contact']/h)**2)
    lay[f'FK_h{h}']=np.exp(-(Fz['dist_fault']/h)**2)
rows=[]
for k in d.index:
    cm=np.asarray(grid.catchment(x=int(d.loc[k,'col']),y=int(d.loc[k,'row']),fdir=fdir,xytype='index'))>0
    nv=max((cm&valid).sum(),1)
    r={'n_px':int(cm.sum()),'valid_frac':(cm&valid).sum()/cm.sum()}
    for thr in [1.0,1.5,2.0]:
        r[f'OH_t{thr}']=((E['OH']>thr)&valid&cm).sum()/nv; r[f'FE_t{thr}']=((E['FE']>thr)&valid&cm).sum()/nv
    r['AS']=float(E['AS'][cm].mean()); r['gran_frac']=float(Gm['gran'][cm].mean()); r['horn_frac']=float(Gm['horn'][cm].mean()); r['fault_dens']=float(Fz['dens'][cm].mean()); r['Dloc']=float(Fz['Dloc'][cm].mean()); r['inter_dens']=float(Fz['inter'][cm].mean())
    for kk,v in lay.items(): r[kk]=float(v[cm].mean())
    rows.append(r)
out=pd.concat([d,pd.DataFrame(rows)],axis=1)
out.to_csv(BASE+'/work/dmu_table.csv',index=False)
print(out.filter(regex='au|OH_t1.5|FE_t1.5|GC_h1000|MC_h1000|gran_frac').corr(method='spearman').round(2))
