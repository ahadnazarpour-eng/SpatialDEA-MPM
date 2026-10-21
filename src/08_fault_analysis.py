"""08_fault_analysis.py -- Orientation, power law, box-counting fractal, local D, association, Fry.

Part of SpatialDEA-MPM (Nazarpour, submitted to Mathematical Geosciences).
Run from the repository root (see README.md). Set SDEA_BASE to override the base directory.
"""
import os, sys
BASE = os.environ.get("SDEA_BASE", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import json, numpy as np, geopandas as gpd, pandas as pd, rasterio
from shapely.geometry import LineString, MultiLineString, Point
from skimage.transform import ProjectiveTransform
from pyproj import Transformer
from rasterio.features import rasterize
from scipy import ndimage, stats
src=np.array([[927,64],[2142,61],[933,1518],[2151,1514]],float); dst=np.array([[49.0,34.0],[49.5,34.0],[49.0,33.5],[49.5,33.5]])
inv=ProjectiveTransform.from_estimate(src,dst); tr=Transformer.from_crs(4326,32639,always_xy=True)
F=json.load(open(BASE+'/work/faults_final_px.json'))
segs=F['auto']+list(F['manual'].values()); src_tag=['auto']*len(F['auto'])+['manual']*len(F['manual'])
geoms=[]
for s in segs:
    ll=inv(np.array(s)); X,Y=tr.transform(ll[:,0],ll[:,1]); geoms.append(LineString(list(zip(X,Y))))
G=gpd.GeoDataFrame({'source':src_tag},geometry=geoms,crs=32639); G['length_m']=G.length
def az(g):
    c=np.array(g.coords); d=c[-1]-c[0]; return np.degrees(np.arctan2(d[0],d[1]))%180
G['azimuth']=G.geometry.apply(az)
G.to_file(BASE+'/work/faults_final.gpkg',driver='GPKG')
ref=rasterio.open(BASE+'/work/gee/shz_dem_nasadem.tif'); T=ref.transform; H,W=ref.shape
# sheet polygon (UTM) for clipping
ll=np.array([[49,34],[49.5,34],[49.5,33.5],[49,33.5]]); X,Y=tr.transform(ll[:,0],ll[:,1])
from shapely.geometry import Polygon; sheet=Polygon(zip(X,Y)); A_km2=sheet.area/1e6
out={'n_segments':len(G),'total_km':G.length_m.sum()/1000,'sheet_km2':A_km2,'density_km_per_km2':G.length_m.sum()/1000/A_km2}
# ---- 1. orientation statistics (axial, length-weighted)
th=np.radians(G.azimuth.values*2); w=G.length_m.values
C=np.sum(w*np.cos(th))/w.sum(); S=np.sum(w*np.sin(th))/w.sum(); R=np.hypot(C,S)
out['mean_azimuth']=float((np.degrees(np.arctan2(S,C))/2)%180); out['R_axial']=float(R)
# sets: histogram 10-deg bins length-weighted
h,_=np.histogram(G.azimuth,bins=np.arange(0,181,10),weights=w); out['rose_bins_km']=(h/1000).round(2).tolist()
# Rayleigh test on doubled angles (unweighted)
n=len(G); Rn=np.hypot(np.cos(th).sum(),np.sin(th).sum())/n; Z=n*Rn**2; out['rayleigh_Z']=float(Z); out['rayleigh_p']=float(np.exp(-Z))
# ---- 2. length distribution: power-law MLE (Clauset et al. 2009) with KS-optimal Lmin
L=np.sort(G.length_m.values)
best=None
for lmin in np.unique(L)[:-10]:
    x=L[L>=lmin]; a=1+len(x)/np.sum(np.log(x/lmin))
    cdf_emp=np.arange(1,len(x)+1)/len(x); cdf_fit=1-(x/lmin)**(1-a); D=np.max(np.abs(cdf_emp-cdf_fit))
    if best is None or D<best[2]: best=(lmin,a,D,len(x))
lmin,a,D,nt=best; out['powerlaw']=dict(Lmin_m=float(lmin),alpha=float(a),alpha_se=float((a-1)/np.sqrt(nt)),KS=float(D),n_tail=int(nt))
# compare with exponential (likelihood ratio)
x=L[L>=lmin]; lam=1/np.mean(x-lmin)
ll_pl=np.sum(np.log((a-1)/lmin)-a*np.log(x/lmin)); ll_ex=np.sum(np.log(lam)-lam*(x-lmin))
out['LR_pl_vs_exp']=float(ll_pl-ll_ex)
# ---- 3. box-counting fractal dimension (global)
res=10.0  # 10 m raster for box counting
xmin,ymin,xmax,ymax=sheet.bounds; Wb=int((xmax-xmin)/res)+1; Hb=int((ymax-ymin)/res)+1
from affine import Affine
Tb=Affine(res,0,xmin,0,-res,ymax)
Fb=rasterize([(g,1) for g in G.geometry],out_shape=(Hb,Wb),transform=Tb,all_touched=True).astype(bool)
Sm=rasterize([(sheet,1)],out_shape=(Hb,Wb),transform=Tb).astype(bool)
def boxcount(Fm,Mm,sizes_px):
    N=[]
    for s in sizes_px:
        hh=(Fm.shape[0]//s)*s; ww=(Fm.shape[1]//s)*s
        f=Fm[:hh,:ww].reshape(hh//s,s,ww//s,s).any(axis=(1,3)); m=Mm[:hh,:ww].reshape(hh//s,s,ww//s,s).mean(axis=(1,3))>0.5
        N.append((f&m).sum())
    return np.array(N)
sizes_m=np.array([125,250,500,1000,2000,4000,8000]); sp=(sizes_m/res).astype(int)
N=boxcount(Fb,Sm,sp)
slope,inter,r,p,se=stats.linregress(np.log(1/sizes_m),np.log(N))
out['boxcount']=dict(sizes_m=sizes_m.tolist(),N=N.tolist(),D=float(slope),D_se=float(se),r2=float(r**2))
# null model: random lines with same length & azimuth distribution, uniformly placed (100 realisations)
rng=np.random.default_rng(0); Dnull=[]
for k in range(100):
    gs=[]
    for g in G.geometry:
        Lm=g.length; az_=np.radians(rng.choice(G.azimuth.values)); 
        while True:
            cx=rng.uniform(xmin,xmax); cy=rng.uniform(ymin,ymax)
            if sheet.contains(Point(cx,cy)): break
        dx,dy=np.sin(az_)*Lm/2,np.cos(az_)*Lm/2; gs.append(LineString([(cx-dx,cy-dy),(cx+dx,cy+dy)]))
    Fr=rasterize([(g,1) for g in gs],out_shape=(Hb,Wb),transform=Tb,all_touched=True).astype(bool)
    Nr=boxcount(Fr,Sm,sp); Dnull.append(stats.linregress(np.log(1/sizes_m),np.log(Nr))[0])
Dnull=np.array(Dnull); out['boxcount_null']=dict(mean=float(Dnull.mean()),sd=float(Dnull.std()),p_lower=float((Dnull<=slope).mean()))
np.save(BASE+'/work/Dnull.npy',Dnull)
# ---- 4. local fractal dimension (moving window 6 km, boxes 125-1000 m), 1 km step
win=6000; step=1000; boxes=np.array([125,250,500,1000]); bp=(boxes/res).astype(int)
xs=np.arange(xmin+win/2,xmax-win/2+1,step); ys=np.arange(ymax-win/2,ymin+win/2-1,-step)
Dloc=np.full((len(ys),len(xs)),np.nan); Dens=np.full_like(Dloc,np.nan)
wp=int(win/res)
for i,yc in enumerate(ys):
    for j,xc in enumerate(xs):
        c0=int((xc-win/2-xmin)/res); r0=int((ymax-(yc+win/2))/res)
        f=Fb[r0:r0+wp,c0:c0+wp]; m=Sm[r0:r0+wp,c0:c0+wp]
        if m.mean()<0.9: continue
        Nn=boxcount(f,m,bp)
        Dens[i,j]=f.sum()*res/1000/((win/1000)**2)
        if (Nn>0).all() and Nn[0]>=8: Dloc[i,j]=stats.linregress(np.log(1/boxes),np.log(Nn))[0]
        else: Dloc[i,j]=np.nan
np.savez(BASE+'/work/fractal_local.npz',D=Dloc,dens=Dens,xs=xs,ys=ys,win=win)
# resample local D and density to 30 m grid (bilinear)
from scipy.interpolate import RegularGridInterpolator
cols,rows=np.meshgrid(np.arange(W)+.5,np.arange(H)+.5); XX=T.c+cols*T.a; YY=T.f+rows*T.e
def rs(Z):
    Zf=np.where(np.isnan(Z),np.nanmedian(Z),Z)
    return RegularGridInterpolator((ys[::-1],xs),Zf[::-1],bounds_error=False,fill_value=None)((YY,XX)).astype('f4')
Dloc30=rs(Dloc); Dens30=rs(Dens)
# ---- 5. fault distance, kernel density, intersections
F30=rasterize([(g,1) for g in G.geometry],out_shape=(H,W),transform=T,all_touched=True).astype(bool)
dist=ndimage.distance_transform_edt(~F30)*30
inter=[]
sidx=G.sindex
for i,g in enumerate(G.geometry):
    for j in sidx.query(g,predicate='intersects'):
        if j>i:
            ip=g.intersection(G.geometry.iloc[j])
            if ip.geom_type=='Point': inter.append(ip)
out['n_intersections']=len(inter)
Ip=np.zeros((H,W)); 
for p in inter:
    r_,c_=ref.index(p.x,p.y)
    if 0<=r_<H and 0<=c_<W: Ip[r_,c_]+=1
Iden=ndimage.gaussian_filter(Ip,sigma=1000/30)*(1e6/900)*2*np.pi*(1000/30)**2/ (2*np.pi*(1000/30)**2)  # per km2 (Gaussian-weighted)
Kden=ndimage.gaussian_filter(F30.astype(float),sigma=1000/30)*(30/1000)*(1e6/900)/1000*1000/30*30/1000  # placeholder, recomputed below
Kden=ndimage.gaussian_filter(F30.astype(float),sigma=1000/30)/0.03  # km of fault per km2 (cell fraction / cell size in km)
np.savez_compressed(BASE+'/work/faults_px.npz',dist_fault=dist.astype('f4'),fault=F30,dens=Kden.astype('f4'),inter=Iden.astype('f4'),Dloc=Dloc30,Dens6=Dens30)
# ---- 6. spatial association with deposits: distance-distribution analysis (Berman 1977) + Monte Carlo
dep=gpd.read_file(BASE+'/data/Deposits.shp').to_crs(32639)
dep['grp']=np.where(dep['type'].isin(['Cu','W','Au_Vein']),'target',np.where(dep['type'].isin(['Talc','Si']),'control','other'))
inside=rasterize([(sheet,1)],out_shape=(H,W),transform=T).astype(bool)
def at(arr,pts): return np.array([arr[ref.index(p.x,p.y)] for p in pts])
assoc={}
for name,arr in [('dist_fault_m',dist),('fault_density',Kden),('intersection_density',Iden),('local_D',Dloc30)]:
    bg=arr[inside]; r={}
    for grp in ['target','control']:
        v=at(arr,dep[dep.grp==grp].geometry)
        # Monte Carlo: mean value at random points (same n) -> p-value
        sims=np.array([np.mean(rng.choice(bg,len(v))) for _ in range(10000)])
        pl=(sims<=v.mean()).mean(); r[grp]=dict(mean=float(v.mean()),median=float(np.median(v)),bg_mean=float(bg.mean()),p_lower=float(pl),p_upper=float((sims>=v.mean()).mean()))
    assoc[name]=r
out['association']=assoc
# distance distribution curves (deposit cumulative vs area cumulative)
dd=np.linspace(0,5000,101); bg=dist[inside]
cur={'d_m':dd.tolist(),'area':[(bg<=x).mean() for x in dd]}
for grp in ['target','control']:
    v=at(dist,dep[dep.grp==grp].geometry); cur[grp]=[(v<=x).mean() for x in dd]
    cur[grp+'_D']=float(np.max(np.array(cur[grp])-np.array(cur['area']))); cur[grp+'_dmax']=float(dd[np.argmax(np.array(cur[grp])-np.array(cur['area']))])
json.dump(cur,open(BASE+'/work/dist_curves.json','w'))
# ---- 7. Fry analysis (target deposits)
P=np.c_[[g.x for g in dep[dep.grp=='target'].geometry],[g.y for g in dep[dep.grp=='target'].geometry]]
fry=np.array([P[i]-P[j] for i in range(len(P)) for j in range(len(P)) if i!=j])
np.save(BASE+'/work/fry.npy',fry)
fa=np.degrees(np.arctan2(fry[:,0],fry[:,1]))%180; fd=np.hypot(fry[:,0],fry[:,1])
th2=np.radians(2*fa); Rf=np.hypot(np.cos(th2).mean(),np.sin(th2).mean()); out['fry']=dict(n_vectors=len(fry),mean_trend=float((np.degrees(np.arctan2(np.sin(th2).mean(),np.cos(th2).mean()))/2)%180),R=float(Rf),
   mean_trend_lt10km=float((np.degrees(np.arctan2(np.sin(th2[fd<10000]).mean(),np.cos(th2[fd<10000]).mean()))/2)%180), n_lt10km=int((fd<10000).sum()))
json.dump(out,open(BASE+'/work/fault_stats.json','w'),indent=1)
print(json.dumps({k:v for k,v in out.items() if k!='rose_bins_km'},indent=1))
