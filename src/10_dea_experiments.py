"""10_dea_experiments.py -- Lambda path, L-curve, model comparison, sensitivity, stability.

Part of SpatialDEA-MPM (Nazarpour, submitted to Mathematical Geosciences).
Run from the repository root (see README.md). Set SDEA_BASE to override the base directory.
"""
import os, sys
BASE = os.environ.get("SDEA_BASE", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import sys; pass
import numpy as np, pandas as pd, json
from dea import *
from scipy.stats import spearmanr
d=pd.read_csv(BASE+'/work/dmu_table.csv'); d['au_c']=np.where(d.au<=0.10,0.05,d.au)
D=d[d.area_up_km2<=50].reset_index()
dep=pd.read_json(BASE+'/work/deposit_labels.json')
pos=D['index'].isin([274,80,263,261,265,129]).values
ctrl=D['index'].isin([100,247,97,259,96,105,244,242,93,122,117,118,119]).values
# proximity-based positives: outlet within R of a target deposit
tg=dep[dep.target][['x','y']].values
dmin=np.min(np.hypot(D['xs'].values[:,None]-tg[:,0],D['ys'].values[:,None]-tg[:,1]),1)
pos3=dmin<=2500
print('positives strict',pos.sum(),'proximity(2.5km)',pos3.sum())
W,b=knn_weights(D[['xs','ys']].values,6)
def Y_of(h=1000,t=1.5): return normalise(D[['au_c',f'OH_t{t}',f'FE_t{t}',f'GC_h{h}',f'MC_h{h}']].values)
Y=Y_of()
out={}
# 1. lambda path (alpha=0.05)
lams=[0,0.1,0.3,1,2,3,5,10,20,30,50,100,300,1000]
path=[];TH={}
for lam in lams:
    th,U,E,st=grdea(Y,W,lam,alpha=0.05); TH[lam]=(th,U)
    path.append(dict(lam=lam,sum_theta=th.sum(),energy=E,auc=auc(th,pos),auc_prox=auc(th,pos3),auc_ctrl=auc(th,pos,ctrl),max=th.max()))
P=pd.DataFrame(path)
# L-curve corner: max curvature of (log E, -log sum) wrt log lam over lam>0
Q=P[P.lam>0].reset_index(drop=True)
x=np.log(Q.energy.values); y=np.log(Q.sum_theta.values); t=np.log(Q.lam.values)
dx=np.gradient(x,t); dy=np.gradient(y,t); ddx=np.gradient(dx,t); ddy=np.gradient(dy,t)
kappa=(dx*ddy-dy*ddx)/((dx**2+dy**2)**1.5); Q['kappa']=kappa
lam_star=float(Q.lam[np.argmax(np.abs(kappa[1:-1]))+1]); print(Q.round(4).to_string()); print('lambda* (L-curve)',lam_star)
P.to_csv(BASE+'/work/lambda_path.csv',index=False); Q.to_csv(BASE+'/work/lcurve.csv',index=False)
# 2. model comparison
th_gr,U_gr,E_gr,_=grdea(Y,W,lam_star,alpha=0.05)
models={'BLEG Au only':D.au_c.values,'Equal weights':Y.mean(1)}
models['BoD (CCR, unit input)'],U_bod=bod(Y)
models['AR-BoD (alpha=0.05)'],U_ar=bod(Y,alpha=0.05)
models['SE-AR-BoD'],_=bod(Y,alpha=0.05,super_eff=True)
models['Common weights (lambda->inf)'],u_cw=common_weights(Y)
models[f'GR-DEA (alpha=0.05, lambda*={lam_star:g})']=th_gr
for nm,col in [('fault proximity','FK_h1000'),('fault density','fault_dens'),('local fractal D','Dloc')]:
    Y6=normalise(D[['au_c','OH_t1.5','FE_t1.5','GC_h1000','MC_h1000',col]].values)
    models[f'AR-BoD + {nm}'],_=bod(Y6,alpha=0.05)
    models[f'GR-DEA + {nm}'],_,_,_=grdea(Y6,W,lam_star,alpha=0.05)
rows=[]
for k,v in models.items():
    a,p=auc_perm_p(v,pos,nperm=10000); a3,p3=auc_perm_p(v,pos3,nperm=10000)
    rows.append(dict(model=k,AUC_strict=a,p_strict=p,AUC_prox=a3,p_prox=p3,AUC_vs_ctrl=auc(v,pos,ctrl),n_eff=int((v>0.9999).sum())))
M=pd.DataFrame(rows); print(M.round(3).to_string()); M.to_csv(BASE+'/work/model_comparison.csv',index=False)
# 3. sensitivity to bandwidth h and anomaly threshold t and alpha
sens=[]
for h in [250,500,1000,2000,3000]:
    for t_ in [1.0,1.5,2.0]:
        Yh=Y_of(h,t_)
        for a in [0,0.05,0.1]:
            th,_=bod(Yh,alpha=a); thg,_,_,_=grdea(Yh,W,lam_star,alpha=a)
            sens.append(dict(h=h,t=t_,alpha=a,AUC_ARBoD=auc(th,pos),AUC_GR=auc(thg,pos),AUCprox_GR=auc(thg,pos3),rho_vs_base=spearmanr(thg,th_gr)[0]))
S=pd.DataFrame(sens); S.to_csv(BASE+'/work/sensitivity.csv',index=False); print(S.round(3).to_string())
# 4. stability under multiplicative noise
rng=np.random.default_rng(1); stab=[]
Yraw=D[['au_c','OH_t1.5','FE_t1.5','GC_h1000','MC_h1000']].values
for sd in [0.05,0.1,0.2]:
    for rep in range(20):
        Yp=normalise(Yraw*np.exp(rng.normal(0,sd,Yraw.shape)))
        thp,_,_,_=grdea(Yp,W,lam_star,alpha=0.05); thb,_=bod(Yp)
        stab.append(dict(sd=sd,rep=rep,rho_GR=spearmanr(thp,th_gr)[0],rho_BoD=spearmanr(thb,models['BoD (CCR, unit input)'])[0],auc_GR=auc(thp,pos)))
St=pd.DataFrame(stab); print(St.groupby('sd').agg(['mean','std']).round(3)); St.to_csv(BASE+'/work/stability.csv',index=False)
# save
R=D[['index','id','xs','ys','area_up_km2','au']].copy(); R['pos']=pos; R['pos_prox']=pos3; R['ctrl']=ctrl
for k,v in models.items(): R[k]=v
for j,c in enumerate(['w_Au','w_OH','w_FE','w_GC','w_MC']): R['GR_'+c]=U_gr[:,j]*Y[:,j]  # virtual weights (contribution)
for j,c in enumerate(['w_Au','w_OH','w_FE','w_GC','w_MC']): R['BoD_'+c]=U_bod[:,j]*Y[:,j]
R.to_csv(BASE+'/work/results_dmu.csv',index=False)
json.dump(dict(lam_star=lam_star,knn_bandwidth=b,n=len(D)),open(BASE+'/work/meta.json','w'))
