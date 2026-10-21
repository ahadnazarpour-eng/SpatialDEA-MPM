"""11_bootstrap_auc.py -- Bootstrap confidence intervals and paired differences.

Part of SpatialDEA-MPM (Nazarpour, submitted to Mathematical Geosciences).
Run from the repository root (see README.md). Set SDEA_BASE to override the base directory.
"""
import os, sys
BASE = os.environ.get("SDEA_BASE", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import sys; pass
import numpy as np, pandas as pd
from dea import auc
R=pd.read_csv(BASE+'/work/results_dmu.csv')
models=[c for c in R.columns if c in ['BLEG Au only','Equal weights','BoD (CCR, unit input)','AR-BoD (alpha=0.05)','Common weights (lambda->inf)'] or c.startswith('GR-DEA (')]
rng=np.random.default_rng(7); B=5000
out=[]
for lab in ['pos','pos_prox']:
    P=np.where(R[lab])[0]; N=np.where(~R[lab])[0]
    boots={m:[] for m in models}
    for b in range(B):
        p=rng.choice(P,len(P)); n=rng.choice(N,len(N))
        for m in models:
            s=R[m].values; boots[m].append(auc(np.r_[s[p],s[n]],np.r_[np.ones(len(p),bool),np.zeros(len(n),bool)]))
    gr=[m for m in models if m.startswith('GR-DEA (')][0]
    for m in models:
        a=np.array(boots[m]); diff=np.array(boots[gr])-a
        out.append(dict(label=lab,model=m,auc=auc(R[m].values,R[lab].values),lo=np.percentile(a,2.5),hi=np.percentile(a,97.5),
                        diff_vs_GR_lo=np.percentile(diff,2.5),diff_vs_GR_hi=np.percentile(diff,97.5),P_GR_better=(diff>0).mean()))
O=pd.DataFrame(out); print(O.round(3).to_string()); O.to_csv(BASE+'/work/bootstrap_auc.csv',index=False)
