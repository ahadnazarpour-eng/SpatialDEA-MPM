"""06a_fault_candidate_mask.py -- Candidate fault pixels by local colour contrast.

Part of SpatialDEA-MPM (Nazarpour, submitted to Mathematical Geosciences).
Run from the repository root (see README.md). Set SDEA_BASE to override the base directory.
"""
import os, sys
BASE = os.environ.get("SDEA_BASE", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from PIL import Image; import numpy as np
from scipy import ndimage
im=np.asarray(Image.open(BASE+'/data/Shazand_geomap.jpg').convert('RGB')).astype(float)
bg=np.stack([ndimage.median_filter(im[...,i],size=9) for i in range(3)],-1)
dl=im-bg
dark=dl.sum(-1)
m=(dark<-45)&((dl[...,0]-dl[...,1])>8)&((dl[...,2]-dl[...,1])>-2)
frame=np.zeros(m.shape,bool); frame[68:1510,932:2145]=True; m&=frame
np.save(BASE+'/work/fault_mask2.npy',m)
Image.fromarray(np.where(m[60:1520,925:2155],0,255).astype('uint8')).save(BASE+'/work/fault_mask2.png'); print(m.sum())
