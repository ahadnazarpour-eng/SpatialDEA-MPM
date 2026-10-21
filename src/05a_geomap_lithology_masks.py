"""05a_geomap_lithology_masks.py -- Colour segmentation of granitoids and hornfels on the scanned map.

Part of SpatialDEA-MPM (Nazarpour, submitted to Mathematical Geosciences).
Run from the repository root (see README.md). Set SDEA_BASE to override the base directory.
"""
import os, sys
BASE = os.environ.get("SDEA_BASE", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from PIL import Image; import numpy as np
from scipy import ndimage
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
im=np.asarray(Image.open(BASE+'/data/Shazand_geomap.jpg').convert('RGB')).astype(float)
R,G,B=im[...,0],im[...,1],im[...,2]
def near(c,tol): return np.sqrt(((im-np.array(c))**2).sum(-1))<tol
gran=near((231,132,110),38)|near((211,114,99),32)
horn=near((201,153,134),22)|near((177,131,120),20)|near((165,104,87),22)
# faults: dark magenta/red-violet thin lines
fault=(R-G>45)&(B-G>5)&(R+G+B<480)&~gran
frame=np.zeros(R.shape,bool); frame[66:1512,930:2147]=True
for nm,mk in [('gran',gran),('horn',horn),('fault',fault)]: print(nm,(mk&frame).mean().round(4))
gran_c=ndimage.binary_opening(ndimage.binary_closing(gran&frame,iterations=3),iterations=2)
gran_c=ndimage.binary_fill_holes(gran_c)
lab,n=ndimage.label(gran_c); sz=ndimage.sum(gran_c,lab,range(1,n+1)); keep=np.isin(lab,np.where(sz>150)[0]+1)
horn_c=ndimage.binary_opening(ndimage.binary_closing(horn&frame,iterations=2),iterations=1)
np.savez_compressed(BASE+'/work/geomap_masks.npz',gran=keep,horn=horn_c,fault=fault&frame)
fig,ax=plt.subplots(1,3,figsize=(21,9))
ax[0].imshow(im.astype('uint8')[60:1520,920:2155]); ax[0].set_title('map')
ov=im.copy()/255*0.4; ov[keep]=[1,0,0]; ov[horn_c&~keep]=[0.6,0.3,0]; ax[1].imshow(ov[60:1520,920:2155]); ax[1].set_title('granitoid (red) / hornfels (brown)')
ax[2].imshow((fault&frame)[60:1520,920:2155],cmap='gray_r'); ax[2].set_title('fault-colour pixels')
fig.tight_layout(); fig.savefig(BASE+'/work/geomap_extract.png',dpi=55)
