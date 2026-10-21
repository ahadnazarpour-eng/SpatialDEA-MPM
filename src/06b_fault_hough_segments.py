"""06b_fault_hough_segments.py -- Skeletonisation and probabilistic Hough segments.

Part of SpatialDEA-MPM (Nazarpour, submitted to Mathematical Geosciences).
Run from the repository root (see README.md). Set SDEA_BASE to override the base directory.
"""
import os, sys
BASE = os.environ.get("SDEA_BASE", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from PIL import Image, ImageDraw; import numpy as np
from skimage.transform import probabilistic_hough_line
from skimage.morphology import skeletonize, binary_closing, disk
m=np.load(BASE+'/work/fault_mask2.npy')
m=binary_closing(m,disk(1))
sk=skeletonize(m)
lines=probabilistic_hough_line(sk,threshold=8,line_length=22,line_gap=6,rng=0)
print('segments',len(lines))
L=np.array([[x0,y0,x1,y1] for (x0,y0),(x1,y1) in lines],float)
ln=np.hypot(L[:,2]-L[:,0],L[:,3]-L[:,1]); print('total length px',ln.sum().round(),'km',(ln.sum()*38/1000).round())
np.save(BASE+'/work/fault_segments_px.npy',L)
im=Image.open(BASE+'/data/Shazand_geomap.jpg').convert('RGB').crop((925,60,2155,1520))
dr=ImageDraw.Draw(im)
for x0,y0,x1,y1 in L: dr.line([(x0-925,y0-60),(x1-925,y1-60)],fill=(0,0,0),width=3)
im.save(BASE+'/work/fault_hough.png')
