"""06c_fault_qc_tiles.py -- Gridded, numbered tiles for visual quality control.

Part of SpatialDEA-MPM (Nazarpour, submitted to Mathematical Geosciences).
Run from the repository root (see README.md). Set SDEA_BASE to override the base directory.
"""
import os, sys
BASE = os.environ.get("SDEA_BASE", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from PIL import Image, ImageDraw, ImageFont
import numpy as np, sys
L=np.load(BASE+'/work/fault_segments_px.npy')   # mask2-based 365 segments
os.makedirs(BASE+'/work/qc',exist_ok=True)
im=Image.open(BASE+'/data/Shazand_geomap.jpg').convert('RGB')
x0,y0,x1,y1=925,60,2155,1520; tw,th=410,365; S=3
try: font=ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf',16)
except: font=None
for ty in range(y0,y1,th):
    for tx in range(x0,x1,tw):
        c=im.crop((tx,ty,min(tx+tw,x1),min(ty+th,y1))).resize((tw*S,th*S),Image.LANCZOS); dr=ImageDraw.Draw(c)
        for gx in range((tx//25+1)*25, tx+tw, 25):
            X=(gx-tx)*S; dr.line([(X,0),(X,c.height)],fill=(120,200,255) if gx%50 else (0,120,255),width=1)
            if gx%50==0: dr.text((X+2,2),str(gx),fill=(0,0,200),font=font)
        for gy in range((ty//25+1)*25, ty+th, 25):
            Y=(gy-ty)*S; dr.line([(0,Y),(c.width,Y)],fill=(120,200,255) if gy%50 else (0,120,255),width=1)
            if gy%50==0: dr.text((2,Y+2),str(gy),fill=(0,0,200),font=font)
        for i,(a,b,cc,d) in enumerate(L):
            if max(a,cc)<tx or min(a,cc)>tx+tw or max(b,d)<ty or min(b,d)>ty+th: continue
            dr.line([((a-tx)*S,(b-ty)*S),((cc-tx)*S,(d-ty)*S)],fill=(0,230,230),width=3)
            mx,my=((a+cc)/2-tx)*S,((b+d)/2-ty)*S; dr.text((mx+4,my-8),str(i),fill=(0,0,0),font=font)
        c.save(f'{BASE}/work/qc/q_{ty}_{tx}.png')
