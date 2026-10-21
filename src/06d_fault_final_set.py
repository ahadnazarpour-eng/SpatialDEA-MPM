"""06d_fault_final_set.py -- QC deletions and manually traced polylines.

Part of SpatialDEA-MPM (Nazarpour, submitted to Mathematical Geosciences).
Run from the repository root (see README.md). Set SDEA_BASE to override the base directory.
"""
import os, sys
BASE = os.environ.get("SDEA_BASE", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

"""Final fault set = auto-detected segments (QC-filtered) + manually traced polylines (image px)."""
import numpy as np, json
L=np.load(BASE+'/work/fault_segments_px.npy')
bad=[229,44,265,59,252,168,261,312,105,97,66,292,267,162,258,282,358,176,19,40,306,124,95,323,311,37,71,48,244,164,106,195,342,143,255,191,346,167,43]
auto=[ [(a,b),(c,d)] for i,(a,b,c,d) in enumerate(L) if i not in bad]
manual={
 # tile NW (Jurassic phyllites)
 'nw_a':[(925,123),(980,217)], 'nw_b':[(925,200),(1035,277)], 'nw_d':[(980,60),(1025,137)],
 'nw_e':[(1022,73),(1038,110),(1055,135),(1065,170)], 'nw_f':[(1035,60),(1080,103)],
 'nw_h':[(1098,60),(1212,160)], 'nw_i':[(1135,60),(1200,150)], 'nw_l':[(1098,200),(1255,347)],
 'nw_m':[(1085,200),(1335,290)], 'nw_o':[(1175,253),(1202,293),(1250,347)], 'nw_bb':[(1035,353),(1098,375)],
 'n_k':[(1512,122),(1575,143),(1632,167)], 'n_l':[(1502,72),(1538,123),(1568,193)],
 'n_p':[(1375,267),(1410,343)], 'n_s':[(1382,342),(1435,350),(1475,360),(1538,383)],
 # Sp.sch2 / Astaneh granodiorite
 'c_e':[(1608,535),(1648,568)], 'c_f':[(1695,535),(1735,565)], 'c_i':[(1392,588),(1435,615)], 'c_h':[(1435,575),(1462,642)],
 # east Jurassic block
 'e_41':[(1675,1003),(1700,947)], 'e_345':[(2000,1102),(2055,1148)], 'e_166':[(1390,1042),(1453,1082)],
 # teal Triassic belt
 's_225':[(1585,1170),(1745,1247)], 's_341':[(1535,1197),(1632,1252)], 's_90':[(1635,1243),(1712,1285)], 's_262':[(1510,1270),(1572,1303)],
 # regional sutures
 'MZRF':[(972,1235),(1035,1295),(1098,1328),(1158,1348),(1202,1368),(1250,1390),(1290,1443),(1322,1482)],
 'MRF':[(938,1385),(1035,1492)],
}
json.dump({'auto':[[list(map(float,p)) for p in s] for s in auto],'manual':{k:[list(map(float,p)) for p in v] for k,v in manual.items()}},open(BASE+'/work/faults_final_px.json','w'))
print('auto kept',len(auto),'removed',len(bad),'manual',len(manual))
