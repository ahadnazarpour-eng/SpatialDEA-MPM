import numpy as np, gzip, io, sys
def load(p):
    b = open(p,'rb').read()
    if p.endswith('.gz'): b = gzip.decompress(b)
    a = np.load(io.BytesIO(b))
    return a
if __name__=='__main__':
    for p in sys.argv[1:]:
        a=load(p); print(p, a.shape, a.dtype.names)
        for n in a.dtype.names:
            v=a[n].astype(float); m=v>-32768; print('  ',n, m.mean().round(3), np.percentile(v[m],[1,50,99]) if m.any() else None)
