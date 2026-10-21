"""13_figures.py -- Journal-format figures for Mathematical Geosciences.

Width 174 mm (single column of the large format), max height 234 mm, Arial-metric font (Liberation Sans)
at 7-8 pt, panel labels (a), (b), ... only (no titles inside figures), 600 dpi TIFF (combination art).
Output: fig_journal/Fig1.tif ... Fig12.tif, plus 300 dpi JPEG copies for embedding in the Word file.
"""
import os, sys, json
import numpy as np, pandas as pd, rasterio, geopandas as gpd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import LightSource, ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, Patch
from matplotlib.gridspec import GridSpec
from scipy.stats import rankdata
from scipy.ndimage import map_coordinates
from PIL import Image
from skimage.transform import ProjectiveTransform
from shapely.geometry import Polygon
from shapely.ops import unary_union

BASE = os.environ.get('SDEA_BASE', os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
W = os.path.join(BASE, 'work'); DATA = os.path.join(BASE, 'data')
OUT = os.path.join(W, 'fig_journal'); os.makedirs(OUT, exist_ok=True)
MM = 1/25.4; FW = 174*MM; HMAX = 234*MM
plt.rcParams.update({'font.family': ['Arial', 'Liberation Sans', 'DejaVu Sans'], 'font.size': 7.5, 'axes.labelsize': 7.5, 'xtick.labelsize': 6.5,
                     'ytick.labelsize': 6.5, 'legend.fontsize': 6.5, 'axes.linewidth': .5, 'xtick.major.width': .5,
                     'ytick.major.width': .5, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'lines.linewidth': 1,
                     'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': True, 'legend.framealpha': .92,
                     'legend.edgecolor': '#999999', 'legend.borderpad': .4, 'legend.handlelength': 1.6})
C = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#7a5af0', '#6b6b6b']

ref = rasterio.open(f'{W}/gee/shz_dem_nasadem.tif'); dem = ref.read(1).astype(float); T = ref.transform
ext = [T.c/1000, (T.c+T.a*ref.width)/1000, (T.f+T.e*ref.height)/1000, T.f/1000]
hs = LightSource(315, 45).hillshade(dem, vert_exag=2, dx=30, dy=30)
dep = gpd.read_file(f'{DATA}/Deposits.shp').to_crs(32639)
dep['grp'] = np.where(dep['type'].isin(['Cu', 'W', 'Au_Vein']), 'W–Cu–Au (target)',
              np.where(dep['type'] == 'Talc', 'Talc (control)', np.where(dep['type'] == 'Si', 'Silica (control)', 'Feldspar')))
S = pd.read_csv(f'{W}/samples_snapped.csv')
R = pd.read_csv(f'{W}/results_dmu.csv'); meta = json.load(open(f'{W}/meta.json'))
grcol = [c for c in R.columns if c.startswith('GR-DEA (')][0]
Fz = gpd.read_file(f'{W}/faults_final.gpkg')
G = np.load(f'{W}/geology_px.npz'); L = rasterio.open(f'{W}/dmu_map.tif').read(1)
tg = dep[dep.grp.str.startswith('W')]
STY = {'W–Cu–Au (target)': ('*', C[1], 45), 'Talc (control)': ('s', C[0], 10), 'Silica (control)': ('^', C[2], 12), 'Feldspar': ('D', C[7], 9)}

def panel(ax, letter, x=-0.02, y=1.015):
    ax.text(x, y, f'({letter})', transform=ax.transAxes, fontsize=8.5, fontweight='bold', ha='right', va='bottom')

def base(ax, alpha=.9, xlabel=True, ylabel=True):
    ax.imshow(hs, cmap='gray', extent=ext, alpha=alpha, vmin=0, vmax=1, interpolation='bilinear')
    ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3]); ax.set_aspect('equal')
    if xlabel: ax.set_xlabel('Easting (km)')
    if ylabel: ax.set_ylabel('Northing (km)')
    ax.set_xticks([320, 330, 340, 350, 360])

def deps(ax, legend=True, loc='lower left', scale=1.0):
    for g, (mk, col, sz) in STY.items():
        q = dep[dep.grp == g]; ax.scatter(q.geometry.x/1000, q.geometry.y/1000, marker=mk, s=sz*scale, c=col, edgecolor='k', lw=.35, zorder=6, label=g)
    if legend: ax.legend(loc=loc)

def score_raster(col):
    lut = np.full(L.max()+1, np.nan); lut[R['index'].values+1] = R[col].values
    out = lut[L]; out[L == 0] = np.nan; return out

def cbar(fig, im, ax, label, **kw):
    cb = fig.colorbar(im, ax=ax, shrink=kw.pop('shrink', .8), pad=kw.pop('pad', .02), aspect=kw.pop('aspect', 25), **kw)
    cb.set_label(label); cb.outline.set_linewidth(.4); cb.ax.tick_params(width=.4, length=2); return cb

def save(fig, n):
    assert fig.get_size_inches()[0] <= FW+1e-3 and fig.get_size_inches()[1] <= HMAX+1e-3, fig.get_size_inches()
    fig.savefig(f'{OUT}/Fig{n}.tif', dpi=600, pil_kwargs={'compression': 'tiff_lzw'})
    fig.savefig(f'{OUT}/Fig{n}_300dpi.jpg', dpi=300, pil_kwargs={'quality': 92})
    plt.close(fig); print('Fig', n)

# ------------------------------------------------------------- Fig 1
Cn = gpd.read_file(f'{DATA}/countries_LSIB.geojson'); Pv = gpd.read_file(f'{DATA}/iran_provinces_GAUL.geojson')
Z = gpd.read_file(f'{W}/zones_schematic.geojson')
fig = plt.figure(figsize=(FW, 3.55)); gs = GridSpec(1, 2, width_ratios=[1.12, 1], wspace=.22, left=.07, right=.99, bottom=.11, top=.95)
ax = fig.add_subplot(gs[0]); ax.set_facecolor('#d8e8f5')
Cn[Cn.name != 'Iran'].plot(ax=ax, color='#efefea', edgecolor='#9a9a9a', lw=.3)
iran = Cn[Cn.name == 'Iran']; iran.plot(ax=ax, color='#fbfaf5', edgecolor='k', lw=.6)
Pv.boundary.plot(ax=ax, color='#b5b5b5', lw=.25)
zmt = Z[Z.name == 'Main Zagros Thrust'].geometry.iloc[0]; sas = Z[Z.name == 'SaSZ–UDMA boundary'].geometry.iloc[0]; ud = Z[Z.name == 'UDMA NE boundary'].geometry.iloc[0]
band = lambda l1, l2: Polygon(list(l1.coords)+list(l2.coords)[::-1])
ir = unary_union([g.buffer(0) for g in iran.geometry])
gpd.GeoSeries([band(zmt, sas).buffer(0).intersection(ir)]).plot(ax=ax, color='#f4b183', alpha=.8, lw=0)
gpd.GeoSeries([band(sas, ud).buffer(0).intersection(ir)]).plot(ax=ax, color='#c5a3d9', alpha=.8, lw=0)
for nm, ls, col, lw in [('Main Zagros Thrust', '-', '#c0392b', 1.1), ('SaSZ–UDMA boundary', '--', '#555', .6), ('UDMA NE boundary', '--', '#555', .6)]:
    g = Z[Z.name == nm].geometry.iloc[0].intersection(ir.buffer(.05))
    for p in getattr(g, 'geoms', [g]): x, y = p.xy; ax.plot(x, y, ls, c=col, lw=lw)
ax.add_patch(Rectangle((49, 33.5), .5, .5, fc='#d62728', ec='k', lw=.6, zorder=8))
ax.annotate('Shazand sheet', xy=(49.3, 33.85), xytext=(54.6, 38.9), fontsize=6.5, arrowprops=dict(arrowstyle='->', lw=.6), zorder=9)
for t, x, y, r in [('Zagros Fold–Thrust Belt', 50.2, 30.4, -38), ('Sanandaj–Sirjan Zone', 52.0, 32.3, -38), ('Urumieh–Dokhtar Arc', 53.5, 33.5, -38),
                   ('Central Iran', 56.8, 33.2, 0), ('Alborz', 52.5, 36.4, 0), ('Kopet Dagh', 58.6, 37.6, 0), ('Makran', 60.5, 26.4, 0)]:
    ax.text(x, y, t, fontsize=5.5, rotation=r, ha='center', va='center', color='#333')
for t, x, y in [('Caspian Sea', 51.0, 39.6), ('Persian Gulf', 51.0, 27.2), ('Gulf of Oman', 58.6, 24.7)]:
    ax.text(x, y, t, fontsize=5.5, color='#2a5d8a', ha='center')
ax.set_xlim(43.5, 63.8); ax.set_ylim(24.2, 40.3); ax.set_aspect(1/np.cos(np.radians(32)))
ax.set_xlabel('Longitude (°E)'); ax.set_ylabel('Latitude (°N)')
ax.legend(handles=[Patch(fc='#f4b183', label='Sanandaj–Sirjan Zone'), Patch(fc='#c5a3d9', label='Urumieh–Dokhtar Magmatic Arc'),
                   Line2D([], [], c='#c0392b', lw=1.1, label='Main Zagros Thrust'), Rectangle((0, 0), 1, 1, fc='#d62728', ec='k', lw=.5, label='Study area')],
          loc='lower left', fontsize=5.8)
panel(ax, 'a')
img = np.asarray(Image.open(f'{DATA}/Shazand_geomap.jpg').convert('RGB'))
srcp = np.array([[927, 64], [2142, 61], [933, 1518], [2151, 1514]], float); dstp = np.array([[49.0, 34.0], [49.5, 34.0], [49.0, 33.5], [49.5, 33.5]])
nx = ny = 1500; gx = 49.0+(np.arange(nx)+.5)/nx*.5; gy = 34.0-(np.arange(ny)+.5)/ny*.5
pt = ProjectiveTransform.from_estimate(dstp, srcp); GX, GY = np.meshgrid(gx, gy); pix = pt(np.c_[GX.ravel(), GY.ravel()])
war = np.stack([map_coordinates(img[..., k], [pix[:, 1], pix[:, 0]], order=1).reshape(ny, nx) for k in range(3)], -1).astype('uint8')
ax = fig.add_subplot(gs[1]); ax.imshow(war, extent=[49, 49.5, 33.5, 34], interpolation='lanczos'); ax.set_aspect(1/np.cos(np.radians(33.75)))
dll = dep.to_crs(4326)
for g, (mk, col, sz) in STY.items():
    q = dll[dll.grp == g]; ax.scatter(q.geometry.x, q.geometry.y, marker=mk, s=sz, c=col, edgecolor='k', lw=.35, zorder=6, label=g)
ax.legend(loc='lower right', fontsize=5.8)
ax.set_xticks([49, 49.125, 49.25, 49.375, 49.5]); ax.set_xticklabels(["49°00′", "49°07.5′", "49°15′", "49°22.5′", "49°30′"])
ax.set_yticks([33.5, 33.625, 33.75, 33.875, 34]); ax.set_yticklabels(["33°30′", "33°37.5′", "33°45′", "33°52.5′", "34°00′"])
panel(ax, 'b'); save(fig, 1)

# ------------------------------------------------------------- Fig 2
fig, ax = plt.subplots(1, 2, figsize=(FW, 4.0), gridspec_kw={'wspace': .12, 'left': .075, 'right': .88, 'bottom': .09, 'top': .96})
base(ax[0]); ov = np.zeros(dem.shape+(4,)); ov[G['horn']] = [.55, .35, .2, .55]; ov[G['gran']] = [.85, .2, .15, .7]
ax[0].imshow(ov, extent=ext)
for g in Fz.geometry: x, y = g.xy; ax[0].plot(np.array(x)/1000, np.array(y)/1000, c='#b0306a', lw=.45)
h = [Rectangle((0, 0), 1, 1, color=(.85, .2, .15, .7)), Rectangle((0, 0), 1, 1, color=(.55, .35, .2, .55)), Line2D([], [], c='#b0306a', lw=.8)]
l1 = ax[0].legend(h, ['Granitoids', 'Hornfels aureoles', f'Fault traces (n = {len(Fz)})'], loc='upper right', fontsize=5.8); ax[0].add_artist(l1)
deps(ax[0], loc='lower left'); panel(ax[0], 'a')
base(ax[1], ylabel=False)
sc = ax[1].scatter(S.x/1000, S.y/1000, c=np.log10(np.where(S.au <= .1, .05, S.au)), cmap='Blues', s=7, edgecolor='k', lw=.2, vmin=-1.3, vmax=1, zorder=5)
deps(ax[1], False)
cax = ax[1].inset_axes([1.03, .15, .035, .7]); cb = fig.colorbar(sc, cax=cax); cb.set_label('log$_{10}$ Au (ppb)'); cb.outline.set_linewidth(.4)
panel(ax[1], 'b'); save(fig, 2)

# ------------------------------------------------------------- Fig 3
E = np.load(f'{W}/evidence_px.npz')
fig, ax = plt.subplots(1, 2, figsize=(FW, 4.15), gridspec_kw={'wspace': .12, 'left': .075, 'right': .93, 'bottom': .09, 'top': .96})
for i, (a, k) in enumerate(zip(ax, ['OH', 'FE'])):
    base(a, ylabel=(i == 0)); v = np.where(E['mask'], E[k], np.nan)
    im = a.imshow(np.where(v > 1.0, v, np.nan), cmap='OrRd', vmin=1, vmax=3, extent=ext, interpolation='nearest')
    deps(a, i == 0); panel(a, 'ab'[i])
cbar(fig, im, list(ax), 'Composite z-score', shrink=.7); save(fig, 3)

# ------------------------------------------------------------- Fig 4 (2 + 1)
mg = np.load(f'{W}/mag_grids.npz'); e2 = [48.995, 49.495, 33.495, 33.995]
fig = plt.figure(figsize=(FW, 6.9)); gs = GridSpec(2, 4, hspace=.28, wspace=.55, left=.07, right=.97, bottom=.06, top=.97)
axs = [fig.add_subplot(gs[0, 0:2]), fig.add_subplot(gs[0, 2:4]), fig.add_subplot(gs[1, 1:3])]
for i, (a, (k, cm, lab)) in enumerate(zip(axs, [('rtp', 'RdBu_r', 'nT'), ('AS', 'Blues', 'nT m$^{-1}$'), ('TDR', 'RdBu_r', 'rad')])):
    im = a.imshow(mg[k], extent=e2, cmap=cm, interpolation='bicubic'); cbar(fig, im, a, lab, shrink=.85)
    if k == 'TDR': a.contour(np.linspace(49, 49.49, 50), np.linspace(33.99, 33.5, 50), mg['TDR'], levels=[0], colors='k', linewidths=.5)
    a.scatter(tg.X, tg.Y, marker='*', s=40, c=C[1], edgecolor='k', lw=.35, zorder=5, label='W–Cu–Au occurrences')
    a.set_xlabel('Longitude (°E)'); a.set_ylabel('Latitude (°N)'); a.set_aspect(1/np.cos(np.radians(33.75))); panel(a, 'abc'[i])
axs[0].legend(loc='lower left'); save(fig, 4)

# ------------------------------------------------------------- Fig 5
fig, ax = plt.subplots(figsize=(4.6, 5.3)); fig.subplots_adjust(left=.13, right=.98, bottom=.08, top=.98); base(ax)
rng = np.random.default_rng(3); lut = rng.random(L.max()+1); cat = np.where(L > 0, lut[L], np.nan)
ax.imshow(cat, cmap=ListedColormap(plt.cm.tab20(np.linspace(0, 1, 20))), alpha=.45, extent=ext, interpolation='nearest')
ax.scatter(S['xs']/1000, S['ys']/1000, s=3, c='k', label='Sample outlet', zorder=5); deps(ax); save(fig, 5)

# ------------------------------------------------------------- Fig 6 (2 x 2)
fs = json.load(open(f'{W}/fault_stats.json')); dc = json.load(open(f'{W}/dist_curves.json'))
fig = plt.figure(figsize=(FW, 6.4)); gs = GridSpec(2, 2, hspace=.38, wspace=.3, left=.08, right=.98, bottom=.08, top=.95)
a = fig.add_subplot(gs[0, 0], projection='polar')
bins = np.radians(np.arange(0, 361, 10)); w = Fz.length_m.values/1000
hh, _ = np.histogram(np.radians(np.r_[Fz.azimuth.values, Fz.azimuth.values+180]), bins=bins, weights=np.r_[w, w])
a.bar(bins[:-1]+np.radians(5), hh, width=np.radians(10), color=C[0], edgecolor='white', lw=.4)
a.set_theta_zero_location('N'); a.set_theta_direction(-1); a.set_yticklabels([]); a.grid(lw=.3); a.spines['polar'].set_linewidth(.4)
a.text(-0.05, 1.04, '(a)', transform=a.transAxes, fontsize=8.5, fontweight='bold')
a = fig.add_subplot(gs[0, 1]); Ls = np.sort(Fz.length_m.values)[::-1]/1000; Nc = np.arange(1, len(Ls)+1)
a.loglog(Ls, Nc, 'o', ms=2.2, c=C[0], mfc='none', mew=.5, label='Fault traces')
pl = fs['powerlaw']; xx = np.logspace(np.log10(pl['Lmin_m']/1000), np.log10(Ls.max()), 50)
a.loglog(xx, pl['n_tail']*(xx/(pl['Lmin_m']/1000))**(1-pl['alpha']), c=C[1], lw=1.2, label=f"Power law, α = {pl['alpha']:.2f} ± {pl['alpha_se']:.2f}")
a.axvline(pl['Lmin_m']/1000, c='#999', ls='--', lw=.6); a.set_xlabel('Trace length L (km)'); a.set_ylabel('Cumulative number N(≥ L)'); a.legend(); panel(a, 'b')
a = fig.add_subplot(gs[1, 0]); bc = fs['boxcount']; r_ = np.array(bc['sizes_m'])/1000; Nb = np.array(bc['N'])
a.loglog(1/r_, Nb, 's', ms=3.5, c=C[0], label='Fault network')
cf = np.polyfit(np.log(1/r_), np.log(Nb), 1); a.loglog(1/r_, np.exp(np.polyval(cf, np.log(1/r_))), c=C[0], lw=1, label=f"D = {bc['D']:.2f} ± {bc['D_se']:.2f}")
a.loglog(1/r_, Nb[0]*(r_/r_[0])**(-fs['boxcount_null']['mean']), c=C[7], ls='--', lw=.9, label=f"Random lines, D = {fs['boxcount_null']['mean']:.2f}")
a.set_xticks(1/r_); a.set_xticklabels([f'{v:g}' for v in 1/r_]); a.minorticks_off()
a.set_xlabel('1/r (km$^{-1}$)'); a.set_ylabel('Occupied boxes N(r)'); a.legend(); panel(a, 'c')
a = fig.add_subplot(gs[1, 1]); dd = np.array(dc['d_m'])/1000
a.plot(dd, np.array(dc['area'])*100, c=C[7], lw=1.2, label='Sheet area')
a.step(dd, np.array(dc['target'])*100, where='post', c=C[1], lw=1.2, label='W–Cu–Au (n = 10)')
a.step(dd, np.array(dc['control'])*100, where='post', c=C[0], lw=1.2, label='Talc + silica (n = 14)')
a.set_xlabel('Distance to nearest fault (km)'); a.set_ylabel('Cumulative percentage'); a.legend(loc='lower right'); panel(a, 'd'); save(fig, 6)

# ------------------------------------------------------------- Fig 7 (2 + 1)
FP = np.load(f'{W}/faults_px.npz'); FL = np.load(f'{W}/fractal_local.npz'); fry = np.load(f'{W}/fry.npy')
fig = plt.figure(figsize=(FW, 7.6)); gs = GridSpec(2, 4, height_ratios=[1.15, 1], hspace=.22, wspace=.5, left=.07, right=.95, bottom=.06, top=.97)
a1 = fig.add_subplot(gs[0, 0:2]); a2 = fig.add_subplot(gs[0, 2:4]); a3 = fig.add_subplot(gs[1, 1:3])
base(a1); im = a1.imshow(np.where(FP['dens'] > .02, FP['dens'], np.nan), cmap='Purples', extent=ext, vmin=0, vmax=1.2, alpha=.85)
for g in Fz.geometry: x, y = g.xy; a1.plot(np.array(x)/1000, np.array(y)/1000, c='k', lw=.25)
cbar(fig, im, a1, 'Fault density (km km$^{-2}$)', shrink=.75); deps(a1); panel(a1, 'a')
base(a2, ylabel=False); xs, ys = FL['xs']/1000, FL['ys']/1000; e3 = [xs[0]-.5, xs[-1]+.5, ys[-1]-.5, ys[0]+.5]
a2.imshow(np.where(np.isnan(FL['D']), 1, np.nan), extent=e3, cmap=ListedColormap(['#d9d9d9']), alpha=.6, interpolation='nearest')
im = a2.imshow(FL['D'], extent=e3, cmap='viridis', vmin=0.6, vmax=1.4, alpha=.85, interpolation='bilinear')
for g in Fz.geometry: x, y = g.xy; a2.plot(np.array(x)/1000, np.array(y)/1000, c='w', lw=.2)
cbar(fig, im, a2, 'Local dimension D', shrink=.75); deps(a2, False)
a2.legend(handles=[Patch(fc='#d9d9d9', label='D not estimated')], loc='lower left'); panel(a2, 'b')
a3.scatter(fry[:, 0]/1000, fry[:, 1]/1000, s=6, c=C[1], edgecolor='k', lw=.2); a3.axhline(0, c='#ccc', lw=.5); a3.axvline(0, c='#ccc', lw=.5)
fr = fs['fry']; th = np.radians(fr['mean_trend'])
a3.plot([-30*np.sin(th), 30*np.sin(th)], [-30*np.cos(th), 30*np.cos(th)], c=C[0], lw=.9, ls='--', label=f"Mean trend {fr['mean_trend']:.0f}° (R = {fr['R']:.2f})")
a3.set_aspect('equal'); a3.set_xlim(-32, 32); a3.set_ylim(-32, 32); a3.set_xlabel('ΔE (km)'); a3.set_ylabel('ΔN (km)'); a3.legend(loc='lower right'); panel(a3, 'c')
save(fig, 7)

# ------------------------------------------------------------- Fig 8 (2 x 2)
cols = ['BLEG Au only', 'BoD (CCR, unit input)', 'AR-BoD (alpha=0.05)', grcol]
fig, ax = plt.subplots(2, 2, figsize=(FW, 7.2), gridspec_kw={'hspace': .04, 'wspace': .06, 'left': .07, 'right': .9, 'bottom': .06, 'top': .97})
for i, (a, c) in enumerate(zip(ax.flat, cols)):
    base(a, xlabel=(i >= 2), ylabel=(i % 2 == 0))
    if i < 2: a.set_xticklabels([])
    if i % 2: a.set_yticklabels([])
    r = score_raster(c); pct = np.full_like(r, np.nan); m = ~np.isnan(r); pct[m] = rankdata(r[m])/m.sum()*100
    im = a.imshow(pct, cmap='YlOrRd', extent=ext, alpha=.85, interpolation='nearest', vmin=0, vmax=100); deps(a, i == 0); panel(a, 'abcd'[i])
cbar(fig, im, list(ax.flat), 'Score percentile (area-weighted)', shrink=.45); save(fig, 8)

# ------------------------------------------------------------- Fig 9 prediction rate
def pr_curve(col):
    r = score_raster(col); m = ~np.isnan(r); vals = np.sort(r[m])[::-1]; dv = []
    for g in tg.geometry:
        rr, cc = ref.index(g.x, g.y); w_ = 33; sub = r[max(rr-w_, 0):rr+w_+1, max(cc-w_, 0):cc+w_+1]
        dv.append(np.nanmax(sub) if np.isfinite(sub).any() else np.nan)
    dv = np.array(dv); ok = np.isfinite(dv); fa = np.array([(vals >= v).mean() for v in dv[ok]])
    xs_ = np.linspace(0, 1, 401); return xs_, [(fa <= x).mean() for x in xs_], ok.sum()
fig, ax = plt.subplots(figsize=(3.6, 3.0)); fig.subplots_adjust(left=.16, right=.97, bottom=.15, top=.97)
for c, col, lab in zip(['BLEG Au only', 'Equal weights', 'BoD (CCR, unit input)', 'AR-BoD (alpha=0.05)', grcol], [C[7], C[3], C[2], C[0], C[1]],
                       ['BLEG Au only', 'Equal weights', 'BoD-DEA', 'AR-BoD', 'GR-DEA (λ* = 30)']):
    xs_, ys_, n = pr_curve(c); ax.plot(np.array(xs_)*100, np.array(ys_)*100, lw=1.2, c=col, label=lab)
ax.plot([0, 100], [0, 100], c='#bbb', lw=.6, ls='--'); ax.set_xlabel('Cumulative mapped area (%)'); ax.set_ylabel('Target occurrences captured (%)'); ax.legend(loc='lower right')
save(fig, 9)

# ------------------------------------------------------------- Fig 10 lambda path (2 + 1)
P = pd.read_csv(f'{W}/lambda_path.csv'); Q = P[P.lam > 0]
fig = plt.figure(figsize=(FW, 5.4)); gs = GridSpec(2, 4, hspace=.45, wspace=.75, left=.08, right=.98, bottom=.09, top=.96)
a0 = fig.add_subplot(gs[0, 0:2]); a1 = fig.add_subplot(gs[0, 2:4]); a2 = fig.add_subplot(gs[1, 1:3])
a0.plot(Q.energy, Q.sum_theta, '-o', c=C[0], lw=1, ms=2.5)
for _, r in Q.iterrows(): a0.annotate(f"{r.lam:g}", (r.energy, r.sum_theta), fontsize=5.5, xytext=(2, 2), textcoords='offset points', color='#444')
s_ = Q[Q.lam == meta['lam_star']]; a0.scatter(s_.energy, s_.sum_theta, s=45, facecolor='none', edgecolor=C[1], lw=1, zorder=5)
a0.set_xscale('log'); a0.set_xlabel('Dirichlet energy $\\mathcal{E}(U)$'); a0.set_ylabel('Total efficiency $\\sum_j\\theta_j$'); panel(a0, 'a')
a1.plot(Q.lam, Q.sum_theta/Q.sum_theta.max(), '-o', c=C[0], lw=1, ms=2.5, label='Total efficiency')
a1.plot(Q.lam, Q.energy/Q.energy.max(), '-s', c=C[1], lw=1, ms=2.5, label='Dirichlet energy')
a1.set_xscale('log'); a1.set_xlabel('λ'); a1.set_ylabel('Normalized value'); a1.legend(); panel(a1, 'b')
a2.plot(Q.lam, Q.auc, '-o', c=C[0], lw=1, ms=2.5, label='Strict labels'); a2.plot(Q.lam, Q.auc_prox, '-s', c=C[2], lw=1, ms=2.5, label='Proximity labels')
a2.axvline(meta['lam_star'], c='#999', ls='--', lw=.6); a2.set_xscale('log'); a2.set_xlabel('λ'); a2.set_ylabel('AUC'); a2.legend(); panel(a2, 'c')
save(fig, 10)

# ------------------------------------------------------------- Fig 11 virtual weights (2 + 2 + 1)
wc = ['GR_w_Au', 'GR_w_OH', 'GR_w_FE', 'GR_w_GC', 'GR_w_MC']
fig = plt.figure(figsize=(FW, 9.1)); gs = GridSpec(3, 4, hspace=.12, wspace=.35, left=.08, right=.88, bottom=.05, top=.975)
pos_ = [gs[0, 0:2], gs[0, 2:4], gs[1, 0:2], gs[1, 2:4], gs[2, 1:3]]
tot = R[wc].sum(axis=1); axs = []
for i, (p_, c) in enumerate(zip(pos_, wc)):
    a = fig.add_subplot(p_); axs.append(a); R['_share'] = R[c]/tot; r = score_raster('_share')
    base(a, xlabel=(i >= 3), ylabel=(i in (0, 2, 4)))
    im = a.imshow(r, cmap='Blues', vmin=0, vmax=.8, extent=ext, interpolation='nearest', alpha=.9); panel(a, 'abcde'[i])
cbar(fig, im, axs, 'Virtual-weight share', shrink=.4); save(fig, 11)

# ------------------------------------------------------------- Fig 12 sensitivity
Sd = pd.read_csv(f'{W}/sensitivity.csv'); Sd = Sd[Sd.t == 1.5]
fig, ax = plt.subplots(1, 2, figsize=(FW, 2.4), gridspec_kw={'wspace': .3, 'left': .08, 'right': .97, 'bottom': .2, 'top': .9})
for i, (a, col) in enumerate(zip(ax, ['AUC_GR', 'AUCprox_GR'])):
    piv = Sd.pivot(index='alpha', columns='h', values=col); im = a.imshow(piv.values, cmap='Blues', vmin=.55, vmax=.85, aspect='auto')
    a.set_xticks(range(len(piv.columns))); a.set_xticklabels(piv.columns); a.set_yticks(range(len(piv.index))); a.set_yticklabels(piv.index)
    a.set_xlabel('Kernel bandwidth h (m)'); a.set_ylabel('Virtual-weight floor α')
    for ii in range(piv.shape[0]):
        for jj in range(piv.shape[1]): a.text(jj, ii, f'{piv.values[ii, jj]:.2f}', ha='center', va='center', fontsize=6.5, color='k' if piv.values[ii, jj] < .75 else 'w')
    cbar(fig, im, a, 'AUC', shrink=1); panel(a, 'ab'[i])
save(fig, 12)
print('done')
