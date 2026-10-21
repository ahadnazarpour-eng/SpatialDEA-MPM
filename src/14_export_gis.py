"""14_export_gis.py -- Export GeoTIFF / shapefile products.

Part of SpatialDEA-MPM (Nazarpour, submitted to Mathematical Geosciences).
Run from the repository root (see README.md). Set SDEA_BASE to override the base directory.

All rasters are written as int16 on the 30 m UTM 39N grid with value = quantity x scale
(the scale is stored in the GeoTIFF tag 'scale_factor' as 1/scale); nodata = -32768.
"""
import os, sys
BASE = os.environ.get("SDEA_BASE", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, rasterio, pandas as pd, geopandas as gpd

W = BASE + '/work'; OUT = W + '/GIS'; os.makedirs(OUT, exist_ok=True)
ref = rasterio.open(W + '/gee/shz_dem_nasadem.tif')
prof = ref.profile; prof.update(compress='deflate', predictor=2, zlevel=9, dtype='int16', nodata=-32768, count=1)

def w(name, arr, scale):
    a = np.where(np.isfinite(arr), np.round(np.asarray(arr, float) * scale), -32768).clip(-32768, 32767).astype('int16')
    with rasterio.open(f'{OUT}/{name}.tif', 'w', **prof) as o:
        o.write(a, 1); o.update_tags(scale_factor=str(1 / scale))

G = np.load(W + '/geology_px.npz'); E = np.load(W + '/evidence_px.npz'); A = np.load(W + '/alteration_z.npz'); F = np.load(W + '/faults_px.npz')
w('granitoids_from_map', G['gran'], 1); w('hornfels_from_map', G['horn'], 1)
w('dist_granitoid_contact_m', G['dist_gran_contact'], 1)
for k in ['A_phyllic', 'A_argillic', 'A_propylitic', 'A_ferric', 'S_ferric', 'S_hydroxyl', 'S_ferrous']:
    w('alt_z_' + k, A[k], 100)
w('alt_OH_composite_z', E['OH'], 100); w('alt_FE_composite_z', E['FE'], 100)
w('mag_analytic_signal', E['AS'], 1e5); w('mag_tilt_derivative', E['TDR'], 100); w('dist_mag_contact_m', E['dist_mag_contact'], 1)
w('fault_density_km_per_km2_x1000', F['dens'], 1000); w('fault_local_fractal_D_x1000', F['Dloc'], 1000)
w('fault_intersection_density_x1e5', F['inter'], 1e5); w('dist_fault_m', F['dist_fault'], 1)

L = rasterio.open(W + '/dmu_map.tif').read(1); w('dmu_catchment_id', L.astype(float), 1)
R = pd.read_csv(W + '/results_dmu.csv'); g = [c for c in R.columns if c.startswith('GR-DEA (')][0]
for col, nm in [(g, 'score_GRDEA'), ('AR-BoD (alpha=0.05)', 'score_ARBoD'), ('BoD (CCR, unit input)', 'score_BoD')]:
    lut = np.full(L.max() + 1, np.nan); lut[R['index'].values + 1] = R[col].values
    r = lut[L]; r[L == 0] = np.nan; w(nm, r, 10000)

gpd.read_file(W + '/faults_final.gpkg').to_file(OUT + '/faults_final.shp')
pts = gpd.GeoDataFrame(R.rename(columns={g: 'GRDEA', 'AR-BoD (alpha=0.05)': 'ARBoD', 'BoD (CCR, unit input)': 'BoD', 'Equal weights': 'EqualW',
                                         'BLEG Au only': 'Au_only', 'Common weights (lambda->inf)': 'CommonW', 'area_up_km2': 'area_km2'}),
                       geometry=gpd.points_from_xy(R['xs'], R['ys']), crs=32639)
pts = pts[['id', 'au', 'area_km2', 'GRDEA', 'ARBoD', 'BoD', 'EqualW', 'CommonW', 'pos', 'pos_prox', 'ctrl', 'geometry']]
pts[['pos', 'pos_prox', 'ctrl']] = pts[['pos', 'pos_prox', 'ctrl']].astype(int)
pts.to_file(OUT + '/dmu_outlets_scores.shp')
print('GIS products written to', OUT)
