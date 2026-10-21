# Input data (place these files here before running)

| File | Content | Source |
|---|---|---|
| `Shazand_BLEG.xlsx` | 295 BLEG samples: No., Sample ID, X (UTM 39N), Y (UTM 39N), Au (ppb) | Geological Survey and Mineral Exploration of Iran (Sobhani Dargah 2013) |
| `Deposits.shp` (+ .dbf .shx .prj) | 25 mineral occurrences of the Shazand sheet (fields `NAME`, `type`, `X`, `Y`) | GSI mineral occurrence database |
| `Mag.shp` (+ .dbf .shx .prj) | Aeromagnetic grid points, 0.01° spacing (field `GRID_CODE`, nT) | National aeromagnetic compilation of Iran |
| `Shazand_geomap.jpg` | Scanned 1:100,000 Shazand geological map, sheet 5857 (2576 × 1781 px) | Sahandi et al. (2006), GSI |
| `countries_LSIB.geojson` | Country boundaries (USDOS/LSIB_SIMPLE/2017), exported from Earth Engine | US Department of State |
| `iran_provinces_GAUL.geojson` | Iranian provinces (FAO/GAUL_SIMPLIFIED_500m/2015/level1) | FAO |

The GSI data are not redistributed with the code; request them from the Geological Survey and Mineral
Exploration of Iran or from the author.

Earth Engine exports (produced by `gee/gee_processing.js`) go in `work/gee/`:
`shz_dem_nasadem.tif` (NASADEM, 30 m, EPSG:32639, extent 314100–361500 E, 3707600–3764000 N),
`shz_aster_tile*.tif` (9-band ASTER TOA reflectance × 10000) and `shz_s2_tile*.tif`
(Sentinel-2 L2A B2, B3, B4, B8, B11, B12 × 10000). One un-tiled export per product also works if it is named
`shz_aster_tile0.tif` / `shz_s2_tile0.tif`.
