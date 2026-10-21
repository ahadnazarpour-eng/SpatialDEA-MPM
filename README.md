# SpatialDEA-MPM

**Online Resource 1** for the article

> Nazarpour A, Mehregan F. *A Constrained Spatial Data Envelopment Analysis Framework with Distance-Decay Kernels
> for Multi-Source Mineral Prospectivity Modeling.* Mathematical Geosciences (submitted).
>
> Corresponding author: Ahad Nazarpour, Department of Geology, Ahv.C., Islamic Azad University, Ahvaz, Iran;
> ahad.nazarpour@gmail.com

The workflow builds catchment decision-making units (DMUs) from BLEG stream-sediment samples. It then derives
remote-sensing, magnetic, lithological and structural evidence, analyses the fault network (orientation,
power-law and fractal scaling, spatial association, Fry analysis), and ranks catchments with benefit-of-the-doubt
DEA, assurance-region DEA, super-efficiency and geographically regularised DEA (GR-DEA).

## Computer code availability

| Item | |
|---|---|
| Name of code | SpatialDEA-MPM |
| Developer | Ahad Nazarpour (ahad.nazarpour@gmail.com) |
| Year first available | 2026 |
| Hardware | Any 64-bit PC; ≥ 8 GB RAM (16 GB recommended) |
| Software | Python ≥ 3.11 (tested with 3.13) and the packages in `requirements.txt`; a Google Earth Engine account for `gee/gee_processing.js` |
| Program language | Python, JavaScript (Earth Engine) |
| Program size | ≈ 90 kB of source code |
| Licence | MIT |
| Availability | Online Resource 1 of the article; https://github.com/ahadnazarpour-eng/SpatialDEA-MPM |

## Installation

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python3 tests/test_theory.py      # numerical checks of Propositions 2–3 and Theorem 1 (< 1 min)
```

## Directory layout

```
SpatialDEA-MPM/
├── data/            input data (see data/README.md; GSI data are not redistributed)
├── gee/             Earth Engine script (ASTER, Sentinel-2, NASADEM, boundaries)
├── src/             numbered pipeline scripts + dea.py (the DEA library)
├── tests/           numerical verification of the mathematical results
├── work/            created at run time: intermediate files, figures (work/fig_journal), GIS products (work/GIS)
└── run_all.sh       runs the whole pipeline
```

## Running

1. Put the input files in `data/` (`data/README.md`).
2. Run `gee/gee_processing.js` in the Earth Engine Code Editor and download the exported GeoTIFFs to `work/gee/`.
3. Run `./run_all.sh`. A single script can be run with `SDEA_BASE=$(pwd) python3 src/<script>.py`.

A full run takes about 30–40 minutes on a laptop. Most of the time is spent in step 10, which solves
several hundred quadratic programmes.

## Pipeline and correspondence with the manuscript

| Step | Script | Manuscript |
|---|---|---|
| 01 | `01_magnetic_derivatives.py` – FFT reduction to the pole, analytic signal, tilt derivative | Sect. 2.4, Fig. 4 |
| 02 | `02_mosaic_gee_tiles.py` – mosaic Earth Engine exports on the 30 m grid | Sect. 2.3 |
| 03 | `03_alteration_indices.py` – ASTER / Sentinel-2 band-ratio indices, NDVI mask | Sect. 2.3 |
| 04 | `04_evidence_layers.py` – hydroxyl and ferric composites, magnetic contacts | Sect. 2.3–2.4, Fig. 3 |
| 05a–b | `05a_geomap_lithology_masks.py`, `05b_georeference_geology.py` – granitoids and hornfels from the scanned map | Sect. 2.5 |
| 06a–d | `06a_fault_candidate_mask.py`, `06b_fault_hough_segments.py`, `06c_fault_qc_tiles.py`, `06d_fault_final_set.py` – fault extraction, visual QC, manual traces | Sect. 2.5, Fig. 2a |
| 07a–c | `07a_catchments.py`, `07b_dmu_map_and_labels.py`, `07c_deposit_table.py` – D8 routing, nested catchments, DMU map, deposit links | Sect. 3.1, Fig. 5 |
| 08 | `08_fault_analysis.py` – axial statistics, Rayleigh test, MLE power law, box-counting D with random-line null model, local D, densities, Monte Carlo association, Fry | Sect. 3.6, 4.2, Table 1, Figs. 6–7 |
| 09 | `09_dmu_evidence_table.py` – catchment-mean outputs, Gaussian distance-decay kernels | Sect. 3.2 |
| — | `dea.py` – BoD (Eq. 3), AR-BoD (Eq. 4), super-efficiency, GR-DEA (Eqs. 5–6), common weights, AUC and permutation tests | Sect. 3.3–3.5 |
| 10 | `10_dea_experiments.py` – λ path and L-curve, model comparison, sensitivity, stability | Sect. 4.3–4.5, Tables 2–3, Figs. 8–10, 12 |
| 11 | `11_bootstrap_auc.py` – stratified bootstrap and paired differences | Table 3 |
| 12 | `12_iran_zones_schematic.py` – schematic structural zones for the location map | Fig. 1a |
| 13 | `13_figures.py` – all figures in journal format (174 mm wide, 600 dpi TIFF) | Figs. 1–12 |
| 14 | `14_export_gis.py` – GeoTIFF and shapefile products | Supplementary data |

## Notes on reproducibility

- Random elements (Hough transform, null models, Monte Carlo tests, bootstrap, noise experiments) use fixed seeds.
- The fault set combines automatic detection with a documented visual QC: the list of deleted segment IDs and
  the manually traced polylines (image pixel coordinates) are hard-coded in `06d_fault_final_set.py`.
- The validation labels (DMU indices of target and control catchments) printed by `07b_dmu_map_and_labels.py`
  are hard-coded in `10_dea_experiments.py`.
- Solvers: HiGHS (linear programmes, via SciPy) and Clarabel (quadratic programmes, via CVXPY).

## Licence

MIT (see `LICENSE`).
