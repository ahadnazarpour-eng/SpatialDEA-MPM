#!/usr/bin/env bash
# Reproduce the full SpatialDEA-MPM workflow (after placing inputs in data/ and Earth Engine exports in work/gee/).
set -euo pipefail
cd "$(dirname "$0")"
export SDEA_BASE="$(pwd)"
mkdir -p work
for s in 01_magnetic_derivatives 02_mosaic_gee_tiles 03_alteration_indices 04_evidence_layers \
         05a_geomap_lithology_masks 05b_georeference_geology \
         06a_fault_candidate_mask 06b_fault_hough_segments 06c_fault_qc_tiles 06d_fault_final_set \
         07a_catchments 07b_dmu_map_and_labels 07c_deposit_table \
         08_fault_analysis 09_dmu_evidence_table 10_dea_experiments 11_bootstrap_auc \
         12_iran_zones_schematic 13_figures 14_export_gis; do
  echo "=== $s ==="; python3 -I "src/$s.py"
done
echo "All steps finished. Results in work/, figures in work/fig_journal/, GIS products in work/GIS/."
