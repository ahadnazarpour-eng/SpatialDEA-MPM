// gee_processing.js -- Earth Engine part of SpatialDEA-MPM (run in the Code Editor, then press Run on the three tasks).
// Download the three GeoTIFFs from Google Drive into work/gee/.
// In the study the same images were fetched in tiles with getDownloadURL; the result is identical.
// Earth Engine processing used in the study (run in the Code Editor).
var AOI = ee.Geometry.Rectangle([49.0,33.5,49.5,34.0]);
// ---- ASTER L1T -> TOA reflectance, summer median 2000-2007 (pre-SWIR failure)
var B = ['B01','B02','B3N','B04','B05','B06','B07','B08','B09'];
var ESUN = [1845.99,1555.74,1119.47,231.25,79.81,74.99,68.66,59.74,56.92];
var astCol = ee.ImageCollection('ASTER/AST_L1T_003').filterBounds(AOI).filterDate('2000-01-01','2008-04-01')
  .filter(ee.Filter.and(ee.Filter.listContains('ORIGINAL_BANDS_PRESENT','B04'),ee.Filter.listContains('ORIGINAL_BANDS_PRESENT','B09'),ee.Filter.listContains('ORIGINAL_BANDS_PRESENT','B01')))
  .filter(ee.Filter.lt('CLOUDCOVER',5)).filter(ee.Filter.calendarRange(6,10,'month'));
var toRefl = function(img){
  var doy = img.date().getRelative('day','year');
  var d = ee.Number(1).subtract(ee.Number(0.01672).multiply(ee.Number(2*Math.PI).multiply(doy.subtract(4)).divide(365).cos()));
  var cosz = ee.Number(90).subtract(img.get('SOLAR_ELEVATION')).multiply(Math.PI/180).cos();
  var bands = B.map(function(b,i){
    var g = ee.Number(img.get('GAIN_COEFFICIENT_'+b)); var dn = img.select(b);
    return dn.subtract(1).multiply(g).multiply(Math.PI).multiply(d.pow(2)).divide(cosz.multiply(ESUN[i])).updateMask(dn.gt(0)).toFloat().rename(b);});
  return ee.Image.cat(bands).copyProperties(img,['system:time_start']);
};
var aster = astCol.map(toRefl).median().clip(AOI);
// ---- Sentinel-2 L2A summer median 2023-2024
var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED').filterBounds(AOI).filterDate('2023-07-01','2024-10-01')
  .filter(ee.Filter.calendarRange(7,9,'month')).filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',3))
  .map(function(i){var s=i.select('SCL'); return i.updateMask(s.neq(3).and(s.neq(8)).and(s.neq(9)).and(s.neq(10)).and(s.neq(11))).select(['B2','B3','B4','B8','B11','B12']);})
  .median().clip(AOI);
var dem = ee.Image('NASA/NASADEM_HGT/001').select('elevation');
var region = ee.Geometry.Rectangle([314100,3707600,361500,3764000],'EPSG:32639',false);
Export.image.toDrive({image: aster.multiply(10000).round().toInt16().unmask(-32768), description:'shz_aster_tile0', fileNamePrefix:'shz_aster_tile0', region: region, crs:'EPSG:32639', scale:30, maxPixels:1e10});
Export.image.toDrive({image: s2.toInt16().unmask(-32768), description:'shz_s2_tile0', fileNamePrefix:'shz_s2_tile0', region: region, crs:'EPSG:32639', scale:30, maxPixels:1e10});
Export.image.toDrive({image: dem.toInt16(), description:'shz_dem_nasadem', fileNamePrefix:'shz_dem_nasadem', region: region, crs:'EPSG:32639', scale:30, maxPixels:1e10});
// ---- Boundaries for the location map (data/countries_LSIB.geojson, data/iran_provinces_GAUL.geojson)
var box = ee.Geometry.Rectangle([38,22,66,42]);
var ctry = ee.FeatureCollection('USDOS/LSIB_SIMPLE/2017').filterBounds(box)
  .map(function(f){return ee.Feature(f.geometry().intersection(box,1000)).set('name',f.get('country_na'));});
var prov = ee.FeatureCollection('FAO/GAUL_SIMPLIFIED_500m/2015/level1').filter(ee.Filter.eq('ADM0_NAME','Iran  (Islamic Republic of)')).select(['ADM1_NAME']);
Export.table.toDrive({collection: ctry, description:'countries_LSIB', fileFormat:'GeoJSON'});
Export.table.toDrive({collection: prov, description:'iran_provinces_GAUL', fileFormat:'GeoJSON'});
