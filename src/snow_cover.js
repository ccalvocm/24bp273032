var geometry = ee.Geometry.Rectangle([-71.71782, -32.28247, -69.809361, -29.0366]);

// Load Landsat 8 TOA imagery for the same date range and area
var landsat8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_TOA')
  .filterDate('2022-09-01', '2022-10-31')
  .filterBounds(geometry)
  .median()
  .clip(geometry);

// Visualization parameters for Landsat 8 (RGB: B4, B3, B2)
var landsatVis = {bands: ['B4', 'B3', 'B2'], min: 0, max: 0.3};

// Add Landsat 8 basemap layer
Map.centerObject(geometry, 8);
Map.addLayer(landsat8, landsatVis, 'Landsat 8 Basemap');

// Load VIIRS surface reflectance
var viirs = ee.ImageCollection('NOAA/VIIRS/001/VNP09GA')
  .filterDate('2022-09-01', '2022-10-31')
  .filterBounds(geometry);

// Create composite
var image = viirs.median().clip(geometry);

// Calculate NDSI using I1 (visible) and I3 (SWIR)
var ndsi = image.normalizedDifference(['I1', 'I3']).rename('NDSI');

// Use the QF2 band (cast to integer) and check Bit 5 (value 32)
var qf2 = image.select('QF2').toInt();
var snowFlag2 = qf2.bitwiseAnd(32).gt(0); // Bit 5 flag from QF2

// Use the QF7 band (cast to integer) and check Bit 0 (value 1) and Bit 5 (value 32)
var qf7 = image.select('QF7').toInt();
var snowFlag7 = qf7.bitwiseAnd(1).gt(0)         // Bit 0: snow present
                  .or(qf7.bitwiseAnd(32).gt(0));   // Bit 5: snow/ice

// Combine the flags from QF2 and QF7
var combinedSnowFlag = snowFlag2;

// Create combined snow cover mask using NDSI threshold and flags
var snowCover = ndsi.gt(0).and(combinedSnowFlag);

// Visualization: add the snow cover mask layer
Map.addLayer(snowCover.updateMask(snowCover), {palette: ['cyan']}, 'Snow Cover');