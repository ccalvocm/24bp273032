// Define area of interest (example coordinates)
var geometry = ee.Geometry.Rectangle([-71.71782, -32.28247, -69.809361	, -29.0366]);

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
var ndsi = image.normalizedDifference(['I1', 'I3'])
  .rename('NDSI');

// Create snow cover mask (threshold 0.4)
var snowCover = ndsi.gt(0.4);

// Visualization
Map.centerObject(geometry, 8);
//Map.addLayer(ndsi, {min: 0.4, max: 1, palette: ['blue', 'white']}, 'NDSI');
Map.addLayer(snowCover.updateMask(snowCover), {palette: ['cyan']}, 'Snow Cover');

