import h5py
import rioxarray as rxr
import numpy as np
import xarray as xr

# Path to VIIRS HDF5
path = '/Users/farrospide/Downloads/VNP10A1F.A2022273.h11v12.002.2023058115916.h5'

# Inspect HDF5 structure
with h5py.File(path, 'r') as f:
    print("File structure:")
    f.visit(print)
    
    # Get snow cover data path
    snow_data = f['HDFEOS/GRIDS/VIIRS_Grid_IMG_2D/Data Fields/Daily_NDSI_Snow_Cover'][:]
    
    # Get geolocation information
    x_coords = f['HDFEOS/GRIDS/VIIRS_Grid_IMG_2D/XDim'][:]
    y_coords = f['HDFEOS/GRIDS/VIIRS_Grid_IMG_2D/YDim'][:]

# Create xarray DataArray with proper coordinates
dataset = rxr.open_rasterio(path, 
    group='HDFEOS/GRIDS/VIIRS_Grid_IMG_2D/Data Fields/Daily_NDSI_Snow_Cover',
    x_dim='XDim',
    y_dim='YDim')

# Create DataArray from snow_data with proper dimensions
dataset = xr.DataArray(
    data=snow_data,
    dims=['y', 'x'],
    coords={
        'y': y_coords,
        'x': x_coords
    }
)

# Define VIIRS sinusoidal projection
viirs_sinu = '+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'
dataset=dataset.rio.write_crs(viirs_sinu, inplace=True)
dataset=dataset.rio.reproject('EPSG:4326')

print("Dataset coords:", dataset.coords)
dataset.plot()

dataset.rio.to_raster('/Users/farrospide/Downloads/snow_cover.tif')