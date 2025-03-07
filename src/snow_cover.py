from earthaccess import Store
import earthaccess
from io import BytesIO
import rioxarray
import h5py
from rioxarray.merge import merge_arrays
import numpy as np
import geopandas as gpd
import pandas as pd
import os
import warnings
from scipy.interpolate import NearestNDInterpolator
from scipy.interpolate import LinearNDInterpolator
import xarray as xr 

# Suppress all warnings
warnings.filterwarnings("ignore")

def get_coords(data_bytes):
    with h5py.File(data_bytes, 'r') as f:
        x_coords = f['HDFEOS/GRIDS/VIIRS_Grid_IMG_2D/XDim'][:]
        y_coords = f['HDFEOS/GRIDS/VIIRS_Grid_IMG_2D/YDim'][:]
    return x_coords, y_coords

def authy():
    # Initialize authentication properly
    auth = earthaccess.login()
    store = Store(auth=auth)
    session = auth.get_session()
    return session

def process_results(results, session):
    sinusoidal_crs = (
        "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 "
        "+a=6371007.181 +b=6371007.181 +units=m +no_defs"
    )
    ds_list = []
    for granule in results:
        try:
            links = granule.data_links()
            if not links:
                continue
            url = links[0]
            response = session.get(url, stream=True)
            response.raise_for_status()
            content = response.content
            # Use memoryview to avoid extra copying
            mem = memoryview(content)
            # Create a single BytesIO instance from the memoryview
            data_bytes = BytesIO(mem)
            # Open the dataset
            ds = rioxarray.open_rasterio(data_bytes, engine='h5netcdf')
            ds = ds.HDFEOS_GRIDS_VIIRS_Grid_IMG_2D_Data_Fields_NDSI_Snow_Cover
            # ds.HDFEOS_GRIDS_VIIRS_Grid_IMG_2D_Data_Fields_Basic_QA

            # Rewind to extract coordinates without creating a new BytesIO object
            data_bytes.seek(0)
            x_dim, y_dim = get_coords(data_bytes)
            ds = ds.assign_coords(x=x_dim, y=y_dim)
            ds = ds.rio.write_crs(sinusoidal_crs, inplace=False)
            ds_list.append(ds)
        except Exception as e:
            print(f"Error processing {url}: {e}")
    return ds_list

def mosaic(ds_list):
    mosaic_da = merge_arrays(ds_list)
    mosaic_utm = mosaic_da.rio.reproject("EPSG:32719")
    return mosaic_utm

def plot(mosaic_utm):
    mosaic_utm.data = mosaic_utm.data.astype(float)
    mosaic_utm.data[mosaic_utm.data > 1] = np.nan
    mosaic_utm.plot(vmin=0, vmax=1)
    return None

def load_basins():
    path=os.path.join('..','geodata',
    'Cuencas_DARH_Coquimbo.gpkg')
    gdf=gpd.read_file(path)
    return gdf

def clip_image(ds, basins, out_path=None):
    """
    Clip the input raster (ds) using the union of the basin geometries.
    
    Parameters:
        ds (xarray.DataArray): The input raster.
        basins (GeoDataFrame): Basin polygons.
        out_path (str, optional): File path to save the clipped image.
        
    Returns:
        Clipped xarray.DataArray.
    """
    # Create a union of all basin polygons
    basins_union = basins.unary_union
    
    # Clip the dataset using the union geometry
    ds_clipped = ds.rio.clip([basins_union], basins.crs, drop=True)

    return ds_clipped

def fill_nosnow(ds,ds_other,no_snow,unknown):
    values=ds.values
    values_other=ds_other.values
    mask=ds.isin(unknown)
    values[mask] = np.where(np.isin(values_other[mask], no_snow), 0, values[mask])
    return ds

def fill_snow(ds,ds_other,rango,unknown):
    values = ds.values
    values_other = ds_other.values

    mask = np.isin(values, unknown) & np.isin(values_other, rango)
    values[mask] = values_other[mask]
    return ds

def filter_snow(ds, no_snow):
    valores_today = ds.values
    valores_today = np.where(np.isin(valores_today, no_snow), 0, valores_today)
    ds.values = valores_today
    return ds

def filter_clouds(ds):
    values=ds.values
    values = np.where(values == 250, np.nan, values)
    ds.values = values
    return ds

def process_all(lista, no_snow, unknown, rango, range_limit=20):
    n = len(lista)
    for i in range(0,n):  # Process ALL elements (0 to n-1)
        ds = lista[i]
        # Apply filter first
        ds = filter_snow(ds, no_snow)
        # Iterate over adjacent entries up to range_limit away
        for j in range(1, range_limit + 1):
            # Process previous entries if available
            if i - j >= 0:
                ds = fill_nosnow(ds, lista[i - j], no_snow, unknown)
                ds = fill_snow(ds, lista[i - j], rango, unknown)
            # Process posterior entries if available
            if i + j < n:
                ds = fill_nosnow(ds, lista[i + j], no_snow, unknown)
                ds = fill_snow(ds, lista[i + j], rango, unknown)
        lista[i] = ds
    return lista

def interpolate_linear(ds):
    """
    Fast interpolation of missing values using SciPy's LinearNDInterpolator.
    
    Parameters:
      ds (rioxarray.DataArray): Input DataArray.
    
    Returns:
      rioxarray.DataArray: Interpolated DataArray.
    """
    # Remove all singleton dimensions with a view (non-copying operation)
    ds_2d = ds.squeeze()
    
    # Access data directly with views instead of copies
    Z = ds_2d.data
    
    # Only proceed with interpolation if there are NaN         
    # Get coordinates - using data instead of values for speed
    x = ds_2d['x'].data
    y = ds_2d['y'].data
    
    # Use efficient masking to identify valid points
    valid_mask = ~np.isnan(Z)
        
    # Optimize point collection for large arrays
    # Pre-allocate arrays instead of using column_stack for large datasets
    # Build meshgrid only once
    X, Y = np.meshgrid(x, y, indexing='xy')
    
    # Pre-allocate points array for valid coordinates
    num_valid = valid_mask.sum()
    pts = np.empty((num_valid, 2), dtype=np.float64)
    pts[:, 0] = X[valid_mask]
    pts[:, 1] = Y[valid_mask]
    
    # Get valid values directly 
    vals = Z[valid_mask]
    
    # Configure interpolator for better performance
    interpolator = LinearNDInterpolator(
        pts, 
        vals,
        fill_value=np.nan,  # Explicitly set fill_value
        rescale=True        # Rescale points for better numerical stability
    )
    
    # Apply interpolation - directly reshape the output for efficiency
    Z_interp = interpolator(X, Y)
    
    # Create output with minimal copying
    ds_2d.data = Z_interp

    return ds_2d

def interpolate(ds):
    """
    Fast interpolation of missing values using SciPy's LinearNDInterpolator.
    
    Parameters:
      ds (rioxarray.DataArray): Input DataArray.
    
    Returns:
      rioxarray.DataArray: Interpolated DataArray.
    """
    # Remove all singleton dimensions with a view (non-copying operation)
    ds_2d = ds.squeeze()
    
    # Access data directly with views instead of copies
    Z = ds_2d.data
    
    # Only proceed with interpolation if there are NaN         
    # Get coordinates - using data instead of values for speed
    x = ds_2d['x'].data
    y = ds_2d['y'].data
    
    # Use efficient masking to identify valid points
    valid_mask = ~np.isnan(Z)
        
    # Optimize point collection for large arrays
    # Pre-allocate arrays instead of using column_stack for large datasets
    # Build meshgrid only once
    X, Y = np.meshgrid(x, y, indexing='xy')
    
    # Pre-allocate points array for valid coordinates
    num_valid = valid_mask.sum()
    pts = np.empty((num_valid, 2), dtype=np.float64)
    pts[:, 0] = X[valid_mask]
    pts[:, 1] = Y[valid_mask]
    
    # Get valid values directly 
    vals = Z[valid_mask]
    
    # Configure interpolator for better performance
    interpolator = NearestNDInterpolator(
        pts, 
        vals
        # NearestNDInterpolator doesn't accept fill_value or rescale parameters
    )
    
    # Apply interpolation - directly reshape the output for efficiency
    Z_interp = interpolator(X, Y)
    
    # Create output with minimal copying
    ds_2d.data = Z_interp

    # fill nans with nearest value

    # Preserve CRS with single call (write_crs returns a new object)
    return ds_2d

def main():
    # Define search parameters
    short_name = "VNP10A1"
    version = "2"
    date_list = pd.date_range('2022-09-01',
    '2022-09-30', freq='D')
    time_dim = len(date_list)
    bounding_box = (-71.71782, -32.28247, -69.809361,
                    -29.0366)  # Global coverage
    no_snow = [211, 237, 239, 251, 252,253]
    unknown = [201,250,254,255]
    rango=list(range(101))
    session = authy()
    basins=load_basins()
    first_ds = None
    for t in date_list:
        temporal = (t, t)  # Query exactly one day
        print(f"Processing date: {t}")

        results = earthaccess.search_data(
            short_name=short_name,
            version=version,
            temporal=temporal,
            bounding_box=bounding_box
        )

        ds_list = process_results(results, session)
        # Create mosaic from the day's granules
        mosaic_utm = mosaic(ds_list)

        # Clip the mosaic using basins
        clipped_ds = clip_image(mosaic_utm, basins)
        
        # Squeeze out the 'band' dimension if present
        if 'band' in clipped_ds.dims:
            clipped_ds = clipped_ds.squeeze('band')
        
        if first_ds is None:
            first_ds = clipped_ds
            y_dim, x_dim = clipped_ds.shape[-2], clipped_ds.shape[-1]
            empty_data = np.full((time_dim, y_dim, x_dim), np.nan, dtype=np.float32)
            combined_ds = xr.DataArray(
                empty_data,
                coords={"time": date_list, "y": clipped_ds["y"], "x": clipped_ds["x"]},
                dims=["time", "y", "x"]
            )
        
        # Reindex clipped_ds to match the coordinates of combined_ds
        clipped_ds = clipped_ds.reindex(y=combined_ds.y, x=combined_ds.x, method='nearest')
        
        # Assign the clipped DataArray to the corresponding time slice
        combined_ds.loc[dict(time=t)] = clipped_ds

    processed_slices = process_all([combined_ds.sel(time=t) for t in date_list], no_snow, unknown, rango, range_limit=20)

    processed_ds = list(map(filter_clouds,processed_slices))
    processed_ds = list(map(interpolate,processed_ds))

    processed_da = xr.concat(processed_ds, dim="time")

if __name__ == "__main__":
    main()
