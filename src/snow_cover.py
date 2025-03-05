from earthaccess import Store
import earthaccess
from io import BytesIO
import rioxarray
import h5py
from rioxarray.merge import merge_arrays
import numpy as np
import geopandas as gpd
import concurrent.futures

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


def process_granule(granule, session, sinusoidal_crs):
    """Process a single granule - used by ThreadPoolExecutor"""
    try:
        # Get URL safely with guard clause
        links = granule.data_links()
        if not links:
            return None
        url = links[0]
        
        # Download data efficiently
        response = session.get(url, stream=True)
        response.raise_for_status()
        content = response.content
        
        # Create BytesIO buffer once
        data_bytes = BytesIO(content)
        
        # Open dataset efficiently
        ds = rioxarray.open_rasterio(
            data_bytes, 
            engine='h5netcdf'
        ).astype('float32')  # Convert to float32 early
        
        # Extract NDSI data
        ds = ds.HDFEOS_GRIDS_VIIRS_Grid_IMG_2D_Data_Fields_NDSI
        
        # Reset pointer for coordinate extraction
        data_bytes.seek(0)
        x_dim, y_dim = get_coords(data_bytes)
        
        # Assign coordinates and CRS in one operation
        ds = ds.assign_coords(x=x_dim, y=y_dim)
        ds = ds.rio.write_crs(sinusoidal_crs, inplace=False)
        
        return ds
        
    except Exception as e:
        print(f"Error processing {url if 'url' in locals() else 'unknown'}: {e}")
        return None

def process_results(results, session):
    """Process multiple granules in parallel"""
    sinusoidal_crs = "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
    
    # Use ThreadPoolExecutor for concurrent downloads (I/O bound)
    with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
        # Submit all tasks
        future_to_granule = {
            executor.submit(process_granule, granule, session, sinusoidal_crs): granule 
            for granule in results
        }
        
        # Collect results as they complete
        ds_list = []
        for future in concurrent.futures.as_completed(future_to_granule):
            ds = future.result()
            if ds is not None:
                ds_list.append(ds)
    
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

def main():
    # Define search parameters
    short_name = "VNP10A1"
    version = "2"
    date_list = ["2022-09-29", "2022-09-30", "2022-10-01"]
    bounding_box = (-71.71782, -32.28247, -69.809361,
                    -29.0366)  # Global coverage
    session = authy()
    basins=load_basins()

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

        # Optionally save the day's mosaic before clipping (commented out)
        # mosaic_utm.to_netcdf(f"mosaic_{t}.nc")

        # Clip the mosaic using basins
        clipped_ds = clip_image(mosaic_utm, basins)

if __name__ == "__main__":
    main()
