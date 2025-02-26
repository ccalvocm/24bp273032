from earthaccess import Store
import earthaccess
from io import BytesIO
import rioxarray
import h5py
from rioxarray.merge import merge_arrays
import numpy as np

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
    sinusoidal_crs = "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
    ds_list = []
    for granule in results:
        try:
            # Get HTTPS URL
            url = granule.data_links()[0]
            
            # Download data using authenticated session
            response = session.get(url, stream=True)
            #   response.raise_for_status()
            
            # Read into memory and inspect with h5py
            data_bytes = BytesIO(response.content)
            data_bytes.seek(0)
            ds = rioxarray.open_rasterio(data_bytes, engine='h5netcdf')
            ds=ds.HDFEOS_GRIDS_VIIRS_Grid_IMG_2D_Data_Fields_NDSI
            x_dim, y_dim = get_coords(data_bytes)
            ds = ds.assign_coords(x=x_dim, y=y_dim)
            ds = ds.rio.write_crs(sinusoidal_crs, inplace=True)
            ds_list.append(ds)

        except Exception as e:
            print(f"Error processing {url}: {e}")
        
    return ds_list
        
def mosaic(ds_list):
        
def mosaic(ds_list):
    mosaic_da = merge_arrays(ds_list)
    mosaic_utm = mosaic_da.rio.reproject("EPSG:32719")
    return mosaic_utm

def main():
    # Define search parameters
    short_name = "VNP10A1"
    version = "2"
    temporal=("2022-09-29","2022-09-30")
    bounding_box = (-71.71782, -32.28247, -69.809361,
                    -29.0366)  # Global coverage

    results = earthaccess.search_data(
        short_name=short_name,
        version=version,
        temporal=temporal,
        bounding_box=bounding_box)
    
    session = authy()
    ds_list = process_results(results, session)
    mosaic_utm = mosaic(ds_list)
    plot(mosaic_utm)

if __name__ == "__main__":
    main()