from earthaccess import Auth, Store
import earthaccess
import xarray as xr
import requests
from io import BytesIO
import rioxarray
import h5py

def get_coords(data_bytes):
    with h5py.File(data_bytes, 'r') as f:
        x_coords = f['HDFEOS/GRIDS/VIIRS_Grid_IMG_2D/XDim'][:]
        y_coords = f['HDFEOS/GRIDS/VIIRS_Grid_IMG_2D/YDim'][:]
    return x_coords, y_coords

# Initialize authentication properly
auth = earthaccess.login()
store = Store(auth=auth)
session = auth.get_session()

# Initialize S3 filesystem
short_name = "VNP10A1"
version = "2"
temporal=("2022-09-29","2022-09-30")
bounding_box = (-71.71782, -32.28247, -69.809361,
                 -29.0366)  # Global coverage

results = earthaccess.search_data(
    short_name=short_name,
    version=version,
    temporal=temporal,
    bounding_box=bounding_box
)

for granule in results:
  try:
      # Get HTTPS URL
      url = granule.data_links()[0]
      
      # Download data using authenticated session
      response = session.get(url, stream=True)
      response.raise_for_status()
      
      # Read into memory and inspect with h5py
      data_bytes = BytesIO(response.content)
      ds = rioxarray.open_rasterio(data_bytes, engine='h5netcdf')
      x_dim, y_dim = get_coords(data_bytes)
      ds = ds.assign_coords(x=x_dim, y=y_dim)
    #   ds.rio.set_spatial_dims(x_dim='x', y_dim='y', inplace=True)
      ds_snow=ds.HDFEOS_GRIDS_VIIRS_Grid_IMG_2D_Data_Fields_NDSI
 
      sinusoidal_crs = "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
      ds_snow = ds_snow.rio.write_crs(sinusoidal_crs)
      
      # Reproject to UTM zone 19S
      ds_utm = ds_snow.rio.reproject("EPSG:32719")
      ds_utm.plot()


  except Exception as e:
      print(f"Error processing {url}: {e}")
    
