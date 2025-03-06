from earthaccess import Store
import earthaccess
from io import BytesIO
import rioxarray
import h5py
from rioxarray.merge import merge_arrays
import numpy as np
import geopandas as gpd
import concurrent.futures
import os

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

def main():
    # Define search parameters
    short_name = "VNP10A1"
    version = "2"
    date_list = ["2022-09-29", "2022-09-30", "2022-10-01"]
    bounding_box = (-71.71782, -32.28247, -69.809361,
                    -29.0366)  # Global coverage
    no_snow = [211, 237, 239, 251, 252,253]
    unknown = [201,250,254,255]
    rango=list(range(101))
    session = authy()
    basins=load_basins()

    lista=[None,None,None]
    for ind,t in enumerate(date_list):
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
        # get clipped_ds resolution
        # get unique values of dataset


        # , 254, 255]

        # Use the .where method to keep original values where condition is False, and set to 0 otherwise.

        lista[ind]=clipped_ds.copy()

        lista[-2]=fill_nosnow(lista[-2],lista[-3],no_snow,unknown)
        lista[-2]=fill_nosnow(lista[-2],lista[-1],no_snow,unknown)
        lista[-2]=fill_snow(lista[-2],lista[-3],rango,unknown)
        lista[-2]=fill_snow(lista[-2],lista[-1],rango,unknown)
        lista[-2]=filter_snow(lista[-2],no_snow)
        lista[-2].rio.to_raster("today_correct_new.tif")






if __name__ == "__main__":
    main()
