import geopandas as gpd
import rasterio
import numpy as np
import pandas as pd
import xarray as xr
from numba import jit, prange
import time
import geopandas as gpd
import xarray as xr
import numpy as np
from rasterio.features import rasterize
from rasterio.transform import from_bounds

# === OPTIMIZED JIT FUNCTIONS ===
@jit(nopython=True, parallel=True, cache=True)
def vectorized_idw_3d(stream_coords_3d, glofas_coords_3d, glofas_q_values, 
                      max_dist_2d, p=2.0):
    """Ultra-fast vectorized 3D IDW computation"""
    n_streams = stream_coords_3d.shape[0]
    n_glofas = glofas_coords_3d.shape[0]
    results = np.zeros(n_streams)
    
    for i in prange(n_streams):
        # Compute all 3D distances at once
        dx = glofas_coords_3d[:, 0] - stream_coords_3d[i, 0]
        dy = glofas_coords_3d[:, 1] - stream_coords_3d[i, 1]
        dz = glofas_coords_3d[:, 2] - stream_coords_3d[i, 2]
        
        # 2D and 3D distances
        d_2d = np.sqrt(dx*dx + dy*dy)
        d_3d = np.sqrt(dx*dx + dy*dy + dz*dz)
        
        # Filter by 2D distance
        valid_mask = d_2d < max_dist_2d
        
        if np.sum(valid_mask) == 0:
            results[i] = 0.0
            continue
        
        # Use 3D distances for weighting
        d_3d_valid = d_3d[valid_mask]
        q_valid = glofas_q_values[valid_mask]
        
        # Avoid division by zero
        d_3d_valid = np.where(d_3d_valid < 1e-10, 1e-10, d_3d_valid)
        
        # IDW weights
        weights = 1.0 / (d_3d_valid ** p)
        weight_sum = np.sum(weights)
        
        if weight_sum > 0:
            weights = weights / weight_sum
            results[i] = np.sum(weights * q_valid)
        else:
            results[i] = np.mean(q_valid)
    
    return results

@jit(nopython=True, cache=True)
def nash_sutcliffe_jit(obs, sim):
    """JIT-compiled Nash-Sutcliffe efficiency"""
    obs_mean = np.mean(obs)
    num = np.sum((obs - sim) ** 2)
    den = np.sum((obs - obs_mean) ** 2)
    return 1.0 - num/den if den != 0 else np.nan

    # Build GloFAS points

def build_glofas_gdf(dis, lats, lons, target_crs):
    """
    Vectorized extraction of valid GloFAS points and discharge values.
    """
    # Get the discharge array as 2D (lat, lon)
    q_arr = dis.values
    # If q_arr has more than 2 dimensions, squeeze or slice as needed
    while q_arr.ndim > 2:
        q_arr = q_arr[0]
    # Build 2D coordinate grids
    lon2d, lat2d = np.meshgrid(lons, lats)
    # Mask valid (non-NaN) discharge values
    valid = ~np.isnan(q_arr)
    xs = lon2d[valid]
    ys = lat2d[valid]
    qs = q_arr[valid]
    # Build GeoDataFrame in one shot
    glofas_gdf = gpd.GeoDataFrame(
        {'Q': qs},
        geometry=gpd.points_from_xy(xs, ys),
        crs="EPSG:4326"
    ).to_crs(target_crs)
    return glofas_gdf

# === 2. Extract elevation at GloFAS points ===
def extract_elevations(glofas_gdf):
    # Vectorized coordinate extraction
    glofas_coords = np.column_stack([glofas_gdf.geometry.x, glofas_gdf.geometry.y])
    
    # Batch elevation sampling
    with rasterio.open(ELEV_RASTER) as src:
        # Use rasterio's vectorized sampling
        elevations = np.array([val[0] for val in src.sample(glofas_coords)])
    
    glofas_gdf["elevation"] = elevations
    return glofas_gdf

# === 3. Load stream segments ===
def load_stream_data(stream_file):
    streams = gpd.read_file(stream_file).to_crs(target_crs)
    if "segment_id" not in streams.columns:
        streams["segment_id"] = streams.index.astype(str)
    
    # Vectorized centroid calculation
    centroids = streams.geometry.centroid
    stream_coords = np.column_stack([centroids.x, centroids.y])
    
    return streams, stream_coords

# === 4. Extract elevation at stream points ===
def extract_stream_elevations():
    with rasterio.open(ELEV_RASTER) as src:
        stream_elev = np.array([val[0] for val in src.sample(stream_coords)])
    return stream_elev

def prepare_3d_coordinates(glofas_gdf, stream_coords, stream_elev, 
                           elev_scale=0.1):
    # Vectorized 3D coordinate creation
    glofas_coords_3d = np.column_stack([
        glofas_gdf.geometry.x.values,
        glofas_gdf.geometry.y.values,
        glofas_gdf["elevation"].values * elev_scale
    ])
    
    stream_coords_3d = np.column_stack([
        stream_coords[:, 0],
        stream_coords[:, 1],
        stream_elev * elev_scale
    ])
    
    return glofas_coords_3d, stream_coords_3d

# === 6. IDW with 3D distance ===
def compute_idw(stream_coords_3d, glofas_coords_3d, glofas_gdf, 
                stream_elev, max_dist=15000, p=2.0):
    print(f"   Processing {len(stream_coords_3d)} stream points with {len(glofas_coords_3d)} GloFAS points")
    
    # Handle NaN elevations efficiently
    valid_stream_mask = ~np.isnan(stream_elev)
    
    # Initialize results
    idw_values = np.zeros(len(stream_coords_3d))
    
    if np.any(valid_stream_mask):
        # Process valid elevations with 3D distance
        valid_indices = np.where(valid_stream_mask)[0]
        valid_stream_coords_3d = stream_coords_3d[valid_indices]
        
        # Vectorized IDW computation
        valid_results = vectorized_idw_3d(
            valid_stream_coords_3d,
            glofas_coords_3d,
            glofas_gdf["Q"].values,
            max_dist,
            p
        )
        
        idw_values[valid_indices] = valid_results
    
    # Handle NaN elevations with 2D fallback (if any)
    invalid_mask = np.isnan(stream_elev)
    if np.any(invalid_mask):
        print(f"   Fallback 2D processing for {np.sum(invalid_mask)} points with NaN elevation")
        
        invalid_indices = np.where(invalid_mask)[0]
        glofas_coords_2d = glofas_coords_3d[:, :2]  # Just x, y
        
        for idx in invalid_indices:
            # Simple 2D IDW for fallback
            target_2d = stream_coords[idx]
            distances_2d = np.sqrt(np.sum((glofas_coords_2d - target_2d)**2, axis=1))
            
            valid_neighbors = distances_2d < max_dist
            if np.any(valid_neighbors):
                d_valid = distances_2d[valid_neighbors]
                q_valid = glofas_gdf["Q"].values[valid_neighbors]
                
                d_valid = np.where(d_valid < 1e-10, 1e-10, d_valid)
                weights = 1.0 / (d_valid ** p)
                weights /= weights.sum()
                idw_values[idx] = np.sum(weights * q_valid)
    
    return idw_values

# === 7. Assign interpolated discharge to stream segments ===
def assign_results(streams,col="Q_assigned_3d_idw",
                   idw_values_3d=None):
    streams[col] = idw_values_3d
    
    # Clean and prepare for saving
    streams_clean = streams.copy()
    if 'centroid' in streams_clean.columns:
        streams_clean = streams_clean.drop(columns=['centroid'])
    
    streams_clean[col] = pd.to_numeric(
        streams_clean[col], errors='coerce'
    )
    
    return streams_clean

# === 8. Save output ===
# Clean the dataframe before saving
def save_results():
    streams_clean.to_file("streams_3d_idw_optimized.gpkg",
                        driver="GPKG")

def gdf2raster(gdf,
                col='Q_assigned_3d_idw',
                res=500, fill=np.nan):
    """    Convert a GeoDataFrame to a raster using a specified column for values.
    Args:
        gdf (GeoDataFrame): Input GeoDataFrame with geometries and values.
        col (str): Column name in GeoDataFrame to use for raster values.
        res (int): Resolution of the output raster in pixels.
        fill (float): Value to fill empty pixels in the raster.
    Returns:
        xr.DataArray: Rasterized data as an xarray DataArray.
    """
    # Ensure the GeoDataFrame has the specified column
    # Ultra-fast optimized
    res, col = 500, col
    x_min, x_max, y_min, y_max = dis_utm.x.min(), dis_utm.x.max(), dis_utm.y.min(), dis_utm.y.max()
    width, height = int(np.ceil((x_max - x_min) / res)), int(np.ceil((y_max - y_min) / res))

    # Fixed rasterize call
    valid = gdf.dropna(subset=[col])
    raster = rasterize(
        shapes=[(g, v) for g, v in zip(valid.geometry, valid[col])], 
        out_shape=(height, width), 
        transform=from_bounds(x_min, y_min, x_max, y_max, width, height),
        fill=0, 
        all_touched=True, 
        dtype='float32'
    )

    # Save
    da=xr.DataArray(raster, coords={'y': np.linspace(y_max, y_min, height), 
                                'x': np.linspace(x_min, x_max, width)}, 
                dims=['y', 'x']).rio.write_crs(gdf.crs)

    print(f"✅ {np.sum(~np.isnan(raster))} pixels")
    return da

# === PARAMETERS ===
glofas_nc = "../Rst/GloFAS_2025_06_13_f.nc"
glofas_nc='../Rst/dis_1980_2018_clip.nc'
var_name = "dis"
time_idx = 0
ELEV_RASTER= "../Rst/dem90fill.tif"
stream_file = "../geodata/riverQ.gpkg"
target_crs = "EPSG:32719"
col='Q_assigned_3d_idw'
p = 2  # IDW power parameter
elev_scale = 0.1  # Scale factor for elevation to match horizontal distance units
max_dist = 15000  # Max distance in meters for considering neighbors

ras_list = []
# === 0. Setup ===
# === 1. Load GloFAS discharge and convert to grid polygons ===
ds = xr.open_dataset(glofas_nc)
times = ds['time'].values
streams, stream_coords = load_stream_data(stream_file)

for time in times:
    dis = ds[var_name].sel(time=time)
    lat_name = [dim for dim in dis.dims if 'lat' in str(dim)][0]
    lon_name = [dim for dim in dis.dims if 'lon' in str(dim)][0]
    lats = dis[lat_name].values
    lons = dis[lon_name].values

    dis.rio.write_crs("EPSG:4326", inplace=True)
    dis_utm=dis.rio.reproject(target_crs)

    # Usage:
    # glofas_gdf = build_glofas_gdf(dis, lats, lons, target_crs)
    glofas_gdf = build_glofas_gdf(dis, lats, lons, target_crs)

    glofas_gdf = extract_elevations(glofas_gdf)

    stream_elev = extract_stream_elevations()
    # === 5. Create 3D coordinates for BallTree ===
    # Scale elevation to match horizontal distance units

    glofas_coords_3d, stream_coords_3d = prepare_3d_coordinates(glofas_gdf, 
                                                                stream_coords, 
                                                                stream_elev)

    idw_values_3d = compute_idw(stream_coords_3d, glofas_coords_3d, glofas_gdf, 
                stream_elev)

    streams_clean = assign_results(streams,col,
                   idw_values_3d)

    ###### tests

    ras = gdf2raster(streams_clean,
                col=col,
                res=500, fill=np.nan)

    ras_list.append(ras)

# === 9. Save results ===
ras_combined = xr.concat(ras_list, dim='time')
ras_combined = ras_combined.assign_coords(
    time=pd.to_datetime(times)).rename({'time': 'time'})
ras_combined.rio.write_crs(target_crs, inplace=True)
ras_combined.to_netcdf("dis_3d_idw_optimized_1980_2018.nc")