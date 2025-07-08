import geopandas as gpd
import rasterio
import numpy as np
import pandas as pd
import xarray as xr
from numba import jit, prange
import geopandas as gpd
import xarray as xr
import numpy as np
from rasterio.features import rasterize
from rasterio.transform import from_bounds
_RASTER_CACHE = {
'bounds': None,      # (x0, x1, y0, y1, width, height)
'transform': None,   # affine.Transform
'coords': None,      # {'y': array, 'x': array}
'geoms': None        # cached list of geometries
}
# === OPTIMIZED JIT FUNCTIONS ===
@jit(nopython=True, parallel=True, cache=True)
def vectorized_idw_3d(stream_coords_3d, glofas_coords_3d, glofas_q_values, 
                      max_dist_2d=15000, p=2.0):
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
        crs="EPSG:32719")
    return glofas_gdf

# === 2. Extract elevation at GloFAS points ===
def extract_elevations(glofas_gdf, ELEV_RASTER):
    # Vectorized coordinate extraction
    glofas_coords = np.column_stack([glofas_gdf.geometry.x, glofas_gdf.geometry.y])
    
    # Batch elevation sampling
    with rasterio.open(ELEV_RASTER) as src:
        # Use rasterio's vectorized sampling
        elevations = np.array([val[0] for val in src.sample(glofas_coords)])
    
    glofas_gdf["elevation"] = elevations
    return glofas_gdf

# === 3. Load stream segments ===
def load_stream_data(stream_file, target_crs="EPSG:32719"):
    streams = gpd.read_file(stream_file).to_crs(target_crs)
    if "segment_id" not in streams.columns:
        streams["segment_id"] = streams.index.astype(str)
    
    # Vectorized centroid calculation
    centroids = streams.geometry.centroid
    stream_coords = np.column_stack([centroids.x, centroids.y])
    
    return streams, stream_coords

# === 4. Extract elevation at stream points ===
def extract_stream_elevations(stream_coords, ELEV_RASTER):
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
                stream_elev, stream_coords, max_dist=15000, p=2.0):
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
def assign_results(streams,col="_Q_",
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

def gdf2raster_cached(gdf, template_2d, col='Q_assigned_3d_idw', res=250, fill=0):
    """Ultra-fast cached rasterize using template_2d for bounds/transform."""
    global _RASTER_CACHE

    # 1) initialize bounds/transform/coords from template_2d
    if _RASTER_CACHE['bounds'] is None:
        x0, x1 = template_2d.x.min().item(), template_2d.x.max().item()
        y0, y1 = template_2d.y.min().item(), template_2d.y.max().item()
        w = int(np.ceil((x1 - x0) / res))
        h = int(np.ceil((y1 - y0) / res))
        tf = from_bounds(x0, y0, x1, y1, w, h)
        ys = np.linspace(y1, y0, h, dtype='float64')
        xs = np.linspace(x0, x1, w, dtype='float64')
        _RASTER_CACHE.update({
            'bounds': (x0, x1, y0, y1, w, h),
            'transform': tf,
            'coords': {'y': ys, 'x': xs}
        })

    x0, x1, y0, y1, w, h = _RASTER_CACHE['bounds']
    tf = _RASTER_CACHE['transform']
    coords = _RASTER_CACHE['coords']

    # 2) cache geometries once
    if _RASTER_CACHE['geoms'] is None:
        _RASTER_CACHE['geoms'] = list(gdf.geometry)

    geoms = _RASTER_CACHE['geoms']
    vals = gdf[col].fillna(fill).to_numpy()

    # 3) generator of (geom, value) pairs
    shapes = ((geom, v) for geom, v in zip(geoms, vals))

    # 4) rasterize
    raster = rasterize(
        shapes=shapes,
        out_shape=(h, w),
        transform=tf,
        fill=fill,
        all_touched=True,
        dtype='float32'
    )

    # 5) return DataArray (defer .rio.write_crs to after concat)
    return xr.DataArray(raster, coords=coords, dims=['y', 'x'])

def interp_glofas(glofas_nc="nc_test.nc"):
    """
    Main function to interpolate GloFAS discharge data onto stream segments.
    This function handles the entire workflow from loading data to saving results.
    """
    # === 0. Setup ===
    print("Starting GloFAS interpolation...")
    # === PARAMETERS ===    
    stream_file = "../geodata/riverQ.gpkg"
    target_crs = "EPSG:32719"
    col='_Q_'
    ELEV_RASTER= "../Rst/dem90fill.tif"

    # === 0. Setup ===
    # === 1. Load GloFAS discharge and convert to grid polygons ===
    ds_clipped = xr.open_dataset(glofas_nc)
    geopath='RegionCoquimbo.geojson'
    shapefile = gpd.read_file(geopath)
    ds_clipped["longitude"] = ds_clipped["longitude"].where(ds_clipped["longitude"] <= 180, ds_clipped["longitude"] - 360)
    ds_clipped = ds_clipped.rio.write_crs("EPSG:4326")
    shapefile = shapefile.to_crs("EPSG:4326")

    # Alternative: Keep all ensemble members in a 5D array
    template_2d = ds_clipped['dis24'].isel(forecast_period=0, forecast_reference_time=0, number=0).rio.reproject(target_crs)

    # Usage:
    # glofas_gdf = build_glofas_gdf(dis, lats, lons, target_crs)
    lat_name = 'y'
    lon_name = 'x'
    lats = template_2d[lat_name].values
    lons = template_2d[lon_name].values

    #### reproject all scenarios

    shape = (
        len(ds_clipped['forecast_period']), 
        len(ds_clipped['forecast_reference_time']), 
        len(ds_clipped['number']),
        len(template_2d.y), 
        len(template_2d.x)
    )

    ds_utm = xr.DataArray(
        np.full(shape, np.nan, dtype=ds_clipped['dis24'].dtype),
        dims=['forecast_period', 'forecast_reference_time', 'number', 'y', 'x'],
        coords={
            'forecast_period': ds_clipped['forecast_period'], 
            'forecast_reference_time': ds_clipped['forecast_reference_time'],
            'number': ds_clipped['number'],
            'y': template_2d.y, 
            'x': template_2d.x
        },
        name='dis24'
    ).rio.write_crs(target_crs)

    # Reproject all slices
    for i in range(len(ds_clipped['forecast_period'])):
        for j in range(len(ds_clipped['forecast_reference_time'])):
            for k in range(len(ds_clipped['number'])):
                slice_2d = ds_clipped['dis24'].isel(forecast_period=i, forecast_reference_time=j, number=k)
                reprojected_2d = slice_2d.rio.reproject_match(template_2d)
                ds_utm[i, j, k, :, :] = reprojected_2d.values

    print("✅ Reprojection with all ensemble members complete!")

    # load streams and retrieve coordinates
    streams, stream_coords = load_stream_data(stream_file, target_crs)

    # sample stream elevations at stream coordinates
    stream_elev = extract_stream_elevations(stream_coords,ELEV_RASTER)

    ras_list = []

    for i, time_ in enumerate(ds_utm['forecast_period']):
        ras_time_list = []
        for j, ensemble in enumerate(ds_utm['number']):
            print(f"Processing time {i+1}/{len(ds_utm['forecast_period'])}, ensemble {j+1}/{len(ds_utm['number'])}...")

            dis_utm = ds_utm.sel(forecast_period=time_, number=ensemble)

            # === 5. Create 3D coordinates for BallTree ===
            # Scale elevation to match horizontal distance units
            glofas_gdf = build_glofas_gdf(dis_utm, lats, lons, target_crs)
            glofas_gdf = extract_elevations(glofas_gdf, ELEV_RASTER)
            
            glofas_coords_3d, stream_coords_3d = prepare_3d_coordinates(glofas_gdf, stream_coords, stream_elev)

            idw_values_3d = compute_idw(stream_coords_3d, glofas_coords_3d, glofas_gdf, stream_elev, stream_coords)

            streams_clean = assign_results(streams, col, idw_values_3d)

            # Create raster
            ras = gdf2raster_cached(streams_clean, template_2d, col=col, res=250, fill=np.nan)
            ras_time_list.append(ras)
        
        # Combine ensemble members for this time step
        ras_time_combined = xr.concat(ras_time_list, dim='number')
        ras_time_combined = ras_time_combined.assign_coords(number=ds_utm['number'].values)
        ras_list.append(ras_time_combined)

    # === 9. Combine all time steps ===
    ras_combined = xr.concat(ras_list, dim='forecast_period')
    ras_combined = ras_combined.assign_coords(forecast_period=ds_utm['forecast_period'].values)

    # Set proper CRS and save
    ras_combined.rio.write_crs(target_crs, inplace=True)

    return ras_combined

def main():
    glofas_nc = "nc_test.nc"  # Path to your GloFAS NetCDF file
    result = interp_glofas(glofas_nc)
    
    # Save the result to a NetCDF file
    output_path = "glofas_interpolated.nc"
    result.to_netcdf(output_path)
    print(f"Results saved to {output_path}")

if __name__ == "__main__":
    main()
    # Run the main function