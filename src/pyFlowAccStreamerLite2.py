import geopandas as gpd
import rasterio
import numpy as np
import xarray as xr
from shapely.geometry import Point
from sklearn.neighbors import BallTree
from collections import defaultdict
import numba
from numba import njit
import pandas as pd
from scipy.spatial import cKDTree
from joblib import Parallel, delayed

# === PARAMETERS ===
glofas_nc = "../Rst/GloFAS_2025_06_13_f.nc"
var_name = "dis24"
time_idx = 0
accum_raster = "../Rst/flow_accumulation.tif"
flow_dir_raster = "../Rst/flowDir.tif"
stream_file = "../geodata/riverQ.gpkg"
target_crs = "EPSG:32719"
p = 2
alpha = 0.5
max_dist = 15000

# === 1. OPTIMIZED: Load GloFAS data with chunking and vectorization ===
print("Loading GloFAS data...")
with xr.open_dataset(glofas_nc, chunks={'latitude': 100, 'longitude': 100}) as ds:
    dis = ds[var_name][time_idx, :, :].isel(forecast_period=0)
    
    # Handle dimensionality
    if dis.ndim == 3:
        dis2d = dis.isel(forecast_reference_time=0)
    else:
        dis2d = dis
    
    # Vectorized coordinate creation
    lon_vals = dis2d.longitude.values
    lat_vals = dis2d.latitude.values
    
    # Create mask first to reduce memory
    dis_vals = dis2d.values
    mask = (dis_vals > 0) & np.isfinite(dis_vals)
    
    if not np.any(mask):
        raise ValueError("No valid discharge values found")
    
    # Only create coordinates for valid points
    lat_indices, lon_indices = np.where(mask)
    coords = np.column_stack((lon_vals[lon_indices], lat_vals[lat_indices]))
    Q_vals = dis_vals[mask]

print(f"Found {len(Q_vals)} valid GloFAS points")

# Create GeoDataFrame more efficiently
glofas_gdf = gpd.GeoDataFrame(
    {'Q': Q_vals},
    geometry=gpd.points_from_xy(coords[:, 0], coords[:, 1]),
    crs="EPSG:4326"
).to_crs(target_crs)

# === 2. OPTIMIZED: Vectorized raster sampling ===
print("Extracting flow accumulation...")
with rasterio.open(accum_raster) as src:
    # Vectorized sampling
    coord_pairs = list(zip(glofas_gdf.geometry.x, glofas_gdf.geometry.y))
    fa_values = np.array([val[0] for val in src.sample(coord_pairs)])
    glofas_gdf["flow_accum"] = np.clip(fa_values, 1, None)

# === 3. Load streams (minimal processing) ===
print("Loading streams...")
streams = gpd.read_file(stream_file).to_crs(target_crs)
if "segment_id" not in streams.columns:
    streams["segment_id"] = streams.index.astype(str)

# Pre-compute centroids
stream_centroids = streams.geometry.centroid
stream_coords = np.column_stack([stream_centroids.x, stream_centroids.y])

# === 4. OPTIMIZED: Flow direction processing with Numba ===
print("Processing flow directions...")
with rasterio.open(flow_dir_raster) as src:
    flow_dir = src.read(1)
    transform = src.transform
    height, width = flow_dir.shape

@njit
def build_reverse_graph_fast(flow_dir):
    """Numba-optimized reverse graph building"""
    height, width = flow_dir.shape
    
    # Use arrays instead of dicts for Numba compatibility
    # Store as flattened indices
    max_connections = height * width * 8  # Conservative estimate
    upstream_cells = np.full(max_connections, -1, dtype=np.int32)
    downstream_cells = np.full(max_connections, -1, dtype=np.int32)
    
    # D8 directions
    dr = np.array([-1, -1, 0, 1, 1, 1, 0, -1], dtype=np.int32)
    dc = np.array([0, 1, 1, 1, 0, -1, -1, -1], dtype=np.int32)
    flags = np.array([64, 128, 1, 2, 4, 8, 16, 32], dtype=np.int32)
    
    conn_count = 0
    
    for i in range(height):
        for j in range(width):
            current_dir = flow_dir[i, j]
            
            for k in range(8):
                if current_dir == flags[k]:
                    ni = i + dr[k]
                    nj = j + dc[k]
                    
                    if 0 <= ni < height and 0 <= nj < width:
                        upstream_idx = i * width + j
                        downstream_idx = ni * width + nj
                        
                        upstream_cells[conn_count] = upstream_idx
                        downstream_cells[conn_count] = downstream_idx
                        conn_count += 1
                    break
    
    return upstream_cells[:conn_count], downstream_cells[:conn_count]

upstream_arr, downstream_arr = build_reverse_graph_fast(flow_dir)

# Convert to efficient lookup structure
print("Building lookup structures...")
reverse_graph = defaultdict(list)
for i in range(len(upstream_arr)):
    reverse_graph[downstream_arr[i]].append(upstream_arr[i])

# === 5. OPTIMIZED: Coordinate mapping with vectorization ===
def coords_to_indices_vectorized(x_coords, y_coords, transform, width, height):
    """Vectorized coordinate to index conversion"""
    inv_transform = ~transform
    cols, rows = inv_transform * (x_coords, y_coords)
    
    # Vectorized bounds checking
    valid_mask = (rows >= 0) & (rows < height) & (cols >= 0) & (cols < width)
    indices = np.full(len(x_coords), -1, dtype=np.int32)
    indices[valid_mask] = (rows[valid_mask].astype(np.int32) * width + 
                          cols[valid_mask].astype(np.int32))
    
    return indices

# Map all coordinates at once
glofas_indices = coords_to_indices_vectorized(
    glofas_gdf.geometry.x.values, 
    glofas_gdf.geometry.y.values,
    transform, width, height
)
glofas_gdf["cell_index"] = glofas_indices

stream_indices = coords_to_indices_vectorized(
    stream_coords[:, 0], stream_coords[:, 1],
    transform, width, height
)

# Filter valid points
valid_glofas = glofas_gdf[glofas_gdf["cell_index"] >= 0].copy()
print(f"Valid GloFAS points: {len(valid_glofas)}")

# === 6. OPTIMIZED: Upstream finding with NumPy ===
@njit
def find_upstream_numba(start_cell, upstream_arr, downstream_arr):
    """Numba-optimized upstream cell finding"""
    visited = set()
    stack = [start_cell]
    
    while stack:
        cell = stack.pop()
        if cell in visited:
            continue
        visited.add(cell)
        
        # Find upstream cells
        for i in range(len(downstream_arr)):
            if downstream_arr[i] == cell:
                stack.append(upstream_arr[i])
    
    return visited

# === 7. OPTIMIZED: Parallel IDW computation ===
def compute_idw_for_stream(i, stream_idx, valid_glofas, upstream_arr, downstream_arr):
    """Compute IDW for a single stream segment"""
    if stream_idx < 0:
        return 0.0
    
    # Find upstream cells
    upstream_cells = find_upstream_numba(stream_idx, upstream_arr, downstream_arr)
    
    if not upstream_cells:
        return 0.0
    
    # Find relevant GloFAS points
    relevant_mask = valid_glofas["cell_index"].isin(upstream_cells)
    relevant_points = valid_glofas[relevant_mask]
    
    if len(relevant_points) == 0:
        return 0.0
    
    # Compute distances
    stream_coord = stream_coords[i:i+1]
    glofas_coords = np.column_stack([relevant_points.geometry.x, relevant_points.geometry.y])
    
    dists = np.linalg.norm(glofas_coords - stream_coord, axis=1)
    valid_dist_mask = dists <= max_dist
    
    if not np.any(valid_dist_mask):
        return 0.0
    
    # Apply filters
    dists = dists[valid_dist_mask]
    q_vals = relevant_points["Q"].values[valid_dist_mask]
    fa_vals = relevant_points["flow_accum"].values[valid_dist_mask]
    
    # IDW calculation
    weights = (fa_vals ** alpha) / (dists ** p)
    total_weight = np.sum(weights)
    
    if total_weight > 0:
        return np.sum(weights * q_vals) / total_weight
    return 0.0

# Parallel computation
print("Computing IDW interpolation...")
n_jobs = min(4, len(streams))  # Limit to 4 cores

idw_values = Parallel(n_jobs=n_jobs, verbose=1)(
    delayed(compute_idw_for_stream)(i, stream_indices[i], valid_glofas, upstream_arr, downstream_arr)
    for i in range(len(streams))
)

# === 8. Assign results ===
streams["Q_assigned_mfd_idw"] = idw_values

# === 9. Save output ===
print("Saving results...")
streams.to_file("streams_mfd_weighted.gpkg")
print("Complete!")