import geopandas as gpd
import rasterio
import numpy as np
import xarray as xr
from shapely.geometry import Point
from sklearn.neighbors import BallTree
from collections import defaultdict
import numba
from numba import njit

# === PARAMETERS ===
glofas_nc = "../Rst/GloFAS_2025_06_13_f.nc"
var_name = "dis24"
time_idx = 0
accum_raster = "../Rst/flow_accumulation.tif"
flow_dir_raster = "../Rst/flowDir.tif"  # MFD flow direction raster
stream_file = "../geodata/riverQ.gpkg"
target_crs = "EPSG:32719"
p = 2  # IDW power parameter
alpha = 0.5  # Flow accumulation weight exponent
max_dist = 15000  # Max distance in meters for considering neighbors

# === 1. Load GloFAS discharge and convert to points ===
ds = xr.open_dataset(glofas_nc)
dis = ds[var_name][time_idx, :, :].isel(forecast_period=0)

# Create mask for valid values
# If dis is 3D (e.g., time, lat, lon), select a time slice or reduce to 2D
if dis.ndim == 3:
    dis2d = dis.isel(forecast_reference_time=0)  # or dis.mean(dim='time')
else:
    dis2d = dis

# Create 2D coordinate grids
lon2d, lat2d = np.meshgrid(dis2d.longitude.values, dis2d.latitude.values)

# Create mask for valid values (e.g., nonzero and finite)
mask = (dis2d.values > 0) & np.isfinite(dis2d.values)  # mask shape matches dis2d.values

# Extract coordinates and values using the mask
coords = np.column_stack((lon2d[mask], lat2d[mask]))
Q_vals = dis2d.values[mask]

glofas_gdf = gpd.GeoDataFrame(
    {'Q': Q_vals},
    geometry=gpd.points_from_xy(coords[:, 0], coords[:, 1]),
    crs="EPSG:4326"
).to_crs(target_crs)

# === 2. Extract flow accumulation at GloFAS points ===
with rasterio.open(accum_raster) as src:
    glofas_gdf["flow_accum"] = [val[0] for val in src.sample(zip(glofas_gdf.geometry.x, glofas_gdf.geometry.y))]
glofas_gdf["flow_accum"] = np.clip(glofas_gdf["flow_accum"].values, 1, None)

# === 3. Load stream segments ===
streams = gpd.read_file(stream_file).to_crs(target_crs)
if "segment_id" not in streams.columns:
    streams["segment_id"] = streams.index.astype(str)
targets = streams.copy()
targets["centroid"] = targets.geometry.centroid

# === 4. Process MFD flow direction and build connectivity graph ===
with rasterio.open(flow_dir_raster) as src:
    flow_dir = src.read(1)
    transform = src.transform
    height, width = flow_dir.shape
    inv_transform = ~transform

# Build reverse graph for upstream traversal (using Numba for speed)
def build_reverse_graph(flow_dir):  # No @njit decorator
    height, width = flow_dir.shape
    reverse_graph = {}
    
    # Your existing code will work without Numba restrictions
    dirs = [(-1, 0), (-1, 1), (0, 1), (1, 1), 
            (1, 0), (1, -1), (0, -1), (-1, -1)]
    dir_flags = [64, 128, 1, 2, 4, 8, 16, 32]
    
    for i in range(height):
        for j in range(width):
            current_dir = flow_dir[i, j]
            for (dr, dc), flag in zip(dirs, dir_flags):
                if current_dir == flag:
                    ni, nj = i + dr, j + dc
                    if 0 <= ni < height and 0 <= nj < width:
                        downstream = (ni, nj)
                        if downstream not in reverse_graph:
                            reverse_graph[downstream] = []
                        reverse_graph[downstream].append((i, j))
                    break
    
    return reverse_graph

# Build the reverse flow graph
graph = build_reverse_graph(flow_dir)

# === 5. Map points to raster indices ===
def coords_to_index(x, y):
    """Convert coordinates to raster index key"""
    col, row = inv_transform * (x, y)
    if 0 <= row < height and 0 <= col < width:
        return int(row) * width + int(col)
    return None

# Map GloFAS points to raster indices
glofas_gdf["cell_index"] = [coords_to_index(pt.x, pt.y) for pt in glofas_gdf.geometry]
valid_glofas = glofas_gdf.dropna(subset=["cell_index"])
cell_to_glofas = defaultdict(list)
for idx, cell_idx in enumerate(valid_glofas["cell_index"]):
    cell_to_glofas[cell_idx].append(idx)

# Map target centroids to raster indices
targets["cell_index"] = [coords_to_index(pt.x, pt.y) for pt in targets.centroid]

# === 6. Find upstream cells with memoization ===
from functools import lru_cache

@lru_cache(maxsize=100000)
def find_upstream_cells(start_cell):
    """Find all upstream cells using BFS with memoization"""
    visited = set()
    stack = [start_cell]
    while stack:
        cell = stack.pop()
        if cell in visited:
            continue
        visited.add(cell)
        if cell in graph:
            stack.extend(graph[cell])
    return visited

# === 7. Perform hydrologically-aware IDW interpolation ===
idw_values = []
glofas_points = np.column_stack((valid_glofas.geometry.x, valid_glofas.geometry.y))
Q_arr = valid_glofas["Q"].values
FA_arr = valid_glofas["flow_accum"].values

for i, target in targets.iterrows():
    if pd.isna(target["cell_index"]):
        idw_values.append(0)
        continue
        
    # Get all upstream cells
    upstream_cells = find_upstream_cells(target["cell_index"])
    
    # Collect all GloFAS indices in upstream cells
    candidate_indices = []
    for cell in upstream_cells:
        candidate_indices.extend(cell_to_glofas.get(cell, []))
    
    if not candidate_indices:
        idw_values.append(0)
        continue
        
    # Filter by distance and calculate weights
    pts = glofas_points[candidate_indices]
    dists = np.linalg.norm(pts - np.array([[target.centroid.x, target.centroid.y]]), axis=1)
    valid_mask = dists <= max_dist
    
    if not np.any(valid_mask):
        idw_values.append(0)
        continue
        
    # Apply distance filter
    valid_indices = np.array(candidate_indices)[valid_mask]
    dists = dists[valid_mask]
    q_vals = Q_arr[valid_indices]
    fa_vals = FA_arr[valid_indices]
    
    # Calculate IDW weights
    weights = (fa_vals ** alpha) / (dists ** p)
    total_weight = np.sum(weights)
    
    if total_weight > 0:
        idw_value = np.sum(weights * q_vals) / total_weight
        idw_values.append(idw_value)
    else:
        idw_values.append(0)

# === 8. Assign interpolated discharge to stream segments ===
streams["Q_assigned_mfd_idw"] = idw_values

# === 9. Save output ===
streams.to_file("streams_mfd_weighted.gpkg")