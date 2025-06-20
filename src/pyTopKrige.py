import geopandas as gpd
import rasterio
import numpy as np
import pandas as pd
import xarray as xr
from shapely.geometry import Point
from sklearn.neighbors import BallTree
from scipy.sparse import csc_matrix, eye
from scipy.sparse.linalg import spsolve
from collections import defaultdict
import networkx as nx
import matplotlib.pyplot as plt

# === PARAMETERS ===
glofas_nc = "../Rst/GloFAS_2025_06_13_f.nc"
var_name = "dis24"
time_idx = 0
accum_raster = "../Rst/flow_accumulation.tif"
flow_dir_raster = "../Rst/flowDir.tif"
stream_file = "../geodata/riverQ.gpkg"
target_crs = "EPSG:32719"
max_dist = 15000  # Max distance in meters for considering neighbors
variogram_model = "exponential"  # Can be "exponential", "gaussian", or "spherical"
nugget = 0.1
range_val = 50000  # Range parameter for variogram (in meters)
sill = 1.0  # Sill parameter for variogram

# === 1. Load GloFAS discharge and convert to points ===
ds = xr.open_dataset(glofas_nc)
dis = ds[var_name][time_idx, :, :].isel(forecast_period=0)
if dis.ndim == 3:
    dis2d = dis.isel(forecast_reference_time=0)  # or dis.mean(dim='time') or dis.isel(forecast_reference_time=0)
else:
    dis2d = dis

# Create 2D coordinate grids using meshgrid
lon2d, lat2d = np.meshgrid(dis2d.longitude.values, dis2d.latitude.values)

# Create mask for valid values (must match dis2d.values shape)
mask = (dis2d.values > 0) & np.isfinite(dis2d.values)

# Now apply the mask to the 2D coordinate grids
coords = np.column_stack((lon2d[mask], lat2d[mask]))

glofas_gdf = gpd.GeoDataFrame(
    {'Q': dis2d.values[mask]},
    geometry=gpd.points_from_xy(coords[:, 0], coords[:, 1]),
    crs="EPSG:4326"
).to_crs(target_crs)

# === 2. Extract flow accumulation at GloFAS points ===
with rasterio.open(accum_raster) as src:
    glofas_gdf["flow_accum"] = [val[0] for val in src.sample(zip(glofas_gdf.geometry.x, glofas_gdf.geometry.y))]
glofas_gdf["flow_accum"] = np.clip(glofas_gdf["flow_accum"].values, 1, None)

# === 3. Load stream segments and create river network graph ===
streams = gpd.read_file(stream_file).to_crs(target_crs)
if "segment_id" not in streams.columns:
    streams["segment_id"] = streams.index.astype(str)

# Create directed graph of river network
G = nx.DiGraph()
for idx, seg in streams.iterrows():
    seg_id = seg['segment_id']
    G.add_node(seg_id, geometry=seg.geometry, centroid=seg.geometry.centroid)
    
    # Add edges based on connectivity (simplified - should use actual topology)
    # In practice, you should use the true network topology from your data
    if idx > 0:
        prev_id = streams.loc[idx-1, 'segment_id']
        G.add_edge(prev_id, seg_id)

# Add distance attributes to edges
for u, v in G.edges():
    u_pt = G.nodes[u]['centroid']
    v_pt = G.nodes[v]['centroid']
    G.edges[u, v]['length'] = u_pt.distance(v_pt)

# === 4. Compute hydrological distances ===
def compute_hydrological_distance(G, node1, node2):
    """Compute hydrological distance between two nodes in the river network"""
    try:
        # Find path along the network
        path = nx.shortest_path(G, node1, node2, weight='length')
        dist = 0
        for i in range(len(path)-1):
            dist += G.edges[path[i], path[i+1]]['length']
        return dist
    except nx.NetworkXNoPath:
        return float('inf')

# === 5. Variogram functions ===
def exponential_variogram(h, nugget, range_val, sill):
    return nugget + sill * (1 - np.exp(-h / range_val))

def gaussian_variogram(h, nugget, range_val, sill):
    return nugget + sill * (1 - np.exp(-(h ** 2) / (range_val ** 2)))

def spherical_variogram(h, nugget, range_val, sill):
    if h == 0:
        return 0
    elif h <= range_val:
        return nugget + sill * (1.5 * (h / range_val) - 0.5 * (h / range_val) ** 3)
    else:
        return nugget + sill

VARIOS = {
    "exponential": exponential_variogram,
    "gaussian": gaussian_variogram,
    "spherical": spherical_variogram
}

# === 6. Top-kriging interpolation ===
# Create BallTree for spatial indexing
tree = BallTree(np.vstack([glofas_gdf.geometry.x, glofas_gdf.geometry.y]).T, metric="euclidean")

# Prepare data structures
seg_ids = list(G.nodes())
topkrige_values = []

for seg_id in seg_ids:
    target_pt = G.nodes[seg_id]['centroid']
    
    # Find nearby GloFAS points
    dists, idxs = tree.query([[target_pt.x, target_pt.y]], k=min(50, len(glofas_gdf)))
    dists = dists[0]
    idxs = idxs[0]
    
    valid = dists < max_dist
    if not np.any(valid):
        topkrige_values.append(0)
        continue
        
    # Get valid points
    valid_idxs = idxs[valid]
    valid_dists = dists[valid]
    points = glofas_gdf.iloc[valid_idxs]
    
    # Create covariance matrix
    n = len(points)
    C = np.zeros((n, n))
    vario_func = VARIOS[variogram_model]
    
    # Fill covariance matrix
    for i in range(n):
        for j in range(i, n):
            # Use hydrological distance if possible, else Euclidean
            h_dist = valid_dists[i] + valid_dists[j]  # Simplified - should use actual hydrological distance
            covar = sill - vario_func(h_dist, nugget, range_val, sill)
            C[i, j] = covar
            C[j, i] = covar
    
    # Add nugget to diagonal
    C += nugget * np.eye(n)
    
    # Create right-hand side vector
    b = np.zeros(n)
    for i in range(n):
        # Use hydrological distance to target
        h_dist_target = valid_dists[i]  # Simplified
        covar_target = sill - vario_func(h_dist_target, nugget, range_val, sill)
        b[i] = covar_target
    
    # Solve kriging system
    try:
        weights = np.linalg.solve(C, b)
        weights /= weights.sum()  # Ensure unbiasedness
        q_pred = np.sum(weights * points["Q"].values)
        topkrige_values.append(q_pred)
    except np.linalg.LinAlgError:
        topkrige_values.append(0)

# === 7. Assign interpolated discharge to stream segments ===
for i, seg_id in enumerate(seg_ids):
    G.nodes[seg_id]["Q_topkrige"] = topkrige_values[i]

# Transfer results to GeoDataFrame
streams["Q_topkrige"] = [G.nodes[seg_id]["Q_topkrige"] for seg_id in streams["segment_id"]]

# === 8. Save output ===
streams.to_file("streams_topkrige.gpkg")

# ...existing code...
# --- Fast vectorized IDW (distance only) ---
from scipy.spatial import cKDTree

# Get coordinates
stream_xy = np.column_stack([streams.geometry.centroid.x, streams.geometry.centroid.y])
glofas_xy = np.column_stack([glofas_gdf.geometry.x, glofas_gdf.geometry.y])
glofas_q = glofas_gdf["Q"].values

# Build KDTree for fast neighbor search
tree = cKDTree(glofas_xy)
dists, idxs = tree.query(stream_xy, k=8, distance_upper_bound=MAX_DIST)  # k=8 nearest neighbors

# Compute IDW for each stream segment
p = 2  # IDW power parameter
idw_vals = np.zeros(len(streams))
for i in range(len(streams)):
    valid = np.isfinite(dists[i]) & (dists[i] < MAX_DIST)
    if not np.any(valid):
        idw_vals[i] = 0
        continue
    d = dists[i][valid]
    q = glofas_q[idxs[i][valid]]
    weights = 1.0 / (d ** p)
    weights /= weights.sum()
    idw_vals[i] = np.sum(weights * q)

streams["Q_idw"] = idw_vals

# --- Benchmark both methods ---
def nash_sutcliffe(obs, sim):
    num = np.sum((obs - sim) ** 2)
    den = np.sum((obs - np.mean(obs)) ** 2)
    return 1 - num / den if den != 0 else np.nan

def evaluate_results(streams, glofas_gdf, pred_col, max_dist=MAX_DIST):
    centroids = streams.geometry.centroid
    tree = cKDTree(np.column_stack([glofas_gdf.geometry.x, glofas_gdf.geometry.y]))
    dists, idxs = tree.query(np.column_stack([centroids.x, centroids.y]))
    valid = dists < max_dist
    y_obs = glofas_gdf.iloc[idxs[valid]]["Q"].values
    y_pred = streams.loc[valid, pred_col].values
    if len(y_obs) == 0:
        return {"n_points": 0, "MAE": np.nan, "RMSE": np.nan, "R2": np.nan, "NSE": np.nan}
    mae = mean_absolute_error(y_obs, y_pred)
    rmse = mean_squared_error(y_obs, y_pred, squared=False)
    r2 = r2_score(y_obs, y_pred)
    nse = nash_sutcliffe(y_obs, y_pred)
    return {"n_points": len(y_obs), "MAE": mae, "RMSE": rmse, "R2": r2, "NSE": nse}

metrics_topkrige = evaluate_results(streams, glofas_gdf, pred_col="Q_topkrige")
metrics_idw = evaluate_results(streams, glofas_gdf, pred_col="Q_idw")

print("\n=== BENCHMARK RESULTS ===")
print("Topological Kriging:", metrics_topkrige)
print("IDW:", metrics_idw)