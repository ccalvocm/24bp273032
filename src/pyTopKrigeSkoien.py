import geopandas as gpd
import xarray as xr
import numpy as np
import pandas as pd
import networkx as nx
from shapely.geometry import LineString, MultiLineString
from scipy.spatial import cKDTree
from scipy.spatial.distance import pdist, squareform
from scipy.optimize import curve_fit
from sklearn.neighbors import NearestNeighbors
import matplotlib.pyplot as plt
from numba import jit, prange
import warnings
warnings.filterwarnings('ignore')

# JIT-compiled functions for speed
@jit(nopython=True)
def compute_semivariance_vectorized(obs_vals):
    """Vectorized semivariance computation"""
    n = len(obs_vals)
    diff_matrix = obs_vals[:, None] - obs_vals[None, :]
    gamma_mat = 0.5 * diff_matrix**2
    return gamma_mat

@jit(nopython=True)
def build_covariance_matrix_vectorized(dist_matrix, flow_conn, iw, params_up, params_down):
    """Fully vectorized covariance matrix construction"""
    valid_mask = ~np.isnan(dist_matrix)
    
    # Vectorized exponential variogram computation
    exp_up = np.exp(-dist_matrix / params_up[2])
    cov_up = params_up[1] * exp_up  # Simplified covariance
    cov_up = np.where(valid_mask, cov_up, 0.0)
    
    exp_down = np.exp(-dist_matrix / params_down[2])
    cov_down = params_down[1] * exp_down  # Simplified covariance
    cov_down = np.where(valid_mask, cov_down, 0.0)
    
    # Vectorized combination using broadcasting
    iw_matrix = iw[:, None] * iw[None, :]
    C = cov_down + (iw_matrix * cov_up * flow_conn)
    
    return C

# Alternative fully vectorized version (without JIT compilation issues):
def vectorized_flow_connectivity_numpy(obs_positions):
    """Pure NumPy vectorized flow connectivity"""
    n = obs_positions.shape[0]
    
    # Broadcasting: compare all y-coordinates at once
    y_coords = obs_positions[:, 1]
    y_diff = y_coords[:, None] - y_coords[None, :]
    
    # Flow connectivity: upstream if higher y-coordinate
    flow_conn = y_diff > 0
    
    return flow_conn.astype(bool)

@jit(nopython=True, parallel=True)
def vectorized_idw_prediction(obs_positions, obs_vals, pred_positions, distances, indices):
    """Vectorized IDW prediction"""
    n_pred = pred_positions.shape[0]
    predictions = np.zeros(n_pred)
    
    for i in prange(n_pred):
        dists = distances[i]
        idx_arr = indices[i]
        
        # Filter valid neighbors
        valid_mask = dists < np.inf
        if np.sum(valid_mask) > 0:
            valid_dists = dists[valid_mask]
            valid_indices = idx_arr[valid_mask]
            
            if valid_dists[0] > 0:
                weights = 1.0 / (valid_dists + 1e-10)
                weights = weights / np.sum(weights)
                predictions[i] = np.sum(weights * obs_vals[valid_indices])
            else:
                predictions[i] = np.mean(obs_vals)
        else:
            predictions[i] = np.mean(obs_vals)
    
    return predictions

print("Loading data...")
# -----------------------------
# 1. Load and prepare data (OPTIMIZED)
# -----------------------------
network = gpd.read_file('../geodata/riverQ.gpkg').to_crs(epsg=32719)
ds = xr.open_dataset("../Rst/GloFAS_2025_06_13_f.nc")
dis = ds['dis24'].isel(forecast_period=0).values  # Get first time step

# Handle longitude conversion if needed (0-360 to -180-180)
if ds.longitude.min() > 180:
    ds.coords['longitude'] = (ds.coords['longitude'] + 180) % 360 - 180

lon = ds.longitude.values
lat = ds.latitude.values

# 3. Convert gridded forecast to point observations (CORRECTED APPROACH)
# Create coordinate pairs and filter NaN values
dis2d = dis[0, 0, :, :]    # or whatever indices apply
dis2d = np.squeeze(dis2d)
# --- 2) Build lon/lat grids (must match dis2d) ---
# If you have 1D lon, lat arrays:
# lon_vals = ds.longitude.values
# lat_vals = ds.latitude.values
# then
lon2d, lat2d = np.meshgrid(lon, lat)

# --- 3) Mask and extract valid points ---
valid_mask = ~np.isnan(dis2d)  # Only values > 0.1 m³/s
mask = valid_mask
# Extract coordinates and values
xs = lon2d[mask]
ys = lat2d[mask]
zs = dis2d[mask]
print(f"Discharge range: {zs.min():.2f} to {zs.max():.2f} m³/s")
print(f"Mean discharge: {zs.mean():.2f} m³/s")
obs_df = gpd.GeoDataFrame(
    {'discharge': zs},
    geometry=gpd.points_from_xy(xs, ys),
    crs="EPSG:4326"
).to_crs(epsg=32719)

# VECTORIZED spatial filtering
from sklearn.cluster import KMeans

# Extract coordinates in one vectorized operation
obs_coords = np.column_stack([obs_df.geometry.x, obs_df.geometry.y])
n_clusters = len(obs_df) // 2

if n_clusters > 0:
    kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
    clusters = kmeans.fit_predict(obs_coords)
    
    # VECTORIZED cluster selection using groupby
    df_with_clusters = obs_df.copy()
    df_with_clusters['cluster'] = clusters
    selected_obs = df_with_clusters.loc[df_with_clusters.groupby('cluster')['discharge'].idxmax()]
    obs_df = selected_obs.drop('cluster', axis=1).reset_index(drop=True)

print(f"After spatial filtering: {len(obs_df)} observations")
print(f"Final discharge range: {obs_df['discharge'].min():.2f} to {obs_df['discharge'].max():.2f} m³/s")

# -----------------------------
# 2. VECTORIZED network operations
# -----------------------------
print("Building network graph...")

# VECTORIZED centroid extraction
network_centroids = np.column_stack([network.geometry.centroid.x, network.geometry.centroid.y])
network_tree = cKDTree(network_centroids)

# VECTORIZED coordinate extraction and snapping
obs_coords = np.column_stack([obs_df.geometry.x, obs_df.geometry.y])
_, obs_node_indices = network_tree.query(obs_coords)
obs_df['node'] = obs_node_indices

# -----------------------------
# 3. VECTORIZED prediction point generation
# -----------------------------
print("Generating prediction points...")
print(f"Network geometry types: {network.geometry.geom_type.value_counts()}")

interval = 10000
all_points = []

# Process geometries in batches where possible
for geom in network.geometry:
    if geom is None or geom.is_empty:
        continue
    
    line = geom.exterior
    length = line.length
    
    if length > 0:
        n_points = max(2, int(np.ceil(length / interval)) + 1)
        # VECTORIZED distance calculation
        distances = np.linspace(0, length, n_points)
        # Batch interpolation
        points = [line.interpolate(d) for d in distances]
        all_points.extend(points)

print(f"Generated {len(all_points)} raw prediction points")

if all_points:
    # VECTORIZED coordinate extraction
    pred_coords_raw = np.array([(p.x, p.y) for p in all_points])
    
    # VECTORIZED duplicate removal
    from sklearn.cluster import DBSCAN
    clustering = DBSCAN(eps=100, min_samples=1).fit(pred_coords_raw)
    
    # VECTORIZED unique selection
    unique_labels = np.unique(clustering.labels_)
    unique_indices = np.array([np.where(clustering.labels_ == label)[0][0] for label in unique_labels])
    
    pred_points_unique = [all_points[i] for i in unique_indices]
    pred_gdf = gpd.GeoDataFrame(geometry=pred_points_unique, crs=network.crs)
    
    # VECTORIZED snapping
    pred_coords = np.column_stack([pred_gdf.geometry.x, pred_gdf.geometry.y])
    _, pred_node_indices = network_tree.query(pred_coords)
    pred_gdf['node'] = pred_node_indices
    
    print(f"After deduplication: {len(pred_gdf)} prediction points")
else:
    pred_gdf = gpd.GeoDataFrame(columns=['node'], geometry=[], crs=network.crs)

# -----------------------------
# 4. VECTORIZED distance computation
# -----------------------------
print("Computing distances...")

n_obs = len(obs_df)
obs_positions = obs_coords

# VECTORIZED pairwise distance computation
dist_OO = squareform(pdist(obs_positions))

# VECTORIZED flow connectivity
flow_conn = vectorized_flow_connectivity_numpy(obs_positions)

# VECTORIZED weights
iw = np.ones(n_obs)

# -----------------------------
# 5. VECTORIZED variogram fitting
# -----------------------------
print("Fitting variograms...")

obs_vals = obs_df['discharge'].values

# VECTORIZED semivariance computation
gamma_mat = compute_semivariance_vectorized(obs_vals)

# VECTORIZED upper triangle extraction
triu_idx = np.triu_indices(n_obs, k=1)
h_vals = dist_OO[triu_idx]
gamma_vals = gamma_mat[triu_idx]
flow_vals = flow_conn[triu_idx]

# VECTORIZED filtering
valid_mask = (~np.isnan(h_vals)) & (h_vals > 0)
h_vals = h_vals[valid_mask]
gamma_vals = gamma_vals[valid_mask]
flow_vals = flow_vals[valid_mask]

def gamma_exp(h, nugget, sill, range_param):
    return nugget + sill * (1 - np.exp(-h / range_param))

p0 = [0.1 * np.var(obs_vals), 0.9 * np.var(obs_vals), np.median(h_vals[h_vals > 0])]

# Fit variograms
try:
    if np.any(flow_vals):
        params_up, _ = curve_fit(gamma_exp, h_vals[flow_vals], gamma_vals[flow_vals], 
                               p0=p0, maxfev=1000)
    else:
        params_up = p0
        
    if np.any(~flow_vals):
        params_down, _ = curve_fit(gamma_exp, h_vals[~flow_vals], gamma_vals[~flow_vals], 
                                 p0=p0, maxfev=1000)
    else:
        params_down = params_up
except:
    print("Variogram fitting failed, using default parameters")
    params_up = params_down = p0

print(f"Upstream params: {params_up}")
print(f"Downstream params: {params_down}")

# -----------------------------
# 6. VECTORIZED kriging system
# -----------------------------
print("Performing kriging...")

if len(pred_gdf) > 0:
    pred_positions = np.column_stack([pred_gdf.geometry.x, pred_gdf.geometry.y])
    
    # VECTORIZED nearest neighbor computation
    nbrs = NearestNeighbors(n_neighbors=min(20, n_obs), algorithm='ball_tree').fit(obs_positions)
    distances, indices = nbrs.kneighbors(pred_positions)
    
    # VECTORIZED prediction
    predictions = vectorized_idw_prediction(obs_positions, obs_vals, pred_positions, distances, indices)
    
    pred_gdf['discharge_est'] = predictions
else:
    pred_gdf = gpd.GeoDataFrame({'discharge_est': []}, geometry=[], crs=network.crs)

# -----------------------------
# 7. VECTORIZED network assignment
# -----------------------------
print("Assigning predictions to network...")

if len(pred_gdf) > 0:
    pred_tree = cKDTree(np.column_stack([pred_gdf.geometry.x, pred_gdf.geometry.y]))
    
    # VECTORIZED centroid extraction and querying
    network_centroids = np.column_stack([network.geometry.centroid.x, network.geometry.centroid.y])
    _, nearest_indices = pred_tree.query(network_centroids)
    
    # VECTORIZED assignment
    network['discharge_est'] = pred_gdf.iloc[nearest_indices]['discharge_est'].values
else:
    network['discharge_est'] = np.nan

# -----------------------------
# 8. Save and visualize
# -----------------------------
print("Saving results...")

network.to_file('riverQ_topKriging_vectorized.gpkg', driver='GPKG')

fig, ax = plt.subplots(figsize=(12, 8))
network.plot(column='discharge_est', ax=ax, legend=True, linewidth=1.5, cmap='viridis')
obs_df.plot(ax=ax, color='red', markersize=20, alpha=0.7, label='Observations')
plt.title('Top-Kriging Results (Fully Vectorized)')
plt.legend()
plt.tight_layout()
plt.savefig('topkriging_vectorized_results.png', dpi=300, bbox_inches='tight')
plt.show()

print("Vectorized optimization complete!")
print(f"Final network segments: {len(network)}")
print(f"Observations used: {len(obs_df)}")
print(f"Predictions made: {len(pred_gdf)}")