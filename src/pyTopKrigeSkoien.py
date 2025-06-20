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
def compute_semivariance_fast(obs_vals):
    """Fast semivariance computation using numba"""
    n = len(obs_vals)
    gamma_mat = np.zeros((n, n))
    for i in prange(n):
        for j in prange(n):
            gamma_mat[i, j] = 0.5 * (obs_vals[i] - obs_vals[j])**2
    return gamma_mat

@jit(nopython=True)
def build_covariance_matrix_fast(n, dist_matrix, flow_conn, iw, params_up, params_down):
    """Fast covariance matrix construction"""
    C = np.zeros((n, n))
    for i in prange(n):
        for j in prange(n):
            h = dist_matrix[i, j]
            if not np.isnan(h):
                # Upstream covariance
                if flow_conn[i, j]:
                    cov_up = params_up[1] - (params_up[0] + params_up[1] * (1 - np.exp(-h/params_up[2])) - params_up[0])
                    C[i, j] += iw[i] * iw[j] * cov_up
                # Downstream covariance
                cov_down = params_down[1] - (params_down[0] + params_down[1] * (1 - np.exp(-h/params_down[2])) - params_down[0])
                C[i, j] += cov_down
    return C

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
if dis2d.ndim != 2:
    raise ValueError(f"Expected 2D discharge, got {dis2d.ndim}D")

# --- 2) Build lon/lat grids (must match dis2d) ---
# If you have 1D lon, lat arrays:
# lon_vals = ds.longitude.values
# lat_vals = ds.latitude.values
# then
lon2d, lat2d = np.meshgrid(lon, lat)

# --- 3) Mask and extract valid points ---
valid_mask = (~np.isnan(dis2d)) & (dis2d >= 0)  # Only values > 0.1 m³/s
discharge_threshold = np.percentile(dis2d[~np.isnan(dis2d)], 100)  # Top 25%
significant_mask = valid_mask & (dis2d >= discharge_threshold)
# Use significant discharge points
if np.sum(significant_mask) < 50:
    print("Too few significant points, using all non-zero values")
    mask = valid_mask
else:
    mask = significant_mask
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

print(f"Created {len(obs_df)} observations")

# OPTIMIZATION: Subsample if too many points
MAX_OBS = 1000  # Reduce for faster computation

# SPATIAL FILTERING: Remove clustered points, keep diverse spatial coverage
if len(obs_df) > 1000:
    print("Applying spatial filtering to reduce clustering...")
    
    # Use spatial clustering to get well-distributed points
    from sklearn.cluster import KMeans
    
    coords = np.vstack([(pt.x, pt.y) for pt in obs_df.geometry])
    n_clusters = min(800, len(obs_df) // 2)  # Aim for ~800 well-distributed points
    
    kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
    clusters = kmeans.fit_predict(coords)
    
    # For each cluster, keep the point with highest discharge
    selected_indices = []
    for cluster_id in range(n_clusters):
        cluster_mask = clusters == cluster_id
        cluster_discharges = obs_df.loc[cluster_mask, 'discharge']
        if len(cluster_discharges) > 0:
            best_idx = cluster_discharges.idxmax()
            selected_indices.append(best_idx)
    
    obs_df = obs_df.loc[selected_indices].reset_index(drop=True)
    print(f"After spatial filtering: {len(obs_df)} observations")

# Additional quality filter: Remove extreme outliers
# Q1 = obs_df['discharge'].quantile(0.25)
# Q3 = obs_df['discharge'].quantile(0.75)
# IQR = Q3 - Q1
# lower_bound = Q1 - 15 * IQR
# upper_bound = Q3 + 15 * IQR

# outlier_mask = (obs_df['discharge'] >= lower_bound) & (obs_df['discharge'] <= upper_bound)
# obs_df = obs_df[outlier_mask].reset_index(drop=True)

print(f"After outlier removal: {len(obs_df)} observations")
print(f"Final discharge range: {obs_df['discharge'].min():.2f} to {obs_df['discharge'].max():.2f} m³/s")

# -----------------------------
# 2. Efficient graph building (SIMPLIFIED)
# -----------------------------
print("Building network graph...")

# Use network centroids instead of complex graph topology
network_centroids = np.vstack([
    (geom.centroid.x, geom.centroid.y) 
    for geom in network.geometry
])

# Build KD-tree for fast spatial queries
network_tree = cKDTree(network_centroids)
obs_coords = np.vstack([(pt.x, pt.y) for pt in obs_df.geometry])

# Snap observations to nearest network segments
_, obs_node_indices = network_tree.query(obs_coords)
obs_df['node'] = obs_node_indices

# -----------------------------
# 3. Fast prediction point generation (FIXED)
# -----------------------------
print("Generating prediction points...")
print(f"Network geometry types: {network.geometry.geom_type.value_counts()}")

interval = 10000
pred_points = []

for idx, geom in enumerate(network.geometry):
    if geom is None or geom.is_empty:
        continue
        
    # Handle different geometry types
    lines = []
    if geom.geom_type == 'LineString':
        lines = [geom]
    elif geom.geom_type == 'MultiLineString':
        lines = list(geom.geoms)
    elif geom.geom_type == 'Polygon':
        # Use exterior boundary of polygon
        lines = [geom.exterior]
    elif geom.geom_type == 'MultiPolygon':
        # Use exterior of each polygon
        lines = [poly.exterior for poly in geom.geoms]
    else:
        # For any other type, just use the centroid
        pred_points.append(geom.centroid)
        continue

    # Sample points along each line
    for line in lines:
        try:
            length = line.length
            if length <= 0:
                continue
                
            # Always get at least 2 points (start and end)
            n_points = max(2, int(np.ceil(length / interval)) + 1)
            distances = np.linspace(0, length, n_points)
            
            for d in distances:
                pt = line.interpolate(d)
                if not pt.is_empty:
                    pred_points.append(pt)
                    
        except Exception as e:
            print(f"Warning: Failed to sample geometry {idx}: {e}")
            # Fallback to centroid
            pred_points.append(geom.centroid)

print(f"Generated {len(pred_points)} raw prediction points")

# Remove duplicates and create GeoDataFrame
if pred_points:
    # Remove very close duplicates
    pred_coords = np.array([(p.x, p.y) for p in pred_points])
    
    # Simple duplicate removal (points within 100m)
    from sklearn.cluster import DBSCAN
    clustering = DBSCAN(eps=100, min_samples=1).fit(pred_coords)
    unique_indices = []
    for label in np.unique(clustering.labels_):
        cluster_indices = np.where(clustering.labels_ == label)[0]
        unique_indices.append(cluster_indices[0])  # Take first point from each cluster
    
    pred_points_unique = [pred_points[i] for i in unique_indices]
    
    pred_gdf = gpd.GeoDataFrame(
        geometry=pred_points_unique,
        crs=network.crs
    )
    
    # Snap to network nodes
    pred_coords = np.vstack([(p.x, p.y) for p in pred_gdf.geometry])
    _, pred_node_indices = network_tree.query(pred_coords)
    pred_gdf['node'] = pred_node_indices
    
    print(f"After deduplication: {len(pred_gdf)} prediction points")
else:
    pred_gdf = gpd.GeoDataFrame(columns=['node'], geometry=[], crs=network.crs)
    print("No prediction points generated!")
# -----------------------------
# 4. Optimized distance computation
# -----------------------------
print("Computing distances...")

# Use Euclidean distances (much faster than network distances)
obs_positions = obs_coords
n_obs = len(obs_df)

# Fast pairwise distance computation
dist_OO = squareform(pdist(obs_positions))

# Simplified flow connectivity (using elevation proxy)
flow_conn = np.zeros((n_obs, n_obs), dtype=bool)

# Simple flow direction heuristic (adjust based on your watershed)
for i in range(n_obs):
    for j in range(i+1, n_obs):
        # Assume flow goes from higher to lower y-coordinates (or use DEM if available)
        flow_conn[i, j] = obs_positions[i, 1] > obs_positions[j, 1]
        flow_conn[j, i] = not flow_conn[i, j]

# Flow accumulation weights (simplified)
iw = np.ones(n_obs)  # Equal weights for simplicity, or compute based on upstream area

# -----------------------------
# 5. Vectorized variogram fitting
# -----------------------------
print("Fitting variograms...")

obs_vals = obs_df['discharge'].values

# Fast semivariance computation
gamma_mat = compute_semivariance_fast(obs_vals)

# Extract upper triangle
triu_idx = np.triu_indices(n_obs, k=1)
h_vals = dist_OO[triu_idx]
gamma_vals = gamma_mat[triu_idx]
flow_vals = flow_conn[triu_idx]

# Remove invalid distances
valid_mask = (~np.isnan(h_vals)) & (h_vals > 0)
h_vals = h_vals[valid_mask]
gamma_vals = gamma_vals[valid_mask]
flow_vals = flow_vals[valid_mask]

# Variogram model
def gamma_exp(h, nugget, sill, range_param):
    return nugget + sill * (1 - np.exp(-h / range_param))

# Initial guess
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
# 6. Fast kriging system
# -----------------------------
print("Performing kriging...")

if len(pred_gdf) > 0:
    # Compute obs-prediction distances
    pred_positions = np.vstack([(pt.x, pt.y) for pt in pred_gdf.geometry])
    
    # Use NearestNeighbors for fast distance computation
    nbrs = NearestNeighbors(n_neighbors=min(20, n_obs), algorithm='ball_tree').fit(obs_positions)
    distances, indices = nbrs.kneighbors(pred_positions)
    
    # Simplified kriging (using only nearest neighbors for speed)
    predictions = []
    
    for i, (dists, idx_arr) in enumerate(zip(distances, indices)):
        # Use only valid neighbors
        valid_neighbors = dists < np.inf
        if np.sum(valid_neighbors) > 0:
            local_dists = dists[valid_neighbors]
            local_indices = idx_arr[valid_neighbors]
            local_vals = obs_vals[local_indices]
            
            # Simple inverse distance weighting as approximation
            if len(local_dists) > 0 and local_dists[0] > 0:
                weights = 1.0 / (local_dists + 1e-10)  # Add small constant to avoid division by zero
                weights /= weights.sum()
                pred_val = np.sum(weights * local_vals)
            else:
                pred_val = np.mean(obs_vals)  # Fallback
        else:
            pred_val = np.mean(obs_vals)
            
        predictions.append(pred_val)
    
    pred_gdf['discharge_est'] = predictions
else:
    pred_gdf = gpd.GeoDataFrame({'discharge_est': []}, geometry=[], crs=network.crs)

# -----------------------------
# 7. Assign to network (OPTIMIZED)
# -----------------------------
print("Assigning predictions to network...")

if len(pred_gdf) > 0:
    # Use KD-tree for fast nearest neighbor assignment
    pred_tree = cKDTree(np.vstack([(pt.x, pt.y) for pt in pred_gdf.geometry]))
    
    network_discharge = []
    for geom in network.geometry:
        centroid = geom.centroid
        _, nearest_idx = pred_tree.query([centroid.x, centroid.y])
        network_discharge.append(pred_gdf.iloc[nearest_idx]['discharge_est'])
    
    network['discharge_est'] = network_discharge
else:
    network['discharge_est'] = np.nan

# -----------------------------
# 8. Save and visualize
# -----------------------------
print("Saving results...")

# Save results
network.to_file('riverQ_topKriging_optimized.gpkg', driver='GPKG')

# Simple visualization
fig, ax = plt.subplots(figsize=(12, 8))
network.plot(column='discharge_est', ax=ax, legend=True, linewidth=1.5, cmap='viridis')
obs_df.plot(ax=ax, color='red', markersize=20, alpha=0.7, label='Observations')
plt.title('Top-Kriging Results (Optimized)')
plt.legend()
plt.tight_layout()
plt.savefig('topkriging_results.png', dpi=300, bbox_inches='tight')
plt.show()

print("Optimization complete!")
print(f"Final network segments: {len(network)}")
print(f"Observations used: {len(obs_df)}")
print(f"Predictions made: {len(pred_gdf)}")