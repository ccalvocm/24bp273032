import geopandas as gpd
import xarray as xr
import numpy as np
from sklearn.cluster import MiniBatchKMeans  # Faster alternative to KMeans
from scipy.spatial import cKDTree
from scipy.spatial.distance import pdist, squareform
from scipy.optimize import curve_fit
from sklearn.neighbors import NearestNeighbors
import matplotlib.pyplot as plt
from numba import jit, prange
import warnings
warnings.filterwarnings('ignore')
import time

start= time.time()
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

# VECTORIZED spatial filtering with optimizations
obs_coords = np.column_stack([obs_df.geometry.x.values, obs_df.geometry.y.values])  # Use .values for faster access
n_clusters = len(obs_df) // 2  # Ensure at least 1 cluster

if len(obs_df) > 1:  # Only cluster if we have multiple points
    # Use MiniBatchKMeans for faster clustering on large datasets
    kmeans = MiniBatchKMeans(n_clusters=n_clusters, 
                           random_state=42, 
                           batch_size=1024,  # Process in chunks
                           n_init=3)  # Reduced from 10
    
    # Precompute squared norms for faster Euclidean distance calculation
    kmeans.fit(obs_coords)
    clusters = kmeans.labels_
    
    # VECTORIZED selection using numpy - faster than groupby
    discharge_vals = obs_df['discharge'].values
    max_indices = np.zeros(n_clusters, dtype=int)
    
    for i in range(n_clusters):
        mask = clusters == i
        if np.any(mask):
            max_indices[i] = np.argmax(discharge_vals[mask]) + np.where(mask)[0][0]
    
    # Filter only clusters that had points (in case some clusters are empty)
    valid_clusters = max_indices[max_indices != 0]
    obs_df = obs_df.iloc[valid_clusters].reset_index(drop=True)

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
import numpy as np
from shapely import LineString, MultiLineString

# Optimized prediction point generation
interval = 10000

# Pre-filter valid geometries in one operation
valid_geoms = [geom for geom in network.geometry 
              if geom is not None and not geom.is_empty and geom.length > 0]

# Vectorized length calculations
lengths = np.array([geom.length for geom in valid_geoms])
n_points = np.maximum(2, (lengths // interval).astype(int) + 1)

# Generate all points in one batch operation
all_points = []
for geom, n in zip(valid_geoms, n_points):
    line = geom.exterior if hasattr(geom, 'exterior') else geom
    distances = np.linspace(0, line.length, n)
    all_points.extend(line.interpolate(d) for d in distances)

print(f"Generated {len(all_points)} raw prediction points")

if all_points:
    # 1. VECTORIZED coordinate extraction with numpy fromiter (faster than list comprehension)
    dtype = np.dtype([('x', 'f8'), ('y', 'f8')])
    pred_coords_raw = np.fromiter(((p.x, p.y) for p in all_points), dtype=dtype)
    pred_coords_raw = pred_coords_raw.view('f8').reshape(-1, 2)
    
    # 2. OPTIMIZED duplicate removal using rounding + unique (faster than DBSCAN for this case)
    # Round coordinates to 1m precision (adjust based on your needs)
    rounded_coords = np.round(pred_coords_raw / 100) * 100
    _, unique_indices = np.unique(rounded_coords, axis=0, return_index=True)
    
    # 3. VECTORIZED selection of unique points
    pred_points_unique = [all_points[i] for i in unique_indices]
    
    # 4. Create GeoDataFrame in one operation
    pred_gdf = gpd.GeoDataFrame(geometry=pred_points_unique, crs=network.crs)
    
    # 5. VECTORIZED coordinate extraction and snapping
    pred_coords = np.column_stack([pred_gdf.geometry.x.values, pred_gdf.geometry.y.values])
    _, pred_node_indices = network_tree.query(pred_coords)
    pred_gdf['node'] = pred_node_indices
    
    print(f"After deduplication: {len(pred_gdf)} prediction points")
else:
    pred_gdf = gpd.GeoDataFrame(geometry=[], crs=network.crs).reindex(columns=['node'])
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

pred_positions = np.column_stack([pred_gdf.geometry.x, pred_gdf.geometry.y])

# VECTORIZED nearest neighbor computation
nbrs = NearestNeighbors(n_neighbors=min(20, n_obs), algorithm='ball_tree').fit(obs_positions)
distances, indices = nbrs.kneighbors(pred_positions)

# VECTORIZED prediction
predictions = vectorized_idw_prediction(obs_positions, obs_vals, pred_positions, distances, indices)

pred_gdf['discharge_est'] = predictions

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

end= time.time()
print(f"Total execution time: {end - start:.2f} seconds")