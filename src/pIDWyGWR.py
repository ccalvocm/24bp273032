import geopandas as gpd
import rasterio
import numpy as np
import pandas as pd
import xarray as xr
from shapely.geometry import Point
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from scipy.spatial import cKDTree
from numba import jit, prange
import time
import warnings
warnings.filterwarnings('ignore')

# === ROBUST NETWORK-BASED GWR FUNCTIONS ===
@jit(nopython=True, parallel=True, cache=True)
def robust_network_gwr_predict(target_coords, target_elevations, obs_coords, obs_elevations, 
                              obs_q_values, bandwidth=10000.0, min_neighbors=5):
    """
    ROBUST Network-based GWR with multiple fallback mechanisms to prevent zeros
    """
    n_targets = target_coords.shape[0]
    n_obs = obs_coords.shape[0]
    predictions = np.zeros(n_targets)
    
    # Global statistics for fallbacks
    global_mean = np.mean(obs_q_values)
    global_median = np.median(obs_q_values)
    
    for i in prange(n_targets):
        target_x, target_y = target_coords[i, 0], target_coords[i, 1]
        target_elev = target_elevations[i]
        
        # Step 1: Calculate distances and weights
        distances = np.zeros(n_obs)
        weights = np.zeros(n_obs)
        
        for j in range(n_obs):
            obs_x, obs_y = obs_coords[j, 0], obs_coords[j, 1]
            obs_elev = obs_elevations[j]
            
            # Euclidean distance
            dx = obs_x - target_x
            dy = obs_y - target_y
            euclidean_dist = np.sqrt(dx*dx + dy*dy)
            distances[j] = euclidean_dist
            
            if euclidean_dist < bandwidth and euclidean_dist > 0:
                # Network-aware distance adjustment
                network_factor = 1.0
                
                if not (np.isnan(obs_elev) or np.isnan(target_elev)):
                    elev_diff = obs_elev - target_elev
                    
                    # Upstream connection (positive elevation difference)
                    if elev_diff > 0:
                        # Favor upstream connections
                        network_factor = 0.6 + 0.4 * np.exp(-elev_diff / 200.0)
                    elif elev_diff < -100:
                        # Penalize far downstream connections
                        network_factor = 1.8
                
                # Combined network distance
                network_dist = euclidean_dist * network_factor
                
                # Gaussian kernel weight
                weights[j] = np.exp(-(network_dist**2) / (2 * bandwidth**2))
        
        # Step 2: Check if we have enough valid neighbors
        valid_weights = weights > 1e-12
        n_valid = np.sum(valid_weights)
        
        if n_valid >= min_neighbors:
            # Try GWR fitting
            try:
                # Normalize weights
                weight_sum = np.sum(weights)
                if weight_sum > 1e-15:
                    normalized_weights = weights / weight_sum
                    
                    # Prepare data for regression
                    X_list = []
                    y_list = []
                    w_list = []
                    
                    for j in range(n_obs):
                        if weights[j] > 1e-12:
                            # Features: [1, rel_x, rel_y, elevation]
                            rel_x = obs_coords[j, 0] - target_x
                            rel_y = obs_coords[j, 1] - target_y
                            obs_elev_safe = obs_elevations[j] if not np.isnan(obs_elevations[j]) else target_elev
                            
                            X_list.append([1.0, rel_x/1000.0, rel_y/1000.0, obs_elev_safe/1000.0])  # Scale for stability
                            y_list.append(obs_q_values[j])
                            w_list.append(np.sqrt(normalized_weights[j]))
                    
                    if len(X_list) >= min_neighbors:
                        # Convert to arrays
                        n_valid_reg = len(X_list)
                        X_reg = np.zeros((n_valid_reg, 4))
                        y_reg = np.zeros(n_valid_reg)
                        w_reg = np.zeros(n_valid_reg)
                        
                        for k in range(n_valid_reg):
                            for l in range(4):
                                X_reg[k, l] = X_list[k][l]
                            y_reg[k] = y_list[k]
                            w_reg[k] = w_list[k]
                        
                        # Weighted design matrix
                        X_weighted = X_reg * w_reg.reshape(-1, 1)
                        y_weighted = y_reg * w_reg
                        
                        # Normal equations: (X'X)β = X'y
                        XTX = np.dot(X_weighted.T, X_weighted)
                        XTy = np.dot(X_weighted.T, y_weighted)
                        
                        # Add regularization for numerical stability
                        reg_lambda = 1e-4
                        XTX_reg = XTX + reg_lambda * np.eye(4)
                        
                        # Check condition number
                        det = np.linalg.det(XTX_reg)
                        
                        if abs(det) > 1e-10:
                            # Solve normal equations
                            beta = np.linalg.solve(XTX_reg, XTy)
                            
                            # Predict at target location (relative coords = 0)
                            target_elev_safe = target_elev if not np.isnan(target_elev) else np.mean(obs_elevations)
                            prediction = beta[0] + beta[3] * (target_elev_safe/1000.0)
                            
                            # Ensure reasonable bounds
                            if prediction > 0 and prediction < 10 * global_mean:
                                predictions[i] = prediction
                                continue
                
                # GWR failed, fallback to weighted average
                weighted_sum = 0.0
                total_weight = 0.0
                
                for j in range(n_obs):
                    if weights[j] > 1e-12:
                        weighted_sum += weights[j] * obs_q_values[j]
                        total_weight += weights[j]
                
                if total_weight > 1e-15:
                    predictions[i] = weighted_sum / total_weight
                    continue
                    
            except:
                pass  # Fall through to next fallback
        
        # Step 3: Fallback to IDW with closest neighbors
        if n_valid >= 3:
            # IDW with valid neighbors
            idw_sum = 0.0
            idw_weights = 0.0
            
            for j in range(n_obs):
                if weights[j] > 1e-12:
                    dist = distances[j]
                    if dist > 0:
                        idw_weight = 1.0 / (dist**2 + 1e-10)
                        idw_sum += idw_weight * obs_q_values[j]
                        idw_weights += idw_weight
            
            if idw_weights > 1e-15:
                predictions[i] = idw_sum / idw_weights
                continue
        
        # Step 4: Final fallback - nearest neighbors
        # Find k-nearest neighbors
        k_nearest = min(5, n_obs)
        nearest_indices = np.argsort(distances)[:k_nearest]
        
        nn_sum = 0.0
        nn_count = 0
        
        for idx_pos in range(k_nearest):
            j = nearest_indices[idx_pos]
            if distances[j] < bandwidth * 2:  # Extended search radius
                nn_sum += obs_q_values[j]
                nn_count += 1
        
        if nn_count > 0:
            predictions[i] = nn_sum / nn_count
        else:
            # Absolute last resort
            predictions[i] = global_median
    
    return predictions

@jit(nopython=True, parallel=True, cache=True)
def simple_network_idw(target_coords, target_elevations, obs_coords, obs_elevations, 
                      obs_q_values, max_dist=15000.0, p=2.0):
    """
    Simplified network-aware IDW as backup method
    """
    n_targets = target_coords.shape[0]
    n_obs = obs_coords.shape[0]
    predictions = np.zeros(n_targets)
    
    global_mean = np.mean(obs_q_values)
    
    for i in prange(n_targets):
        target_x, target_y = target_coords[i, 0], target_coords[i, 1]
        target_elev = target_elevations[i]
        
        total_weight = 0.0
        weighted_sum = 0.0
        
        for j in range(n_obs):
            dx = obs_coords[j, 0] - target_x
            dy = obs_coords[j, 1] - target_y
            dist = np.sqrt(dx*dx + dy*dy)
            
            if dist < max_dist and dist > 0:
                # Network adjustment
                network_factor = 1.0
                if not (np.isnan(obs_elevations[j]) or np.isnan(target_elev)):
                    elev_diff = obs_elevations[j] - target_elev
                    if elev_diff > 0:  # Upstream
                        network_factor = 0.8
                    elif elev_diff < -100:  # Far downstream
                        network_factor = 1.5
                
                adjusted_dist = dist * network_factor
                weight = 1.0 / (adjusted_dist**p)
                
                weighted_sum += weight * obs_q_values[j]
                total_weight += weight
        
        if total_weight > 1e-15:
            predictions[i] = weighted_sum / total_weight
        else:
            predictions[i] = global_mean
    
    return predictions

@jit(nopython=True, cache=True)
def nash_sutcliffe_jit(obs, sim):
    """JIT-compiled Nash-Sutcliffe efficiency"""
    obs_mean = np.mean(obs)
    num = np.sum((obs - sim) ** 2)
    den = np.sum((obs - obs_mean) ** 2)
    return 1.0 - num/den if den != 0 else np.nan

# === PARAMETERS ===
glofas_nc = "../Rst/GloFAS_2025_06_13_f.nc"
var_name = "dis24"
time_idx = 0
ELEV_RASTER= "../Rst/dem90fill.tif"
stream_file = "../geodata/riverQ.gpkg"
target_crs = "EPSG:32719"
p = 2  # IDW power parameter
elev_scale = 0.1  # Scale factor for elevation to match horizontal distance units
max_dist = 15000  # Max distance in meters for considering neighbors

# GWR-specific parameters
gwr_bandwidth = 15000.0  # Increased bandwidth
min_neighbors = 6        # Reduced minimum neighbors

# === 1. Load GloFAS discharge and convert to grid polygons ===
ds = xr.open_dataset(glofas_nc)
dis = ds[var_name][time_idx, :, :].isel(forecast_period=0)
lat_name = [dim for dim in dis.dims if 'lat' in dim][0]
lon_name = [dim for dim in dis.dims if 'lon' in dim][0]
lats = dis[lat_name].values
lons = dis[lon_name].values

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

# Usage:
# glofas_gdf = build_glofas_gdf(dis, lats, lons, target_crs)
glofas_gdf = build_glofas_gdf(dis, lats, lons, target_crs)
# === 2. Extract elevation at GloFAS points ===

# Extract elevations
glofas_coords = np.column_stack([glofas_gdf.geometry.x, glofas_gdf.geometry.y])
with rasterio.open(ELEV_RASTER) as src:
    glofas_gdf["elevation"] = [val[0] for val in src.sample(glofas_coords)]

print(f"Loaded {len(glofas_gdf)} GloFAS observation points")
print(f"GloFAS Q range: {glofas_gdf['Q'].min():.3f} - {glofas_gdf['Q'].max():.3f}")
print(f"GloFAS Q mean: {glofas_gdf['Q'].mean():.3f}")

# === 2. Load stream network ===
streams = gpd.read_file(stream_file).to_crs(target_crs)
if "segment_id" not in streams.columns:
    streams["segment_id"] = streams.index.astype(str)

centroids = streams.geometry.centroid
stream_coords = np.column_stack([centroids.x, centroids.y])

with rasterio.open(ELEV_RASTER) as src:
    stream_elev = np.array([val[0] for val in src.sample(stream_coords)])

print(f"Loaded {len(streams)} stream segments")

# === 3. Data preparation ===
# Remove any zero or negative discharge values that might cause issues
valid_q_mask = glofas_gdf['Q'] > 0.001  # Minimum threshold
if np.sum(~valid_q_mask) > 0:
    print(f"Removing {np.sum(~valid_q_mask)} GloFAS points with Q <= 0.001")
    glofas_gdf = glofas_gdf[valid_q_mask].copy()
    glofas_coords = np.column_stack([glofas_gdf.geometry.x, glofas_gdf.geometry.y])

# Subsample if needed
max_obs_points = 600
if len(glofas_gdf) > max_obs_points:
    print(f"Subsampling GloFAS points: {len(glofas_gdf)} → {max_obs_points}")
    # Stratified sampling to maintain spatial distribution
    sample_indices = np.random.choice(len(glofas_gdf), max_obs_points, replace=False)
    obs_coords = glofas_coords[sample_indices]
    obs_elevations = glofas_gdf['elevation'].values[sample_indices]
    obs_q_values = glofas_gdf['Q'].values[sample_indices]
else:
    obs_coords = glofas_coords
    obs_elevations = glofas_gdf['elevation'].values
    obs_q_values = glofas_gdf['Q'].values

print(f"Using {len(obs_q_values)} observation points")
print(f"Obs Q range: {np.min(obs_q_values):.3f} - {np.max(obs_q_values):.3f}")

# === 4. Try both methods ===
print("Trying Robust Network GWR...")
gwr_start = time.time()

gwr_predictions = robust_network_gwr_predict(
    stream_coords,
    stream_elev,
    obs_coords,
    obs_elevations,
    obs_q_values,
    bandwidth=gwr_bandwidth,
    min_neighbors=min_neighbors
)

print(f"GWR completed in {time.time() - gwr_start:.2f}s")

# Check for zeros and apply backup method if needed
zero_count = np.sum(gwr_predictions <= 0)
if zero_count > len(gwr_predictions) * 0.1:  # If >10% are zero/negative
    print(f"⚠️ GWR produced {zero_count} zeros/negatives ({zero_count/len(gwr_predictions)*100:.1f}%)")
    print("Applying backup Network IDW...")
    
    idw_start = time.time()
    idw_predictions = simple_network_idw(
        stream_coords,
        stream_elev,
        obs_coords,
        obs_elevations,
        obs_q_values,
        max_dist=max_dist,
        p=2.0
    )
    print(f"Backup IDW completed in {time.time() - idw_start:.2f}s")
    
    # Use IDW for zero/negative predictions
    final_predictions = np.where(gwr_predictions <= 0, idw_predictions, gwr_predictions)
    method_used = "GWR + IDW backup"
else:
    final_predictions = gwr_predictions
    method_used = "Pure GWR"

print(f"Method used: {method_used}")
print(f"Final zero count: {np.sum(final_predictions <= 0)}")

# === 5. Assign results ===
streams["Q_network_gwr"] = final_predictions

# Clean and save
streams_clean = streams.copy()
if 'centroid' in streams_clean.columns:
    streams_clean = streams_clean.drop(columns=['centroid'])

streams_clean["Q_network_gwr"] = pd.to_numeric(streams_clean["Q_network_gwr"], errors='coerce')
streams_clean.to_file("streams_robust_network_gwr.gpkg", driver="GPKG")

# === 6. Evaluation ===
def evaluate_results():
    s = streams.copy()
    centroids = s.geometry.centroid
    s = s.set_geometry(centroids)
    s["Q_network_gwr"] = s["Q_network_gwr"].fillna(0)
    
    g = glofas_gdf.copy()
    g = g.set_crs(s.crs, allow_override=True)
    
    joined = gpd.sjoin_nearest(
        s[["Q_network_gwr", s.geometry.name]],
        g[["Q", g.geometry.name]],
        how="inner",
        distance_col="dist",
        max_distance=max_dist
    )
    
    if joined.empty:
        return {"n_points": 0, "MAE": np.nan, "RMSE": np.nan, "R2": np.nan, "NSE": np.nan}
    
    y_pred = joined["Q_network_gwr"].to_numpy()
    y_obs = joined["Q"].to_numpy()
    mask = np.isfinite(y_pred) & np.isfinite(y_obs) & (y_pred > 0) & (y_obs > 0)
    y_pred, y_obs = y_pred[mask], y_obs[mask]
    
    if len(y_obs) == 0:
        return {"n_points": 0, "MAE": np.nan, "RMSE": np.nan, "R2": np.nan, "NSE": np.nan}
    
    mae = mean_absolute_error(y_obs, y_pred)
    rmse = mean_squared_error(y_obs, y_pred, squared=False)
    r2 = r2_score(y_obs, y_pred)
    nse = nash_sutcliffe_jit(y_obs, y_pred)
    
    return {"n_points": len(y_obs), "MAE": mae, "RMSE": rmse, "R2": r2, "NSE": nse}

metrics = evaluate_results()
end_time = time.time()

# === 7. Results ===
print("\n" + "="*60)
print("🌊 ROBUST NETWORK GWR RESULTS")
print("="*60)

print(f"Method: {method_used}")

print("\n=== BENCHMARK RESULTS ===")
print("Robust Network GWR:", metrics)

print("\n=== DISCHARGE STATISTICS ===")
q_pred = streams_clean["Q_network_gwr"].dropna()
print(f"Predicted discharge statistics:")
print(f"  Count: {len(q_pred)}")
print(f"  Non-zero: {np.sum(q_pred > 0)} ({np.sum(q_pred > 0)/len(q_pred)*100:.1f}%)")
print(f"  Min: {q_pred.min():.4f}")
print(f"  Max: {q_pred.max():.4f}")
print(f"  Mean: {q_pred.mean():.4f}")
print(f"  Median: {np.median(q_pred):.4f}")
print(f"  Std Dev: {q_pred.std():.4f}")

print(f"\n=== PARAMETERS ===")
print(f"Bandwidth: {gwr_bandwidth}m")
print(f"Min neighbors: {min_neighbors}")
print(f"Observation points: {len(obs_q_values)}")

print("\n🌊 Robust Network GWR Complete!")