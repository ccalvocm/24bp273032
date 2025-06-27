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
from rasterio.features import rasterize
from rasterio.transform import from_bounds

# JIT-compiled functions for speed
@jit(nopython=True)
def compute_semivariance_vectorized(obs_vals):
    """Vectorized semivariance computation for Top-kriging"""
    n = len(obs_vals)
    diff_matrix = obs_vals[:, None] - obs_vals[None, :]
    gamma_mat = 0.5 * diff_matrix**2
    return gamma_mat

@jit(nopython=True)
def exp_covariance(h, nugget, sill, range_param):
    """Exponential covariance function as in Skøien et al."""
    return np.where(h == 0, sill, (sill - nugget) * np.exp(-h / range_param))

@jit(nopython=True)
def build_topkriging_covariance(dist_matrix, flow_conn, params_up, params_down, area_weights):
    """Build Top-kriging covariance matrix following Skøien et al. methodology"""
    n = dist_matrix.shape[0]
    C = np.zeros((n, n))
    
    nugget_up, sill_up, range_up = params_up[0], params_up[1], params_up[2]
    nugget_down, sill_down, range_down = params_down[0], params_down[1], params_down[2]
    
    for i in range(n):
        for j in range(n):
            h = dist_matrix[i, j]
            if not np.isnan(h):
                # Downstream covariance (base)
                cov_down = exp_covariance(h, nugget_down, sill_down, range_down)
                
                # Upstream covariance (flow-connected)
                if flow_conn[i, j]:
                    cov_up = exp_covariance(h, nugget_up, sill_up, range_up)
                    # Top-kriging formulation: weighted by area ratios
                    C[i, j] = area_weights[i] * area_weights[j] * cov_up + cov_down
                else:
                    C[i, j] = cov_down
    
    return C

print("Loading data...")


def gdf2raster(gdf,
                col='Q_assigned_3d_idw',
                res=500, fill=np.nan,dis_utm=None):
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

# -----------------------------
# 1. Load and prepare data (OPTIMIZED)
# -----------------------------
# === PARAMETERS ===
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

network = gpd.read_file(stream_file).to_crs(epsg=32719)


# Handle longitude conversion if needed (0-360 to -180-180)
if ds.lon.min() > 180:
    ds.coords['lon'] = (ds.coords['lon'] + 180) % 360 - 180

lon = ds.lon.values
lat = ds.lat.values

for time in times:
    print(f"Processing time: {time}")
    # Extract discharge data for the current time
    dis = ds[var_name].sel(time=time)
    # 3. Convert gridded forecast to point observations (CORRECTED APPROACH)
    # Create coordinate pairs and filter NaN values
    # Ensure we have a 2D discharge array by squeezing any singleton dims
    dis2d = np.squeeze(dis.values)
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

    print(f"Created {len(obs_df)} observations")

    # SPATIAL FILTERING: Remove clustered points, keep diverse spatial coverage
    print("Applying spatial filtering to reduce clustering...")

    # Use spatial clustering to get well-distributed points
    from sklearn.cluster import KMeans

    coords = np.vstack([(pt.x, pt.y) for pt in obs_df.geometry])
    n_clusters = len(obs_df)//2   # Aim for ~800 well-distributed points

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

    print(f"After outlier removal: {len(obs_df)} observations")
    print(f"Final discharge range: {obs_df['discharge'].min():.2f} to {obs_df['discharge'].max():.2f} m³/s")

    # -----------------------------
    # 2. Efficient graph building (SIMPLIFIED)
    # -----------------------------
    print("Building network graph...")

    # Create simplified network graph
    G = nx.DiGraph()
    network_centroids = np.vstack([
        (geom.centroid.x, geom.centroid.y) 
        for geom in network.geometry
    ])

    # Add nodes with attributes
    for i, (x, y) in enumerate(network_centroids):
        G.add_node(i, pos=(x, y), area=1.0)  # Simplified area weights

    # Add edges based on proximity (simplified topology)
    network_tree = cKDTree(network_centroids)
    obs_coords = np.vstack([(pt.x, pt.y) for pt in obs_df.geometry])

    # Snap observations to network
    _, obs_node_indices = network_tree.query(obs_coords)
    obs_df['node'] = obs_node_indices

    # -----------------------------
    # 3. Fast prediction point generation (FIXED)
    # -----------------------------
    print("Generating prediction points...")
    print(f"Network geometry types: {network.geometry.geom_type.value_counts()}")

    interval = 5000  # More dense sampling for Top-kriging
    pred_points = []

    for idx, geom in enumerate(network.geometry):
        if geom is None or geom.is_empty:
            continue
        
        # Use exterior boundary for polygon geometries (stream segments)
        line = geom.exterior
        length = line.length
        
        if length > 0:
            n_points = max(3, int(np.ceil(length / interval)) + 1)
            distances = np.linspace(0, length, n_points)
            
            for d in distances:
                pt = line.interpolate(d)
                pred_points.append(pt)

    print(f"Generated {len(pred_points)} prediction points")

    if pred_points:
        pred_coords = np.array([(p.x, p.y) for p in pred_points])
        
        # Remove duplicates
        from sklearn.cluster import DBSCAN
        clustering = DBSCAN(eps=50, min_samples=1).fit(pred_coords)
        unique_indices = []
        for label in np.unique(clustering.labels_):
            cluster_indices = np.where(clustering.labels_ == label)[0]
            unique_indices.append(cluster_indices[0])
        
        pred_points_unique = [pred_points[i] for i in unique_indices]
        pred_gdf = gpd.GeoDataFrame(geometry=pred_points_unique, crs=network.crs)
        
        # Snap to network
        pred_coords = np.vstack([(p.x, p.y) for p in pred_gdf.geometry])
        _, pred_node_indices = network_tree.query(pred_coords)
        pred_gdf['node'] = pred_node_indices
        
        print(f"After deduplication: {len(pred_gdf)} prediction points")
    else:
        pred_gdf = gpd.GeoDataFrame(columns=['node'], geometry=[], crs=network.crs)

    # -----------------------------
    # 4. Compute stream distances and flow connectivity (Top-kriging specific)
    # -----------------------------
    print("Computing stream distances and flow connectivity...")

    n_obs = len(obs_df)
    obs_positions = obs_coords

    # Use Euclidean distances as proxy for stream distances
    dist_OO = squareform(pdist(obs_positions))

    # Flow connectivity matrix (key component of Top-kriging)
    flow_conn = np.zeros((n_obs, n_obs), dtype=bool)

    # Simplified flow direction: upstream if higher y-coordinate
    for i in range(n_obs):
        for j in range(n_obs):
            if i != j:
                # Point i flows to point j if i has higher y-coordinate (elevation proxy)
                flow_conn[i, j] = obs_positions[i, 1] > obs_positions[j, 1]

    # Compute drainage area weights (simplified)
    area_weights = np.ones(n_obs)  # Equal weights - replace with actual drainage areas if available
    for i in range(n_obs):
        # Count upstream points as proxy for drainage area
        upstream_count = np.sum(flow_conn[:, i])
        area_weights[i] = np.sqrt(upstream_count + 1)  # Simple area weighting

    # Normalize weights
    area_weights = area_weights / np.max(area_weights)

    print(f"Computed flow connectivity for {n_obs} observations")

    # -----------------------------
    # 5. Fit directional variograms (Top-kriging)
    # -----------------------------
    print("Fitting directional variograms for Top-kriging...")

    obs_vals = obs_df['discharge'].values
    gamma_mat = compute_semivariance_vectorized(obs_vals)

    # Extract upper triangle
    triu_idx = np.triu_indices(n_obs, k=1)
    h_vals = dist_OO[triu_idx]
    gamma_vals = gamma_mat[triu_idx]
    flow_vals = flow_conn[triu_idx]

    # Filter valid distances
    valid_mask = (~np.isnan(h_vals)) & (h_vals > 0)
    h_vals = h_vals[valid_mask]
    gamma_vals = gamma_vals[valid_mask]
    flow_vals = flow_vals[valid_mask]

    def gamma_exp(h, nugget, sill, range_param):
        return nugget + sill * (1 - np.exp(-h / range_param))

    # Initial parameters
    p0 = [0.1 * np.var(obs_vals), 0.9 * np.var(obs_vals), np.median(h_vals)]

    # Fit upstream and downstream variograms
    try:
        if np.any(flow_vals):
            params_up, _ = curve_fit(gamma_exp, h_vals[flow_vals], gamma_vals[flow_vals], 
                                p0=p0, maxfev=2000)
        else:
            params_up = p0
            
        if np.any(~flow_vals):
            params_down, _ = curve_fit(gamma_exp, h_vals[~flow_vals], gamma_vals[~flow_vals], 
                                    p0=p0, maxfev=2000)
        else:
            params_down = params_up
    except:
        print("Variogram fitting failed, using default parameters")
        params_up = params_down = p0

    print(f"Upstream variogram: nugget={params_up[0]:.3f}, sill={params_up[1]:.3f}, range={params_up[2]:.1f}")
    print(f"Downstream variogram: nugget={params_down[0]:.3f}, sill={params_down[1]:.3f}, range={params_down[2]:.1f}")

    # -----------------------------
    # 6. Top-kriging prediction system
    # -----------------------------
    print("Performing Top-kriging...")

    if len(pred_gdf) > 0:
        pred_positions = np.vstack([(pt.x, pt.y) for pt in pred_gdf.geometry])
        
        # Build observation-observation covariance matrix
        C_oo = build_topkriging_covariance(dist_OO, flow_conn, 
                                        np.array(params_up), np.array(params_down), 
                                        area_weights)
        
        # Build kriging system matrix (with unbiasedness constraint)
        K = np.zeros((n_obs + 1, n_obs + 1))
        K[:n_obs, :n_obs] = C_oo
        K[:n_obs, n_obs] = 1.0  # Constraint row
        K[n_obs, :n_obs] = 1.0  # Constraint column
        K[n_obs, n_obs] = 0.0   # Corner
        
        # Regularize for numerical stability
        for i in range(n_obs):
            K[i, i] += 1e-6
        
        try:
            K_inv = np.linalg.inv(K)
            print("Successfully inverted kriging system matrix")
        except:
            print("Matrix inversion failed, using pseudoinverse")
            K_inv = np.linalg.pinv(K)
        
        # Compute predictions
        predictions = []
        
        for pred_pos in pred_positions:
            # Compute distances from prediction point to observations
            pred_dists = np.sqrt(np.sum((obs_positions - pred_pos)**2, axis=1))
            
            # Compute flow connectivity from obs to prediction point
            pred_flow_conn = np.array([obs_positions[i, 1] > pred_pos[1] for i in range(n_obs)])
            
            # Build covariance vector (obs to prediction point)
            c_op = np.zeros(n_obs)
            for i in range(n_obs):
                h = pred_dists[i]
                # Downstream covariance
                cov_down = exp_covariance(h, params_down[0], params_down[1], params_down[2])
                
                if pred_flow_conn[i]:
                    # Upstream covariance
                    cov_up = exp_covariance(h, params_up[0], params_up[1], params_up[2])
                    c_op[i] = area_weights[i] * cov_up + cov_down
                else:
                    c_op[i] = cov_down
            
            # Augment with constraint
            rhs = np.zeros(n_obs + 1)
            rhs[:n_obs] = c_op
            rhs[n_obs] = 1.0
            
            # Solve for weights
            weights = K_inv @ rhs
            
            # Compute prediction
            pred_val = np.sum(weights[:n_obs] * obs_vals)
            predictions.append(pred_val)
        
        pred_gdf['discharge_est'] = predictions
        
        print(f"Top-kriging completed. Prediction range: {np.min(predictions):.2f} to {np.max(predictions):.2f}")

    else:
        pred_gdf = gpd.GeoDataFrame({'discharge_est': []}, geometry=[], crs=network.crs)

    # -----------------------------
    # 7. Assign predictions to network
    # -----------------------------
    print("Assigning Top-kriging predictions to network...")

    if len(pred_gdf) > 0:
        pred_tree = cKDTree(np.vstack([(pt.x, pt.y) for pt in pred_gdf.geometry]))
        
        network_discharge = []
        for geom in network.geometry:
            centroid = geom.centroid
            _, nearest_idx = pred_tree.query([centroid.x, centroid.y])
            network_discharge.append(pred_gdf.iloc[nearest_idx]['discharge_est'])
        
        network['discharge_est'] = network_discharge
    else:
        network['discharge_est'] = np.nan

    ras = gdf2raster(network,
                col='discharge_est',
                res=500, fill=np.nan,
                dis_utm=dis.rio.write_crs('EPSG:4326').rio.reproject(target_crs))

    ras_list.append(ras)

# === 9. Save results ===
ras_combined = xr.concat(ras_list, dim='time')
ras_combined = ras_combined.assign_coords(
    time=pd.to_datetime(times)).rename({'time': 'time'})
ras_combined.rio.write_crs(target_crs, inplace=True)
ras_combined.to_netcdf("dis_Top_Krige_1980_2018.nc")