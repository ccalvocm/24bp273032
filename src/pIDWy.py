import geopandas as gpd
import rasterio
import rasterio.mask
import numpy as np
import pandas as pd
import xarray as xr
from shapely.geometry import box, Point
from sklearn.neighbors import BallTree

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

# === 1. Load GloFAS discharge and convert to grid polygons ===
ds = xr.open_dataset(glofas_nc)
dis = ds[var_name][time_idx, :, :].isel(forecast_period=0)
lat_name = [dim for dim in dis.dims if 'lat' in dim][0]
lon_name = [dim for dim in dis.dims if 'lon' in dim][0]
lats = dis[lat_name].values
lons = dis[lon_name].values

# Build GloFAS points
glofas_points = []
for i, y in enumerate(lats):
    for j, x in enumerate(lons):
        q_val = dis.isel({lat_name: i, lon_name: j}).values.item()
        if not np.isnan(q_val):
            pt = Point(x, y)
            glofas_points.append({"geometry": pt, "Q": q_val})

glofas_gdf = gpd.GeoDataFrame(glofas_points, crs="EPSG:4326").to_crs(target_crs)

# === 2. Extract elevation at GloFAS points ===
with rasterio.open(ELEV_RASTER) as src:
    glofas_gdf["elevation"] = [
        float(next(src.sample([(pt.x, pt.y)]))) for pt in glofas_gdf.geometry
    ]

# === 3. Load stream segments ===
streams = gpd.read_file(stream_file).to_crs(target_crs)
if "segment_id" not in streams.columns:
    streams["segment_id"] = streams.index.astype(str)

# Use stream centroids as interpolation targets
targets = streams.copy()
targets["centroid"] = targets.geometry.centroid

# === 4. Extract elevation at stream points ===
stream_coords = np.column_stack([targets.centroid.x, targets.centroid.y])

# Extract elevation at stream centroids
with rasterio.open(ELEV_RASTER) as src:
    stream_elev = np.array([val[0] for val in src.sample(stream_coords)])

# === 5. Create 3D coordinates for BallTree ===
# Scale elevation to match horizontal distance units
glofas_coords_3d = np.column_stack([
    glofas_gdf.geometry.x,
    glofas_gdf.geometry.y,
    glofas_gdf["elevation"] * elev_scale
])

stream_coords_3d = np.column_stack([
    stream_coords[:, 0],
    stream_coords[:, 1],
    stream_elev * elev_scale
])

# Build BallTree with 3D coordinates
tree_3d = BallTree(glofas_coords_3d, metric="euclidean")

# === 6. IDW with 3D distance ===
print("Performing IDW with 3D distance...")

idw_values_3d = []
for i in range(len(stream_coords_3d)):
    # Skip if target elevation is NaN
    if np.isnan(stream_elev[i]):
        # Fallback to 2D distance for this point
        glofas_coords_2d = np.column_stack([glofas_gdf.geometry.x, glofas_gdf.geometry.y])
        tree_2d = BallTree(glofas_coords_2d, metric="euclidean")
        d, ids = tree_2d.query([stream_coords[i]], k=len(glofas_gdf))
        d = d[0]
        ids = ids[0]
    else:
        # Use 3D distance
        d, ids = tree_3d.query([stream_coords_3d[i]], k=len(glofas_gdf))
        d = d[0]
        ids = ids[0]
    
    # Filter by 2D distance for consistency (optional)
    coords_2d_glofas = np.column_stack([glofas_gdf.geometry.x, glofas_gdf.geometry.y])
    coords_2d_neighbors = coords_2d_glofas[ids]
    d_2d = np.sqrt(np.sum((coords_2d_neighbors - stream_coords[i]) ** 2, axis=1))
    valid = d_2d < max_dist
    
    if not np.any(valid):
        idw_values_3d.append(0)
        continue
    
    # Use 3D distances for weighting
    d_3d = d[valid]
    ids_valid = ids[valid]
    
    # Get discharge values for valid neighbors
    q = glofas_gdf.iloc[ids_valid]["Q"].values
    
    # Calculate IDW weights using 3D distance
    weights = 1.0 / (d_3d ** p)
    weights /= weights.sum()
    
    # Interpolated discharge
    q_interp = np.sum(weights * q)
    idw_values_3d.append(q_interp)

# === 7. Assign interpolated discharge to stream segments ===
streams["Q_assigned_3d_idw"] = idw_values_3d

# === 8. Save output ===
# Clean the dataframe before saving
streams_clean = streams.copy()
if 'centroid' in streams_clean.columns:
    streams_clean = streams_clean.drop(columns=['centroid'])

streams_clean["Q_assigned_3d_idw"] = pd.to_numeric(streams_clean["Q_assigned_3d_idw"], errors='coerce')
streams_clean.to_file("streams_3d_idw.gpkg", driver="GPKG")

import numpy as np
import geopandas as gpd
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

def nash_sutcliffe(obs, sim):
    """
    Nash–Sutcliffe Efficiency:
      NSE = 1 - sum((obs - sim)^2) / sum((obs - mean(obs))^2)
    """
    num = np.sum((obs - sim) ** 2)
    den = np.sum((obs - np.mean(obs)) ** 2)
    return 1 - num/den if den != 0 else np.nan

def evaluate_by_spatial_join(
    streams,
    glofas_gdf,
    pred_col="Q_assigned_3d_idw",
    obs_col="Q",
    max_dist=None
):
    """
    Spatial‐join benchmark between streams (predictions) and GloFAS points (observations).
    """
    # ensure the columns exist
    for df, name, col in [
        (streams, "streams", pred_col),
        (glofas_gdf, "glofas_gdf", obs_col)
    ]:
        if col not in df.columns:
            raise KeyError(f"'{col}' not in {name}.columns: {df.columns.tolist()}")

    # prepare centroids for streams
    s = streams.copy()
    s["centroid"] = s.geometry.centroid
    s = s.set_geometry("centroid")
    # fill NaN predictions with zero
    s[pred_col] = s[pred_col].fillna(0)

    # ensure same CRS
    g = glofas_gdf.copy()
    g = g.set_crs(s.crs, allow_override=True)

    # perform nearest‐neighbor spatial join
    joined = gpd.sjoin_nearest(
        s[[pred_col, s.geometry.name]],
        g[[obs_col, g.geometry.name]],
        how="inner",
        distance_col="dist",
        max_distance=max_dist
    )
    if joined.empty:
        print(f"No matches within {max_dist}; aborting benchmark.")
        return {"n_points": 0, "MAE": np.nan, "RMSE": np.nan, "R2": np.nan, "NSE": np.nan}

    # extract arrays and filter finite
    y_pred = joined[pred_col].to_numpy()
    y_obs  = joined[obs_col].to_numpy()
    mask = np.isfinite(y_pred) & np.isfinite(y_obs)
    y_pred, y_obs = y_pred[mask], y_obs[mask]
    if len(y_obs) == 0:
        print("No valid matched pairs after filtering.")
        return {"n_points": 0, "MAE": np.nan, "RMSE": np.nan, "R2": np.nan, "NSE": np.nan}

    # compute metrics
    mae  = mean_absolute_error(y_obs, y_pred)
    rmse = mean_squared_error(y_obs, y_pred, squared=False)
    r2   = r2_score(y_obs, y_pred)
    nse  = nash_sutcliffe(y_obs, y_pred)
    return {"n_points": len(y_obs), "MAE": mae, "RMSE": rmse, "R2": r2, "NSE": nse}

# === 9. Benchmark 3D IDW ===
print("\n=== BENCHMARK RESULTS ===")

metrics_3d_idw = evaluate_by_spatial_join(
    streams,
    glofas_gdf,
    pred_col="Q_assigned_3d_idw",
    obs_col="Q",
    max_dist=max_dist
)
print("3D IDW (using 3D distance):", metrics_3d_idw)

# === 10. Additional statistics ===
print("\n=== ADDITIONAL STATISTICS ===")
q_pred = streams["Q_assigned_3d_idw"].dropna()
print(f"Predicted discharge statistics:")
print(f"  Min: {q_pred.min():.4f}")
print(f"  Max: {q_pred.max():.4f}")
print(f"  Mean: {q_pred.mean():.4f}")
print(f"  Std Dev: {q_pred.std():.4f}")

print(f"\nElev_scale parameter: {elev_scale}")
print(f"This means 1m elevation difference = {elev_scale}m horizontal distance equivalent")