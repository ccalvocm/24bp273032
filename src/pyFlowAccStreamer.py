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
accum_raster = "../Rst/flow_accumulation.tif"
stream_file = "../geodata/riverQ.gpkg"
target_crs = "EPSG:32719"
p = 2  # IDW power parameter
alpha = 0.5  # Flow accumulation weight exponent
max_dist = 15000  # Max distance in meters for considering neighbors

# === 1. Load GloFAS discharge and convert to grid polygons ===
ds = xr.open_dataset(glofas_nc)
dis = ds[var_name][time_idx, :, :].isel(forecast_period=0)
lat_name = [dim for dim in dis.dims if 'lat' in dim][0]
lon_name = [dim for dim in dis.dims if 'lon' in dim][0]
lats = dis[lat_name].values
lons = dis[lon_name].values
dy = np.abs(lats[1] - lats[0])
dx = np.abs(lons[1] - lons[0])

# Build GloFAS points
glofas_points = []
for i, y in enumerate(lats):
    for j, x in enumerate(lons):
        q_val = dis.isel({lat_name: i, lon_name: j}).values.item()
        if not np.isnan(q_val):
            pt = Point(x, y)
            glofas_points.append({"geometry": pt, "Q": q_val})

glofas_gdf = gpd.GeoDataFrame(glofas_points, crs="EPSG:4326").to_crs(target_crs)

# === 2. Extract flow accumulation at GloFAS points ===
with rasterio.open(accum_raster) as src:
    glofas_gdf["flow_accum"] = [
        float(next(src.sample([(pt.x, pt.y)]))) for pt in glofas_gdf.geometry
    ]

# Avoid zero or negative accumulation
glofas_gdf["flow_accum"] = glofas_gdf["flow_accum"].clip(lower=1)

# === 3. Load stream segments ===
streams = gpd.read_file(stream_file).to_crs(target_crs)
if "segment_id" not in streams.columns:
    streams["segment_id"] = streams.index.astype(str)

# Use stream centroids as interpolation targets
targets = streams.copy()
targets["centroid"] = targets.geometry.centroid

# === 4. Build BallTree for efficient IDW ===
tree = BallTree(np.vstack([glofas_gdf.geometry.x, glofas_gdf.geometry.y]).T, metric="euclidean")
query_pts = np.vstack([targets.centroid.x, targets.centroid.y]).T
dists, idxs = tree.query(query_pts, k=len(glofas_gdf))

# === 5. Perform IDW with flow accumulation weights ===
idw_values = []
for i, (d, ids) in enumerate(zip(dists, idxs)):
    valid = d < max_dist
    if not np.any(valid):
        idw_values.append(0)
        continue
    d = d[valid]
    ids = ids[valid]

    fa = glofas_gdf.iloc[ids]["flow_accum"].values
    q = glofas_gdf.iloc[ids]["Q"].values

    weights = (fa ** alpha) / (d ** p)
    weights /= weights.sum()
    q_interp = np.sum(weights * q)
    idw_values.append(q_interp)

# === 6. Assign interpolated discharge to stream segments ===
streams["Q_assigned_idw_weighted"] = idw_values

# === 7. Save output ===
streams.to_file("streams_idw_weighted.gpkg")

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
    pred_col="Q_assigned_mfd_idw",
    obs_col="Q",
    max_dist=None
):
    """
    Spatial‐join benchmark between streams (predictions) and GloFAS points (observations).
    streams: GeoDataFrame with pred_col and geometry
    glofas_gdf: GeoDataFrame with obs_col and geometry
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

# === after your IDW assignment ===
# streams["Q_assigned_mfd_idw"] = idw_values

# run benchmark
metrics_idw = evaluate_by_spatial_join(
    streams,
    glofas_gdf,
    pred_col="Q_assigned_idw_weighted",
    obs_col="Q",
    max_dist=max_dist
)
print("IDW benchmark metrics:", metrics_idw)