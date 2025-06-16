import geopandas as gpd
import rasterio
import numpy as np
import xarray as xr
from shapely.geometry import Point
from sklearn.neighbors import BallTree
import networkx as nx

# === PARAMETERS ===
glofas_nc = "../Rst/GloFAS_2025_06_13_f.nc"
var_name = "dis24"
time_idx = 0
accum_raster = "../Rst/flow_accumulation.tif"
stream_file = "../geodata/riverQ.gpkg"
target_crs = "EPSG:32719"
max_dist = 15000
variogram_model = "exponential"
nugget = 0.1
range_val = 50000
sill = 1.0

# === 1. Load GloFAS discharge and convert to points ===
ds = xr.open_dataset(glofas_nc)
dis = ds[var_name][time_idx, :, :].isel(forecast_period=0)
if dis.ndim == 3:
    dis2d = dis.isel(forecast_reference_time=0)
else:
    dis2d = dis

lon2d, lat2d = np.meshgrid(dis2d.longitude.values, dis2d.latitude.values)
mask = (dis2d.values > 0) & np.isfinite(dis2d.values)
coords = np.column_stack((lon2d[mask], lat2d[mask]))

glofas_gdf = gpd.GeoDataFrame(
    {'Q': dis2d.values[mask]},
    geometry=gpd.points_from_xy(coords[:, 0], coords[:, 1]),
    crs="EPSG:4326"
).to_crs(target_crs)

with rasterio.open(accum_raster) as src:
    glofas_gdf["flow_accum"] = [val[0] for val in src.sample(zip(glofas_gdf.geometry.x, glofas_gdf.geometry.y))]
glofas_gdf["flow_accum"] = np.clip(glofas_gdf["flow_accum"].values, 1, None)

# === 2. Load stream segments and create river network graph ===
streams = gpd.read_file(stream_file).to_crs(target_crs)
if "segment_id" not in streams.columns:
    streams["segment_id"] = streams.index.astype(str)

G = nx.DiGraph()
for idx, seg in streams.iterrows():
    seg_id = seg['segment_id']
    G.add_node(seg_id, geometry=seg.geometry, centroid=seg.geometry.centroid)
    if idx > 0:
        prev_id = streams.loc[idx-1, 'segment_id']
        G.add_edge(prev_id, seg_id)

for u, v in G.edges():
    u_pt = G.nodes[u]['centroid']
    v_pt = G.nodes[v]['centroid']
    G.edges[u, v]['length'] = u_pt.distance(v_pt)

# === 3. Variogram functions ===
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

# === 4. Top-kriging interpolation ===
tree = BallTree(np.vstack([glofas_gdf.geometry.x, glofas_gdf.geometry.y]).T, metric="euclidean")
seg_ids = list(G.nodes())
topkrige_values = []

for seg_id in seg_ids:
    target_pt = G.nodes[seg_id]['centroid']
    dists, idxs = tree.query([[target_pt.x, target_pt.y]], k=min(50, len(glofas_gdf)))
    dists = dists[0]
    idxs = idxs[0]
    valid = dists < max_dist
    if not np.any(valid):
        topkrige_values.append(0)
        continue

    valid_idxs = idxs[valid]
    valid_dists = dists[valid]
    points = glofas_gdf.iloc[valid_idxs]
    n = len(points)
    C = np.zeros((n, n))
    vario_func = VARIOS[variogram_model]

    for i in range(n):
        for j in range(i, n):
            h_dist = valid_dists[i] + valid_dists[j]
            covar = sill - vario_func(h_dist, nugget, range_val, sill)
            C[i, j] = covar
            C[j, i] = covar

    C += nugget * np.eye(n)
    b = np.zeros(n)
    for i in range(n):
        h_dist_target = valid_dists[i]
        covar_target = sill - vario_func(h_dist_target, nugget, range_val, sill)
        b[i] = covar_target

    try:
        weights = np.linalg.solve(C, b)
        weights /= weights.sum()
        q_pred = np.sum(weights * points["Q"].values)
        topkrige_values.append(q_pred)
    except np.linalg.LinAlgError:
        topkrige_values.append(0)

for i, seg_id in enumerate(seg_ids):
    G.nodes[seg_id]["Q_topkrige"] = topkrige_values[i]

streams["Q_topkrige"] = [G.nodes[seg_id]["Q_topkrige"] for seg_id in streams["segment_id"]]
streams.to_file("streams_topkrige_full.gpkg")



###benchmark
import numpy as np
import geopandas as gpd
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

def nash_sutcliffe(obs, sim):
    num = np.sum((obs - sim)**2)
    den = np.sum((obs - np.mean(obs))**2)
    return 1 - num/den if den != 0 else np.nan

# --- Benchmark via spatial join ---
def evaluate_by_spatial_join(
    streams,
    glofas_gdf,
    pred_col="Q_assigned_topkriging",
    obs_col="Q",
    max_dist=None  # in projection units, e.g. meters; None = no distance filter
):
    # 1) prepare GeoDataFrames
    s = streams.copy()
    s["centroid"] = s.geometry.centroid
    s = s.set_geometry("centroid")
    
    g = glofas_gdf.copy()
    # ensure both are same CRS
    g = g.set_crs(s.crs, allow_override=True)
    
    # 2) spatial‐join nearest
    joined = gpd.sjoin_nearest(
        s[[pred_col]],
        g[[obs_col, "geometry"]],
        how="inner",
        distance_col="dist",
        max_distance=max_dist
    )
    
    if joined.empty:
        print("No matched points within max_dist.")
        return {"n_points": 0, "MAE": np.nan, "RMSE": np.nan, "R2": np.nan, "NSE": np.nan}
    
    # 3) extract arrays
    y_pred = joined[pred_col].values
    y_obs  = joined[obs_col].values
    
    # drop NaNs/infs if any
    mask = np.isfinite(y_pred) & np.isfinite(y_obs)
    y_pred = y_pred[mask]
    y_obs  = y_obs[mask]
    
    if len(y_obs) == 0:
        print("All matched values are NaN/infinite.")
        return {"n_points": 0, "MAE": np.nan, "RMSE": np.nan, "R2": np.nan, "NSE": np.nan}
    
    # 4) compute metrics
    mae  = mean_absolute_error(y_obs, y_pred)
    rmse = mean_squared_error(y_obs, y_pred, squared=False)
    r2   = r2_score(y_obs, y_pred)
    nse  = nash_sutcliffe(y_obs, y_pred)
    
    return {"n_points": len(y_obs), "MAE": mae, "RMSE": rmse, "R2": r2, "NSE": nse}

# --- Example usage ---

metrics = evaluate_by_spatial_join(
    streams,
    glofas_gdf,
    pred_col="Q_assigned_topkriging",  # your predicted discharge column
    obs_col="Q",                       # GloFAS observed discharge column
    max_dist=1000                      # e.g. only join if within 1 km
)
print(metrics)