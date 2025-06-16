import geopandas as gpd
import rasterio
import numpy as np
import xarray as xr
from shapely.geometry import Point
from sklearn.neighbors import BallTree
import networkx as nx

# === PARAMETERS ===
GLOFAS_NC = "../Rst/GloFAS_2025_06_13_f.nc"
VAR_NAME = "dis24"
time_idx = 0
accum_raster = "../Rst/flow_accumulation.tif"
STREAM_FILE = "../geodata/riverQ.gpkg"
ELEV_RASTER= "../Rst/dem90fill.tif"
TARGET_CRS = "EPSG:32719"
TIME_IDX = 0
max_dist = 15000
variogram_model = "exponential"
nugget = 0.1
range_val = 50000
sill = 1.0

# === 1. Load GloFAS discharge and convert to points ===
streams = gpd.read_file(STREAM_FILE).to_crs(TARGET_CRS)

# === 2. LOAD GLOFAS DATA AND CONVERT TO POINTS ===
ds = xr.open_dataset(GLOFAS_NC)
dis = ds[VAR_NAME][TIME_IDX, :, :]

# If 3D, reduce to 2D
while dis2d.ndim > 2:
    dis2d = dis2d.isel({dim: 0 for dim in dis2d.dims if dim not in ['latitude', 'longitude']})

lon2d, lat2d = np.meshgrid(dis2d.longitude.values, dis2d.latitude.values)
mask = (dis2d.values > 0) & np.isfinite(dis2d.values)
coords = np.column_stack((lon2d[mask], lat2d[mask]))
Q_vals = dis2d.values[mask]

glofas_gdf = gpd.GeoDataFrame(
    {'Q': Q_vals},
    geometry=gpd.points_from_xy(coords[:, 0], coords[:, 1]),
    crs="EPSG:4326"
).to_crs(TARGET_CRS)

# === 3. EXTRACT COORDINATES AND VALUES ===
coords = np.column_stack([glofas_gdf.geometry.x, glofas_gdf.geometry.y])
values = glofas_gdf['Q'].values
stream_coords = np.column_stack([streams.geometry.centroid.x, streams.geometry.centroid.y])

# === 4. SAMPLE ELEVATION ===
with rasterio.open(ELEV_RASTER) as src:
    pts_elev = np.array([val[0] for val in src.sample(coords)])
    stream_elev = np.array([val[0] for val in src.sample(stream_coords)])

print(f"GloFAS elevation NaNs: {np.isnan(pts_elev).sum()} / {len(pts_elev)}")
print(f"Stream elevation NaNs: {np.isnan(stream_elev).sum()} / {len(stream_elev)}")

# === 5. CLEAN DATA ===
mask = (~np.isnan(values)) & (~np.isnan(pts_elev)) & (values > 0)
coords_clean = coords[mask]
values_clean = values[mask]
elev_clean = pts_elev[mask]

print(f"Valid kriging points: {len(values_clean)}")
print(f"Q variance: {np.var(values_clean):.3f}")
print(f"Elevation variance: {np.var(elev_clean):.3f}")

# === 6. UNIVERSAL KRIGING WITH EXTERNAL DRIFT ===
kriging_success = False
streams['Q_assigned_topkriging'] = np.nan

if len(values_clean) < 10:
    print("Too few valid points for Universal Kriging. Assigning mean Q.")
    streams['Q_assigned_topkriging'] = np.mean(values_clean)
else:
    for model in ['exponential', 'gaussian', 'spherical']:
        try:
            print(f"Trying Universal Kriging with {model} variogram...")
            uk = UniversalKriging(
                coords_clean[:, 0], coords_clean[:, 1], values_clean,
                variogram_model=model,
                drift_terms=['external_Z'],
                external_drift=elev_clean,
                verbose=False
            )
            stream_mask = ~np.isnan(stream_elev)
            z_uk, ss_uk = uk.execute(
                'points',
                stream_coords[stream_mask, 0],
                stream_coords[stream_mask, 1],
                stream_elev[stream_mask]
            )
            streams.loc[stream_mask, 'Q_assigned_topkriging'] = z_uk
            print(f"Universal Kriging successful with {model}.")
            kriging_success = True
            break
        except Exception as e:
            print(f"Universal Kriging failed with {model}: {e}")

if not kriging_success and len(values_clean) >= 10:
    print("Falling back to Ordinary Kriging.")
    try:
        ok = OrdinaryKriging(
            coords_clean[:, 0], coords_clean[:, 1], values_clean,
            variogram_model='exponential',
            verbose=False
        )
        z_ok, ss_ok = ok.execute('points', stream_coords[:, 0], stream_coords[:, 1])
        streams['Q_assigned_topkriging'] = z_ok
        print("Ordinary Kriging fallback successful.")
    except Exception as e:
        print(f"Ordinary Kriging failed: {e}")
        streams['Q_assigned_topkriging'] = np.mean(values_clean)

# === 7. OPTIONAL: ELEVATION-WEIGHTED IDW AS ALTERNATIVE ===
def elevation_weighted_idw(coords_known, vals_known, elev_known, 
                          coords_pred, elev_pred, k=8, power=2, elev_weight=0.5):
    from scipy.spatial import cKDTree
    tree = cKDTree(coords_known)
    dist, idx = tree.query(coords_pred, k=min(k, len(coords_known)))
    if dist.ndim == 1:
        dist = dist.reshape(1, -1)
        idx = idx.reshape(1, -1)
    results = []
    for i in range(len(coords_pred)):
        d = np.maximum(dist[i], 1e-12)
        w_dist = 1.0 / (d ** power)
        elev_diff = np.abs(elev_known[idx[i]] - elev_pred[i])
        w_elev = 1.0 / (1.0 + elev_weight * elev_diff)
        w_total = w_dist * w_elev
        w_norm = w_total / np.sum(w_total)
        pred = np.sum(w_norm * vals_known[idx[i]])
        results.append(pred)
    return np.array(results)

streams['Q_elev_idw'] = elevation_weighted_idw(
    coords_clean, values_clean, elev_clean,
    stream_coords, stream_elev,
    k=8, power=2, elev_weight=0.5
)

# === 8. BENCHMARKING ===
def nash_sutcliffe(obs, sim):
    num = np.sum((obs - sim) ** 2)
    den = np.sum((obs - np.mean(obs)) ** 2)
    return 1 - num / den if den != 0 else np.nan

def evaluate_by_spatial_join(
    streams,
    glofas_gdf,
    pred_col="Q_assigned_topkriging",
    obs_col="Q",
    max_dist=MAX_DIST
):
    s = streams.copy()
    s["centroid"] = s.geometry.centroid
    s = s.set_geometry("centroid")
    s[pred_col] = s[pred_col].fillna(0)
    g = glofas_gdf.copy()
    g = g.set_crs(s.crs, allow_override=True)
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
    y_pred = joined[pred_col].to_numpy()
    y_obs = joined[obs_col].to_numpy()
    mask = np.isfinite(y_pred) & np.isfinite(y_obs)
    y_pred, y_obs = y_pred[mask], y_obs[mask]
    if len(y_obs) == 0:
        print("No valid matched pairs after filtering.")
        return {"n_points": 0, "MAE": np.nan, "RMSE": np.nan, "R2": np.nan, "NSE": np.nan}
    mae = mean_absolute_error(y_obs, y_pred)
    rmse = mean_squared_error(y_obs, y_pred, squared=False)
    r2 = r2_score(y_obs, y_pred)
    nse = nash_sutcliffe(y_obs, y_pred)
    return {"n_points": len(y_obs), "MAE": mae, "RMSE": rmse, "R2": r2, "NSE": nse}

print("\n=== BENCHMARK RESULTS ===")
metrics_topo = evaluate_by_spatial_join(
    streams, glofas_gdf,
    pred_col="Q_assigned_topkriging",
    obs_col="Q", max_dist=MAX_DIST
)
print("Topological Kriging:", metrics_topo)

metrics_elev_idw = evaluate_by_spatial_join(
    streams, glofas_gdf,
    pred_col="Q_elev_idw",
    obs_col="Q", max_dist=MAX_DIST
)
print("Elevation-weighted IDW:", metrics_elev_idw)

# === 9. SAVE OUTPUT ===
streams.to_file("streams_topkrige_output.gpkg")