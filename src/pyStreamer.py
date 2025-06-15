import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.enums import Resampling
import xarray as xr
from shapely.geometry import Point, box
from scipy.spatial import cKDTree
from pykrige.ok import OrdinaryKriging

# === PARAMETERS ===
glofas_nc     = "../Rst/GloFAS_2025_06_13_f.nc"
var_name      = "dis24"
time_idx      = 0
forecast_idx  = 0
stream_file   = "../geodata/riverQ.gpkg"
target_crs    = "EPSG:32719"  # Projected CRS for metric interpolation

# === 1. Load GloFAS discharge slice ===
ds = xr.open_dataset(glofas_nc)
dis = ds[var_name][time_idx, :, :].isel(forecast_period=0)
# Detect dims
lat_name = [d for d in dis.dims if 'lat' in d][0]
lon_name = [d for d in dis.dims if 'lon' in d][0]
lats = dis[lat_name].values
lons = dis[lon_name].values

# === 2. Create point GeoDataFrame of non-zero GloFAS discharge ===
points = []
for i, lat in enumerate(lats):
    for j, lon in enumerate(lons):
        q = float(dis.isel(latitude=i, longitude=j).values.item())
        if np.isfinite(q) and q > 0:
            points.append({'geometry': Point(lon, lat), 'Q': q})

glofas_pts = gpd.GeoDataFrame(points, crs="EPSG:4326").to_crs(target_crs)

# Prepare interpolation data
coords = np.vstack([glofas_pts.geometry.x, glofas_pts.geometry.y]).T
values = glofas_pts['Q'].values

# === 3. Load streams and compute centroids ===
streams = gpd.read_file(stream_file).to_crs(target_crs)
streams['centroid'] = streams.geometry.centroid
stream_coords = np.vstack([streams.centroid.x, streams.centroid.y]).T

# === 4. IDW interpolation ===
def idw(xy_query, xy_known, val_known, k=6, power=2):
    tree = cKDTree(xy_known)
    # Remove n_jobs parameter - it's not supported in older scipy versions
    dist, idx = tree.query(xy_query, k=k)
    dist = np.where(dist == 0, 1e-6, dist)
    w = 1 / dist**power
    w_sum = np.sum(w, axis=1)
    return np.sum(w * val_known[idx], axis=1) / w_sum

streams['Q_idw'] = idw(stream_coords, coords, values, k=8, power=2)

# === 5. Ordinary Kriging ===
OK = OrdinaryKriging(
    coords[:,0], coords[:,1], values,
    variogram_model='exponential', verbose=False, enable_plotting=False
)
z, ss = OK.execute('points', stream_coords[:,0], stream_coords[:,1])
streams['Q_krig'] = z

# === 6. Output ===
# Fix 3: Try Shapefile instead
# If you have multiple geometry columns, keep only one
if 'geometry' in streams.columns:
    streams = streams.set_geometry('geometry')
# Drop any other geometry-like columns
for col in streams.columns:
    if col != 'geometry' and streams[col].dtype.name == "geometry":
        streams = streams.drop(columns=[col])
streams.to_file("../Rst/riverQ_id_krig.gpkg", driver="GPKG")

# benchmark
coords = np.vstack([glofas_pts.geometry.x, glofas_pts.geometry.y]).T
values = glofas_pts['Q'].values
n = len(values)
# === 6. Benchmark via Leave-One-Out CV ===
pred_idw = np.zeros(n)
pred_krig = np.zeros(n)
for i in range(n):
    # Training sets
    mask = np.ones(n, dtype=bool)
    mask[i] = False
    coords_tr = coords[mask]
    values_tr = values[mask]
    # IDW
    pred_idw[i] = idw(coords_tr, values_tr, coords[i:i+1], k=8, power=2)[0]
    # Kriging
    OK_lo = OrdinaryKriging(
        coords_tr[:,0], coords_tr[:,1], values_tr,
        variogram_model='exponential', verbose=False, enable_plotting=False
    )
    zi, _ = OK_lo.execute('points', coords[i,0], coords[i,1])
    pred_krig[i] = zi[0]

# Error metrics
rmse_idw  = np.sqrt(np.mean((pred_idw - values)**2))
mae_idw   = np.mean(abs(pred_idw - values))
rmse_krig = np.sqrt(np.mean((pred_krig - values)**2))
mae_krig  = np.mean(abs(pred_krig - values))

print("Benchmark Results:")
print(f"IDW:    RMSE={rmse_idw:.3f}, MAE={mae_idw:.3f}")
print(f"Kriging:RMSE={rmse_krig:.3f}, MAE={mae_krig:.3f}")