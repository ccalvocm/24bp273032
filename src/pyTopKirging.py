import numpy as np
import xarray as xr
import geopandas as gpd
import rasterio
from rasterio.warp import reproject, Resampling
from shapely.geometry import Point, LineString
from pykrige.ok import OrdinaryKriging
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt

def top_kriging_interpolation():
    """Python implementation of top-kriging for GloFAS discharge forecast"""
    
    # 1. Load libraries (already imported above)
    
    # 2. Read the stream network
    network = gpd.read_file("../geodata/riverQ.gpkg")
    network = network.to_crs(epsg=32719)  # UTM Zone 19S
    
    # 3. Read GloFAS NetCDF forecast
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
    mask = ~np.isnan(dis2d)
    xs = lon2d[mask]
    ys = lat2d[mask]
    zs = dis2d[mask]

    # --- 4) Create GeoDataFrame ---
    obs_df = gpd.GeoDataFrame(
        {'discharge': zs},
        geometry=gpd.points_from_xy(xs, ys),
        crs="EPSG:4326"
    ).to_crs(epsg=32719)

    
    # Create buffers around observation points (1000m radius)
    obs_buffers = obs_df.copy()
    obs_buffers['geometry'] = obs_df.geometry.buffer(1000)
    
    # 5. Generate prediction locations along network
    # Sample points along the network (every 1000m)
    pred_points = []
    for line in network.geometry:
        if line.geom_type == 'LineString':
            length = line.length
            for distance in np.arange(0, length, 1000):
                pred_points.append(line.interpolate(distance))
    
    pred_gdf = gpd.GeoDataFrame(geometry=pred_points, crs=network.crs)
    
    # 6. Prepare data for kriging
    # Get observation points and values
    obs_points = np.array([[p.x, p.y] for p in obs_df.geometry])
    obs_values = obs_df['discharge'].values
    
    # Get prediction points
    pred_points = np.array([[p.x, p.y] for p in pred_gdf.geometry])
    
    # 7. Perform Ordinary Kriging (as Python doesn't have direct top-kriging equivalent)
    # We'll implement network distance-based kriging
    
    # Calculate network distances (simplified - in practice use a network distance matrix)
    # This is a placeholder - actual network distance calculation would use a graph
    if obs_points.ndim == 1 and obs_points.size % 2 == 0:
        obs_points = obs_points.reshape(-1, 2)
    if pred_points.ndim == 1 and pred_points.size % 2 == 0:
        pred_points = pred_points.reshape(-1, 2)
    distances = cdist(obs_points, pred_points, metric='euclidean')
    
    # Fit variogram model (using exponential model)
    variogram_model = {
        'model': 'exponential',
        'psill': 1.0,
        'range': 50000,
        'nugget': 0.1
    }
    
    # Create and fit the kriging model
    ok = OrdinaryKriging(
        obs_points[:, 0], obs_points[:, 1], obs_values,
        variogram_model=variogram_model['model'],
        variogram_parameters=[
            variogram_model['nugget'],
            variogram_model['psill'],
            variogram_model['range']
        ],
        coordinates_type='euclidean'
    )
    
    # 8. Perform kriging prediction
    z, ss = ok.execute('points', pred_points[:, 0], pred_points[:, 1])
    
    # 9. Join predictions back to network
    pred_gdf['estimated'] = z
    
    # Find nearest prediction point for each network segment
    network['estimated'] = network.geometry.apply(
        lambda x: pred_gdf.iloc[pred_gdf.distance(x).argmin()]['estimated']
    )
    
    # 10. Plot results
    fig, ax = plt.subplots(figsize=(12, 8))
    network.plot(ax=ax, color='lightblue', linewidth=2)
    network.plot(ax=ax, column='estimated', linewidth=4, legend=True)
    plt.title("Discharge Estimation Along River Network")
    plt.show()
    
    # 11. Save output
    network.to_file("network_discharge_topkrige.gpkg", driver="GPKG")
    
    return network

# Run the function
result = top_kriging_interpolation()