import pandas as pd
import geopandas as gpd
import os
import xarray as xr
import rioxarray as rxr
import numpy as np
import pickle
from scipy.interpolate import griddata
from scipy.stats import rankdata
from scipy.interpolate import interp1d

# Compute NSE between ds_nat and df_pivot
def nash_sutcliffe_efficiency(observed, simulated, min_n=6000):
    """Calculate Nash-Sutcliffe Efficiency only if sufficient data points"""
    # Remove NaN values
    mask = ~(np.isnan(observed) | np.isnan(simulated))
    obs = observed[mask]
    sim = simulated[mask]
    
    # Check if we have enough data points
    if len(obs) < min_n:
        return np.nan
    
    # NSE formula
    numerator = np.sum((obs - sim) ** 2)
    denominator = np.sum((obs - np.mean(obs)) ** 2)
    
    if denominator == 0:
        return np.nan
    
    nse = 1 - (numerator / denominator)
    return nse

def train_quantile_delta_mapping(observed, simulated, n_quantiles=100, min_n=1000):
    """
    Train Quantile Delta Mapping (QDM) correction factors.

    Returns a dict with:
      - quantile_levels
      - obs_quantiles
      - sim_quantiles
      - delta_quantiles (obs_quantiles - sim_quantiles)
      - obs_mean, sim_mean
      - n_training_points
    """
    # Drop NaNs
    mask = ~(np.isnan(observed) | np.isnan(simulated))
    obs = observed[mask]
    sim = simulated[mask]

    # Require minimum data
    if len(obs) < min_n:
        return None

    # Compute quantiles
    quantile_levels = np.linspace(0.01, 0.99, n_quantiles)
    obs_q = np.quantile(obs, quantile_levels)
    sim_q = np.quantile(sim, quantile_levels)

    # Store deltas
    delta_q = obs_q - sim_q

    return {
        'quantile_levels': quantile_levels,
        'obs_quantiles':    obs_q,
        'sim_quantiles':    sim_q,
        'delta_quantiles':  delta_q,
        'obs_mean':         np.mean(obs),
        'sim_mean':         np.mean(sim),
        'n_training_points': len(obs)
    }


def interpolate_correction_factors_to_grid(
    station_coords, station_models, grid_coords,
    method='rbf', smooth_sigma=2.0, max_change_factor=3.0
):
    """
    Interpolate station QDM factors to every grid point.
    Adds outlier capping, spatial smoothing, and builds delta_quantiles.
    """
    from scipy.interpolate import Rbf, griddata
    from scipy.ndimage import gaussian_filter

    stations = list(station_models.keys())
    xy = np.array([station_coords[s] for s in stations])
    first = station_models[stations[0]]
    ql = first['quantile_levels']
    nq = len(ql)
    ngrid = len(grid_coords)

    # cap station‐level extremes before interpolation
    for s in stations:
        sim_q = station_models[s]['sim_quantiles']
        obs_q = station_models[s]['obs_quantiles']
        with np.errstate(divide='ignore', invalid='ignore'):
            r = obs_q / (sim_q + 1e-10)
        r = np.clip(r, 1.0/max_change_factor, max_change_factor)
        station_models[s]['obs_quantiles'] = r * sim_q

    # prepare arrays
    grid_obs = np.zeros((ngrid, nq))
    grid_sim = np.zeros((ngrid, nq))

    # interpolate each quantile
    for i in range(nq):
        vals_o = np.array([station_models[s]['obs_quantiles'][i] for s in stations])
        vals_s = np.array([station_models[s]['sim_quantiles'][i] for s in stations])

        try:
            if method == 'rbf':
                r_o = Rbf(xy[:,0], xy[:,1], vals_o, function='multiquadric', smooth=0.1)
                r_s = Rbf(xy[:,0], xy[:,1], vals_s, function='multiquadric', smooth=0.1)
                grid_obs[:,i] = r_o(grid_coords[:,0], grid_coords[:,1])
                grid_sim[:,i] = r_s(grid_coords[:,0], grid_coords[:,1])
            else:
                grid_obs[:,i] = griddata(xy, vals_o, grid_coords, method='nearest')
                grid_sim[:,i] = griddata(xy, vals_s, grid_coords, method='nearest')
        except:
            # fallback to nearest
            grid_obs[:,i] = griddata(xy, vals_o, grid_coords, method='nearest')
            grid_sim[:,i] = griddata(xy, vals_s, grid_coords, method='nearest')

    # build delta_quantiles
    grid_delta = grid_obs - grid_sim

    # interpolate means
    obs_m = np.array([station_models[s]['obs_mean'] for s in stations])
    sim_m = np.array([station_models[s]['sim_mean'] for s in stations])
    if method == 'rbf':
        rm = Rbf(xy[:,0], xy[:,1], obs_m, function='multiquadric', smooth=0.1)
        sm = Rbf(xy[:,0], xy[:,1], sim_m, function='multiquadric', smooth=0.1)
        grid_obs_m = rm(grid_coords[:,0], grid_coords[:,1])
        grid_sim_m = sm(grid_coords[:,0], grid_coords[:,1])
    else:
        grid_obs_m = griddata(xy, obs_m, grid_coords, method='nearest')
        grid_sim_m = griddata(xy, sim_m, grid_coords, method='nearest')

    # spatial smoothing if on regular grid
    if smooth_sigma > 0:
        ux, uy = np.unique(grid_coords[:,0]), np.unique(grid_coords[:,1])
        shape = (len(uy), len(ux))
        if len(grid_coords) == np.prod(shape):
            for arr in (grid_obs, grid_sim, grid_delta):
                tmp = arr.copy()
                for i in range(nq):
                    tmp[:,i] = gaussian_filter(arr[:,i].reshape(shape),
                                               sigma=smooth_sigma, mode='nearest').ravel()
                arr[:] = tmp
            grid_obs_m = gaussian_filter(grid_obs_m.reshape(shape),
                                         sigma=smooth_sigma, mode='nearest').ravel()
            grid_sim_m = gaussian_filter(grid_sim_m.reshape(shape),
                                         sigma=smooth_sigma, mode='nearest').ravel()

    return {
        'quantile_levels':   ql,
        'obs_quantiles':     grid_obs,
        'sim_quantiles':     grid_sim,
        'delta_quantiles':   grid_delta,
        'obs_means':         grid_obs_m,
        'sim_means':         grid_sim_m
    }


def apply_grid_corrections_to_forecast_optimized(forecast_ds, correction_data, chunk_size=50000):
    """
    Apply Quantile Delta Mapping to forecast: corrected = f + δ(f)
    Only non-zero, non-NaN pixels are modified.
    """
    grid_f = correction_data['grid_correction_factors']
    grid_coords = correction_data['grid_coords']
    obs_q = grid_f['obs_quantiles']
    sim_q = grid_f['sim_quantiles']
    dlt_q = grid_f['delta_quantiles']
    ql   = grid_f['quantile_levels']

    # extract forecast array & coords
    if isinstance(forecast_ds, xr.DataArray):
        var = forecast_ds
    else:
        name = list(forecast_ds.data_vars)[0]
        var  = forecast_ds[name]

    if 'lat' in var.dims:
        xs, ys = var.lon.values, var.lat.values
    else:
        xs, ys = var.x.values, var.y.values

    # build kdtree
    xx, yy = np.meshgrid(xs, ys, indexing='xy')
    pts = np.column_stack([xx.ravel(), yy.ravel()])
    from scipy.spatial import cKDTree
    tree = cKDTree(grid_coords)
    _, idx = tree.query(pts, k=1)
    idx2 = idx.reshape(len(ys), len(xs))

    def _corr2d(data2d):
        flat = data2d.ravel()
        ni   = flat.copy()
        m = (~np.isnan(flat)) & (flat > 0)
        vals = flat[m]
        nn   = idx[m]
        # get sim_q and delta_q for each point
        sq = sim_q[nn]
        dq = dlt_q[nn]
        out = flat.copy()
        for i, v in enumerate(vals):
            try:
                d = np.interp(v, sq[i], dq[i],
                              left=dq[i,0], right=dq[i,-1])
                out_idx = np.where(m)[0][i]
                out[out_idx] = max(v + d, 0.0)
            except:
                pass
        out = np.maximum(out, 0.0)  # ensure no negative values
        return out.reshape(data2d.shape)

    # apply per time slice
    ds_corr = forecast_ds.copy(deep=True)
    arr     = var.values
    if 'time' in var.dims:
        for t in range(arr.shape[0]):
            arr[t] = _corr2d(arr[t])
    else:
        arr = _corr2d(arr)
    # assign back
    if isinstance(ds_corr, xr.DataArray):
        ds_corr.values = arr
    else:
        ds_corr[list(ds_corr.data_vars)[0]].values = arr

    return ds_corr

# =============================================================================
# MAIN WORKFLOW: CREATE HISTORICAL CORRECTION FACTORS
# =============================================================================

print("=== CREATING HISTORICAL CORRECTION FACTORS ===")

file=os.path.join('..','output','caudal_diario_historico_4.txt')
df=pd.read_csv(file, sep=',', encoding='latin1')

# metadata
metadata=df.copy()
metadata=pd.pivot_table(metadata, index='NOMBRE ESTACION', values=['UTM_ESTE', 'UTM_NORTE'], aggfunc='first')
gdf= gpd.GeoDataFrame(metadata, geometry=gpd.points_from_xy(metadata['UTM_ESTE'], metadata['UTM_NORTE']), crs='EPSG:32719')
gdf=gdf.to_crs('EPSG:32719')
stream_names=['RIO','ESTERO','QUEBRADA']
gdf_nat=gdf[gdf.index.str.contains('|'.join(stream_names))].copy()
gdf_nat.to_file(os.path.join('..',
                             'geodata',
                             'estaciones_DGA.geojson'),
                              driver='GeoJSON')

# data
df_pivot=pd.pivot_table(df, index='FECHA', columns='NOMBRE ESTACION', values='Caudal_diario', aggfunc='mean')
dates=pd.to_datetime(df_pivot.index, format='%d/%m/%Y', errors='coerce', utc=True, dayfirst=True, yearfirst=False, exact=True, infer_datetime_format=False)
df_pivot.index=dates
df_pivot=df_pivot[df_pivot.columns[df_pivot.columns.isin(gdf_nat.index)]]
df_pivot.to_csv(os.path.join('..',
                            'output',
                            'caudal_diario_historico_4.csv'),
                            encoding='utf-8',
                            index_label='FECHA',
                            date_format='%Y-%m-%dT%H:%M:%S.%fZ',
                            float_format='%.2f')
print("✅ Observed data processed and saved")

# Load historical model data (1980-2018)
path_nc='dis_3d_idw_optimized_1980_2018.nc'
ds=xr.open_dataset(path_nc, chunks={'time': 1, 'lat': 500, 'lon': 500})
ds.rio.write_crs("EPSG:32719", inplace=True)

# sample all ds times over gdf_nat points
gdf_nat = gdf_nat.to_crs(ds.rio.crs)

# Fix: Use DataArrays instead of tuples
xs = gdf_nat.geometry.x.to_numpy()
ys = gdf_nat.geometry.y.to_numpy()

# Create DataArrays with station names as coordinates
station_names = gdf_nat.index.tolist()
x_sel = xr.DataArray(xs, dims="points", coords={"points": station_names})
y_sel = xr.DataArray(ys, dims="points", coords={"points": station_names})

# Select using DataArrays
ds_nat = ds.sel(x=x_sel, y=y_sel, method="nearest")

# Convert to dataframe
ds_nat = ds_nat.to_dataframe().reset_index()

# Fix datetime conversion
ds_nat['FECHA'] = pd.to_datetime(ds_nat['time'])
ds_nat = ds_nat.pivot(index='FECHA', columns='points', values='__xarray_dataarray_variable__')

print("✅ Model data sampled at station locations")

# Save sampled model data
ds_nat.to_csv(os.path.join('..',
                          'output',
                          'dis_3d_idw_optimized_1980_2018.csv'),
                          encoding='utf-8',
                          index_label='FECHA',
                          date_format='%Y-%m-%dT%H:%M:%S.%fZ',
                          float_format='%.2f')

# FIX: Align timezone awareness
ds_nat.index = pd.to_datetime(ds_nat.index).tz_localize('UTC')

print(f"ds_nat dates: {ds_nat.index.min()} to {ds_nat.index.max()}")
print(f"df_pivot dates: {df_pivot.index.min()} to {df_pivot.index.max()}")

# Find common dates and stations
common_dates = ds_nat.index.intersection(df_pivot.index)
common_stations = ds_nat.columns.intersection(df_pivot.columns)
print(f"Common dates: {len(common_dates)}")
print(f"Common stations: {len(common_stations)}")

if len(common_dates) > 0 and len(common_stations) > 0:
    # Subset to common dates and stations
    ds_common = ds_nat.loc[common_dates, common_stations]
    df_common = df_pivot.loc[common_dates, common_stations]
    
    # =============================================================================
    # TRAIN CORRECTION FACTORS AT EACH STATION
    # =============================================================================
    print("\n=== TRAINING CORRECTION FACTORS ===")
    
    correction_factors_stations = {}
    nse_results = {}
    station_coords_dict = {}
    
    for station in common_stations:
        observed = df_common[station].values
        simulated = ds_common[station].values
        
        # Calculate NSE for evaluation
        nse = nash_sutcliffe_efficiency(observed, simulated)
        nse_results[station] = nse
        
        # Train correction factors
        correction_model = train_quantile_delta_mapping(observed, simulated, n_quantiles=100)
        
        if correction_model is not None:
            correction_factors_stations[station] = correction_model
            
            # Store station coordinates
            station_coords_dict[station] = (gdf_nat.loc[station, 'geometry'].x, 
                                          gdf_nat.loc[station, 'geometry'].y)
            
            print(f"✅ {station}: NSE = {nse:.3f}, Correction factors trained")
        else:
            print(f"❌ {station}: NSE = {nse:.3f}, Insufficient data for correction")
    
    print(f"\nTrained correction factors for {len(correction_factors_stations)} stations")
    
    # =============================================================================
    # CREATE MODEL GRID COORDINATES
    # =============================================================================
    print("\n=== CREATING GRID COORDINATES ===")
    
    # Get grid structure from the model dataset
    coord1_vals = ds.y.values  # or lat
    coord2_vals = ds.x.values  # or lon
    grid_x, grid_y = np.meshgrid(coord2_vals, coord1_vals, indexing='xy')
    grid_coords = np.column_stack([grid_x.ravel(), grid_y.ravel()])
    
    print(f"Model grid: {len(coord1_vals)} x {len(coord2_vals)} = {len(grid_coords)} points")
    
    # =============================================================================
    # INTERPOLATE TO FULL GRID
    # =============================================================================
    print("\n=== INTERPOLATING TO FULL GRID ===")
    
    grid_correction_factors = interpolate_correction_factors_to_grid(
        station_coords_dict, 
        correction_factors_stations, 
        grid_coords, 
        method='rbf',             # Changed from 'linear' to 'rbf' to avoid TIN artifacts
        smooth_sigma=2.0,         # Added spatial smoothing
        max_change_factor=3.0     # Added outlier capping
    )
    
    # =============================================================================
    # SAVE HISTORICAL CORRECTION FACTORS
    # =============================================================================
    print("\n=== SAVING CORRECTION FACTORS ===")
    
    correction_data = {
        'grid_correction_factors': grid_correction_factors,
        'grid_coords': grid_coords,
        'grid_shape': (len(coord1_vals), len(coord2_vals)),
        'station_correction_factors': correction_factors_stations,
        'station_coords': station_coords_dict,
        'nse_results': nse_results,
        'training_period': f"{ds_nat.index.min().date()} to {ds_nat.index.max().date()}",
        'creation_date': pd.Timestamp.now(),
        'n_stations': len(correction_factors_stations),
        'n_grid_points': len(grid_coords)
    }
    
    # Save correction factors
    with open('../output/historical_correction_factors.pkl', 'wb') as f:
        pickle.dump(correction_data, f)
    
    # Save NSE results
    nse_df = pd.DataFrame(list(nse_results.items()), columns=['Station', 'NSE'])
    nse_df.to_csv(os.path.join('..', 'output', 'nse_results.csv'), index=False)
    
    # Save summary
    summary_stats = {
        'Training Period': correction_data['training_period'],
        'Stations with Corrections': len(correction_factors_stations),
        'Total Stations Evaluated': len(common_stations),
        'Grid Points': len(grid_coords),
        'Mean NSE': np.nanmean(list(nse_results.values())),
        'Stations with NSE > 0.5': np.sum([nse > 0.5 for nse in nse_results.values() if not np.isnan(nse)])
    }
    
    summary_df = pd.DataFrame(list(summary_stats.items()), columns=['Metric', 'Value'])
    summary_df.to_csv('../output/correction_factors_summary.csv', index=False)
    
    print("✅ HISTORICAL CORRECTION FACTORS CREATED AND SAVED!")
    print(f"   📁 Correction factors: ../output/historical_correction_factors.pkl")
    print(f"   📁 NSE results: ../output/nse_results.csv") 
    print(f"   📁 Summary: ../output/correction_factors_summary.csv")
    print(f"\n🎯 These factors can now be applied to any future forecast!")
    
else:
    print("❌ No common dates or stations found!")
    print("Check date ranges overlap:")
    print(f"ds_nat: {ds_nat.index.min()} to {ds_nat.index.max()}")
    print(f"df_pivot: {df_pivot.index.min()} to {df_pivot.index.max()}")

print("\n=== HISTORICAL CORRECTION FACTORS CREATION COMPLETE ===")

def main():
    """Main function to apply corrections to GloFAS forecast"""
    
    print("=== APPLYING BIAS CORRECTION TO GLOFAS FORECAST ===")
    
    # File paths
    forecast_path = os.path.join('..','Rst','GloFAS_2019_01_01_f.nc')
    correction_factors_path = '../output/historical_correction_factors.pkl'
    
    # Load forecast data
    print(f"Loading forecast: {forecast_path}")
    forecast_ds = xr.open_dataset(forecast_path)
    print(f"✅ Forecast loaded successfully")
    print(f"   Variables: {list(forecast_ds.data_vars)}")
    print(f"   Grid size: {forecast_ds.dims}")

    # Load correction factors
    print(f"Loading correction factors: {correction_factors_path}")
    with open(correction_factors_path, 'rb') as f:
        correction_data = pickle.load(f)
        print(f"✅ Correction factors loaded successfully")
        print(f"   Training period: {correction_data['training_period']}")
        print(f"   Number of stations: {correction_data['n_stations']}")
        print(f"   Grid points: {correction_data['n_grid_points']}")
    
    # Ensure forecast has CRS information
    if not hasattr(forecast_ds, 'rio') or forecast_ds.rio.crs is None:
        print("Adding CRS information to forecast...")
        forecast_ds.rio.write_crs("EPSG:32719", inplace=True)
    
    # Convert to UTM if needed (to match correction factors)
    if str(forecast_ds.rio.crs) != "EPSG:32719":
        print("Reprojecting forecast to UTM 19S...")
        forecast_ds = forecast_ds['__xarray_dataarray_variable__'].rio.reproject("EPSG:32719")
    
    # Apply corrections
    print("\nApplying bias corrections...")
    corrected_ds = apply_grid_corrections_to_forecast_optimized(forecast_ds, correction_data)
    
    # Save corrected forecast
    output_path = forecast_path.replace('.nc', '_bias_corrected.nc')
    print(f"\nSaving corrected forecast: {output_path}")
    
    # Add metadata
    corrected_ds.attrs['bias_correction_applied'] = 'True'  # String instead of boolean
    corrected_ds.attrs['correction_training_period'] = correction_data['training_period']
    corrected_ds.attrs['correction_creation_date'] = str(correction_data['creation_date'])
    corrected_ds.attrs['correction_stations_used'] = str(correction_data['n_stations'])  # Convert to string

    corrected_ds.to_netcdf(output_path,encoding={
        '__xarray_dataarray_variable__': {'dtype': 'float32', 'zlib': True, 'complevel': 5}
    })
    
    print("✅ BIAS CORRECTION COMPLETE!")
    print(f"   📁 Original forecast: {forecast_path}")
    print(f"   📁 Corrected forecast: {output_path}")
    
    # Calculate some summary statistics
    print("\n=== CORRECTION SUMMARY ===")
    var_name = list(forecast_ds.data_vars)[0]
    
    original_mean = float(forecast_ds[var_name].mean())
    corrected_mean = float(corrected_ds[var_name].mean())
    
    print(f"Original mean discharge: {original_mean:.2f}")
    print(f"Corrected mean discharge: {corrected_mean:.2f}")
    print(f"Mean change: {((corrected_mean/original_mean - 1) * 100):.1f}%")
    
    return corrected_ds

if __name__ == "__main__":
    corrected_forecast = main()