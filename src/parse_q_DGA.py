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
def nash_sutcliffe_efficiency(observed, simulated, min_n=12000):
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

def train_quantile_delta_mapping(observed, simulated, n_quantiles=100):
    """
    Train quantile delta mapping correction factors
    
    Parameters:
    -----------
    observed : array-like
        Historical observed data at station
    simulated : array-like  
        Historical simulated data at same station
    n_quantiles : int
        Number of quantiles for mapping
    
    Returns:
    --------
    dict : Correction mapping or None if insufficient data
    """
    # Remove NaN values
    mask = ~(np.isnan(observed) | np.isnan(simulated))
    obs_clean = observed[mask]
    sim_clean = simulated[mask]
    
    # Need sufficient data points
    if len(obs_clean) < 1000:  # Minimum threshold
        return None
    
    # Calculate quantiles
    quantile_levels = np.linspace(0.01, 0.99, n_quantiles)
    
    obs_quantiles = np.quantile(obs_clean, quantile_levels)
    sim_quantiles = np.quantile(sim_clean, quantile_levels)
    
    # Create correction mapping
    correction_factors = {
        'quantile_levels': quantile_levels,
        'obs_quantiles': obs_quantiles,
        'sim_quantiles': sim_quantiles,
        'obs_mean': np.mean(obs_clean),
        'sim_mean': np.mean(sim_clean),
        'n_training_points': len(obs_clean)
    }
    
    return correction_factors

def interpolate_correction_factors_to_grid(station_coords, station_corrections, grid_coords, method='linear', smooth_sigma=2.0, max_change_factor=3.0):
    """
    Interpolate station-based correction factors to full model grid with improvements
    
    Parameters:
    -----------
    station_coords : dict
        {station_name: (x, y)} coordinates of stations
    station_corrections : dict
        {station_name: correction_factors} from training
    grid_coords : array
        Grid coordinates as [(x1,y1), (x2,y2), ...]
    method : str
        Interpolation method ('linear', 'cubic', 'nearest')
    smooth_sigma : float
        Gaussian smoothing parameter (0 = no smoothing)
    max_change_factor : float
        Maximum allowed correction factor (3.0 = max 3x change)
        
    Returns:
    --------
    dict : Grid-based correction factors with artifacts reduced
    """
    from scipy.ndimage import gaussian_filter
    
    print(f"Interpolating {len(station_corrections)} station corrections to {len(grid_coords)} grid points...")
    
    # Extract station coordinates and data
    valid_stations = list(station_corrections.keys())
    station_xy = np.array([station_coords[station] for station in valid_stations])
    
    # Get quantile structure from first station
    first_correction = list(station_corrections.values())[0]
    quantile_levels = first_correction['quantile_levels']
    n_quantiles = len(quantile_levels)
    
    # IMPROVEMENT 3: Add outlier detection and capping for station data
    print("Detecting and capping outlier correction factors...")
    for station in valid_stations:
        obs_q = station_corrections[station]['obs_quantiles']
        sim_q = station_corrections[station]['sim_quantiles']
        
        # Calculate correction ratios and cap extreme values
        with np.errstate(divide='ignore', invalid='ignore'):
            ratios = obs_q / (sim_q + 1e-10)
            
        # Cap extreme correction factors
        ratios = np.clip(ratios, 1.0/max_change_factor, max_change_factor)
        
        # Recalculate corrected observed quantiles
        station_corrections[station]['obs_quantiles'] = ratios * sim_q
    
    # Initialize arrays for grid correction factors
    grid_correction_factors = {
        'quantile_levels': quantile_levels,
        'obs_quantiles': np.zeros((len(grid_coords), n_quantiles)),
        'sim_quantiles': np.zeros((len(grid_coords), n_quantiles))
    }
    
    # IMPROVEMENT 1: Use linear interpolation with nearest neighbor fallback
    print(f"Using {method} interpolation with nearest neighbor fallback...")
    
    # Interpolate each quantile level
    for i, q_level in enumerate(quantile_levels):
        # Extract observed and simulated quantiles at this level from all stations
        obs_q_values = np.array([station_corrections[station]['obs_quantiles'][i] 
                                for station in valid_stations])
        sim_q_values = np.array([station_corrections[station]['sim_quantiles'][i] 
                                for station in valid_stations])
        
        # Primary interpolation using specified method
        try:
            grid_obs_q = griddata(station_xy, obs_q_values, grid_coords, 
                                 method=method, fill_value=np.nan)
            grid_sim_q = griddata(station_xy, sim_q_values, grid_coords, 
                                 method=method, fill_value=np.nan)
            
            # Fill NaN values with nearest neighbor as fallback
            if np.any(np.isnan(grid_obs_q)):
                nan_mask = np.isnan(grid_obs_q)
                grid_obs_q_nearest = griddata(station_xy, obs_q_values, grid_coords[nan_mask], 
                                             method='nearest')
                grid_sim_q_nearest = griddata(station_xy, sim_q_values, grid_coords[nan_mask], 
                                             method='nearest')
                grid_obs_q[nan_mask] = grid_obs_q_nearest
                grid_sim_q[nan_mask] = grid_sim_q_nearest
                
        except Exception as e:
            print(f"Warning: {method} interpolation failed for quantile {i}, using nearest neighbor. Error: {e}")
            grid_obs_q = griddata(station_xy, obs_q_values, grid_coords, method='nearest', fill_value=np.nan)
            grid_sim_q = griddata(station_xy, sim_q_values, grid_coords, method='nearest', fill_value=np.nan)
        
        # Store in structure
        grid_correction_factors['obs_quantiles'][:, i] = grid_obs_q
        grid_correction_factors['sim_quantiles'][:, i] = grid_sim_q
    
    # Also interpolate means with same approach
    obs_means = np.array([station_corrections[station]['obs_mean'] for station in valid_stations])
    sim_means = np.array([station_corrections[station]['sim_mean'] for station in valid_stations])
    
    try:
        grid_obs_means = griddata(station_xy, obs_means, grid_coords, method=method, fill_value=np.nan)
        grid_sim_means = griddata(station_xy, sim_means, grid_coords, method=method, fill_value=np.nan)
        
        # Fill NaN values with nearest neighbor
        if np.any(np.isnan(grid_obs_means)):
            nan_mask = np.isnan(grid_obs_means)
            grid_obs_means[nan_mask] = griddata(station_xy, obs_means, grid_coords[nan_mask], method='nearest')
            grid_sim_means[nan_mask] = griddata(station_xy, sim_means, grid_coords[nan_mask], method='nearest')
            
    except Exception as e:
        print(f"Warning: {method} interpolation failed for means, using nearest neighbor. Error: {e}")
        grid_obs_means = griddata(station_xy, obs_means, grid_coords, method='nearest', fill_value=np.nan)
        grid_sim_means = griddata(station_xy, sim_means, grid_coords, method='nearest', fill_value=np.nan)
    
    grid_correction_factors['obs_means'] = grid_obs_means
    grid_correction_factors['sim_means'] = grid_sim_means
    
    # IMPROVEMENT 2: Apply spatial smoothing to reduce artifacts
    if smooth_sigma > 0:
        print(f"Applying spatial smoothing (sigma={smooth_sigma})...")
        
        # Need to determine grid shape for smoothing
        # Estimate grid shape from coordinates
        unique_x = np.unique(grid_coords[:, 0])
        unique_y = np.unique(grid_coords[:, 1])
        grid_shape = (len(unique_y), len(unique_x))
        
        if len(grid_coords) == np.prod(grid_shape):
            # Regular grid - can apply smoothing
            for i in range(n_quantiles):
                # Reshape to 2D grid for smoothing
                obs_2d = grid_correction_factors['obs_quantiles'][:, i].reshape(grid_shape)
                sim_2d = grid_correction_factors['sim_quantiles'][:, i].reshape(grid_shape)
                
                # Apply Gaussian smoothing
                obs_smooth = gaussian_filter(obs_2d, sigma=smooth_sigma, mode='nearest')
                sim_smooth = gaussian_filter(sim_2d, sigma=smooth_sigma, mode='nearest')
                
                # Put back
                grid_correction_factors['obs_quantiles'][:, i] = obs_smooth.ravel()
                grid_correction_factors['sim_quantiles'][:, i] = sim_smooth.ravel()
            
            # Smooth means too
            obs_means_2d = grid_correction_factors['obs_means'].reshape(grid_shape)
            sim_means_2d = grid_correction_factors['sim_means'].reshape(grid_shape)
            
            obs_means_smooth = gaussian_filter(obs_means_2d, sigma=smooth_sigma, mode='nearest')
            sim_means_smooth = gaussian_filter(sim_means_2d, sigma=smooth_sigma, mode='nearest')
            
            grid_correction_factors['obs_means'] = obs_means_smooth.ravel()
            grid_correction_factors['sim_means'] = sim_means_smooth.ravel()
            
            print("✅ Spatial smoothing applied")
        else:
            print("⚠️  Irregular grid detected - skipping spatial smoothing")
    
    # IMPROVEMENT 3: Final check for extreme correction factors
    print("Final check for extreme correction factors...")
    for i in range(n_quantiles):
        obs_q = grid_correction_factors['obs_quantiles'][:, i]
        sim_q = grid_correction_factors['sim_quantiles'][:, i]
        
        # Calculate and cap correction ratios
        with np.errstate(divide='ignore', invalid='ignore'):
            ratios = obs_q / (sim_q + 1e-10)
            
        # Cap extreme ratios
        ratios_capped = np.clip(ratios, 1.0/max_change_factor, max_change_factor)
        
        # Count how many were capped
        n_capped = np.sum(~np.isclose(ratios, ratios_capped, equal_nan=True))
        if n_capped > 0:
            print(f"  Quantile {i}: Capped {n_capped} extreme correction factors")
        
        # Update with capped values
        grid_correction_factors['obs_quantiles'][:, i] = ratios_capped * sim_q
    
    print("✅ Grid interpolation complete with artifact reduction!")
    return grid_correction_factors

def apply_quantile_delta_mapping(forecast_values, correction_factors):
    """
    Apply quantile delta mapping correction to forecast values
    
    Parameters:
    -----------
    forecast_values : array-like or scalar
        Raw forecast values to be corrected
    correction_factors : dict
        Trained correction factors from historical period
        
    Returns:
    --------
    array or scalar : Bias-corrected forecast values
    """
    # Convert to numpy array if needed
    forecast_values = np.asarray(forecast_values)
    
    # Handle scalar case
    if forecast_values.ndim == 0:
        forecast_values = forecast_values.reshape(1)
        return_scalar = True
    else:
        return_scalar = False
    
    # Remove NaN values for processing
    mask = ~np.isnan(forecast_values)
    if not np.any(mask):
        return forecast_values if not return_scalar else forecast_values[0]
    
    clean_forecast = forecast_values[mask]
    corrected_clean = np.zeros_like(clean_forecast)
    
    # Get correction mapping
    quantile_levels = correction_factors['quantile_levels']
    obs_quantiles = correction_factors['obs_quantiles'] 
    sim_quantiles = correction_factors['sim_quantiles']
    
    # Handle edge cases
    if len(clean_forecast) == 0:
        return forecast_values if not return_scalar else forecast_values[0]
    
    # Create interpolation functions
    try:
        from scipy.interpolate import interp1d
        
        # Map simulated quantiles to observed quantiles
        sim_to_obs_interp = interp1d(
            sim_quantiles, obs_quantiles, 
            kind='linear', bounds_error=False, 
            fill_value=(obs_quantiles[0], obs_quantiles[-1])
        )
        
        # Apply correction
        corrected_clean = sim_to_obs_interp(clean_forecast)
        
        # Ensure no negative values for streamflow
        corrected_clean = np.maximum(corrected_clean, 0.0)
        
    except Exception as e:
        print(f"Warning: Correction failed, returning original values. Error: {e}")
        corrected_clean = clean_forecast
    
    # Put corrected values back
    corrected_forecast = forecast_values.copy()
    corrected_forecast[mask] = corrected_clean
    
    return corrected_forecast[0] if return_scalar else corrected_forecast

def apply_grid_corrections_to_forecast_optimized(forecast_ds, correction_data, chunk_size=50000):
    """
    Optimized version using vectorized operations and pre-computed nearest neighbors
    """
    print("Applying bias corrections to forecast (optimized)...")
    
    # Get grid info
    grid_correction_factors = correction_data['grid_correction_factors']
    grid_coords = correction_data['grid_coords']
    
    # Handle both Dataset and DataArray cases
    if isinstance(forecast_ds, xr.DataArray):
        forecast_data = forecast_ds
        discharge_var = forecast_ds.name or 'dis24'
    else:
        var_names = list(forecast_ds.data_vars)
        discharge_var = var_names[0]
        forecast_data = forecast_ds[discharge_var]
    
    # Get forecast coordinates
    if 'lat' in forecast_data.dims:
        y_vals = forecast_data.lat.values
        x_vals = forecast_data.lon.values
    elif 'y' in forecast_data.dims:
        y_vals = forecast_data.y.values
        x_vals = forecast_data.x.values
    else:
        raise ValueError("Cannot find spatial dimensions")
    
    print(f"Forecast grid: {len(y_vals)} x {len(x_vals)} = {len(y_vals)*len(x_vals)} points")
    
    # PRE-COMPUTE NEAREST NEIGHBORS (ONCE!)
    print("Pre-computing nearest neighbors...")
    forecast_x, forecast_y = np.meshgrid(x_vals, y_vals, indexing='xy')
    forecast_coords = np.column_stack([forecast_x.ravel(), forecast_y.ravel()])
    
    from scipy.spatial import cKDTree  # Much faster than cdist
    grid_tree = cKDTree(grid_coords)
    _, nearest_indices = grid_tree.query(forecast_coords, k=1)
    nearest_indices = nearest_indices.reshape(len(y_vals), len(x_vals))
    
    # PRE-EXTRACT CORRECTION FACTORS
    obs_quantiles = grid_correction_factors['obs_quantiles']  # Shape: (n_grid, n_quantiles)
    sim_quantiles = grid_correction_factors['sim_quantiles']
    quantile_levels = grid_correction_factors['quantile_levels']
    
    print("Creating interpolation functions...")
    
    # Create corrected dataset
    corrected_ds = forecast_ds.copy(deep=True)
    
    # VECTORIZED CORRECTION FUNCTION
    def apply_correction_vectorized(forecast_2d):
        """Apply corrections to entire 2D array at once"""
        flat_forecast = forecast_2d.ravel()
        flat_nearest = nearest_indices.ravel()
        corrected_flat = np.zeros_like(flat_forecast)
        
        # Get valid (non-NaN) points
        valid_mask = ~np.isnan(flat_forecast)
        if not np.any(valid_mask):
            return forecast_2d
        
        valid_forecast = flat_forecast[valid_mask]
        valid_nearest = flat_nearest[valid_mask]
        
        # Get correction factors for all valid points at once
        obs_q_matrix = obs_quantiles[valid_nearest, :]  # Shape: (n_valid, n_quantiles)
        sim_q_matrix = sim_quantiles[valid_nearest, :]  # Shape: (n_valid, n_quantiles)
        
        # Vectorized interpolation
        corrected_valid = np.zeros_like(valid_forecast)
        
        for i in range(len(valid_forecast)):
            # Linear interpolation for each point
            try:
                corrected_valid[i] = np.interp(
                    valid_forecast[i], 
                    sim_q_matrix[i, :], 
                    obs_q_matrix[i, :],
                    left=obs_q_matrix[i, 0],
                    right=obs_q_matrix[i, -1]
                )
            except:
                corrected_valid[i] = valid_forecast[i]
        
        # Ensure no negative values
        corrected_valid = np.maximum(corrected_valid, 0.0)
        
        # Put back into full array
        corrected_flat[valid_mask] = corrected_valid
        corrected_flat[~valid_mask] = flat_forecast[~valid_mask]  # Keep NaN values
        
        return corrected_flat.reshape(forecast_2d.shape)
    
    # Process time steps
    if 'time' in forecast_data.dims:
        n_times = len(forecast_data.time)
        
        for t_idx in range(n_times):
            print(f"Processing time step {t_idx+1}/{n_times}")
            
            # Get forecast values for this time step
            forecast_2d = forecast_data.isel(time=t_idx).values
            
            # Apply vectorized correction
            corrected_2d = apply_correction_vectorized(forecast_2d)
            
            # Update the corrected dataset
            if isinstance(forecast_ds, xr.DataArray):
                corrected_ds.values[t_idx, :, :] = corrected_2d
            else:
                corrected_ds[discharge_var].values[t_idx, :, :] = corrected_2d
    else:
        # Single time slice
        print("Processing single time slice")
        forecast_2d = forecast_data.values
        corrected_2d = apply_correction_vectorized(forecast_2d)
        
        if isinstance(forecast_ds, xr.DataArray):
            corrected_ds.values = corrected_2d
        else:
            corrected_ds[discharge_var].values = corrected_2d
    
    print("✅ Bias correction complete!")
    return corrected_ds

# =============================================================================
# MAIN WORKFLOW: CREATE HISTORICAL CORRECTION FACTORS
# =============================================================================

print("=== CREATING HISTORICAL CORRECTION FACTORS ===")

file='/Users/carlos/Downloads/caudal_diario_historico_4.txt'
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
path_nc='/Users/carlos/Downloads/dis_3d_idw_optimized_1980_2018.nc'
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
        method='linear',           # Changed from 'nearest' to 'linear'
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
    forecast_path = '/Users/carlos/Documents/GitHub/24bp273032/Rst/GloFAS_2019_01_01_f.nc'
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

    corrected_ds.to_netcdf(output_path, encoding=encoding)
    
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