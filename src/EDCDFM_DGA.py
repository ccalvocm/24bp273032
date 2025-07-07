import numpy as np
import pandas as pd
import xarray as xr
import os
import geopandas as gpd
import pickle
import scipy.stats

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

def edcdfm(raw_x, raw_cdf, train_cdf, ref_cdf, train_cdf_min=1e-6):
    """
    Calculate multipliers using an adapted version of the EDCDFm technique

    This routine implements part of the PresRat bias correction method from
    Pierce et al. (2015; http://dx.doi.org/10.1175/JHM-D-14-0236.1), which is
    itself an extension of the Equidistant quantile matching (EDCDFm) technique
    of Li et al. (2010; http://dx.doi.org/10.1029/94JD00483). The part that is
    implemented here is the amended form of EDCDFm that determines
    multiplicative changes in the quantiles of a CDF.

    In particular, if the value `raw_x` falls at quantile `u_t` (in `raw_cdf`),
    then the bias-corrected value is the value in `ref_cdf` at `u_t`
    (`ref_x`) multiplied by the model-predicted change at `u_t` evaluated as
    a ratio (i.e., model future (or `raw_x`) / model historical (or
    `ref_x`)). Thus, the bias-corrected value is `raw_x` multiplied by
    `ref_x`/`train_x`. Here we only return the multiplier
    `ref_x`/`train_x`. This method preserves the model-predicted median
    (not mean) change evaluated multiplicatively. Additional corrections
    are required to preserve the mean change. Inclusion of these additional
    corrections constitutes the PresRat method.

    Parameters
    ----------
    raw_x : pandas.Series
        Series of raw values that will be used to determine the quantile `u_t`
    raw_cdf : pandas.Series
        Series of raw values that represents the CDF that is used to
        determine the non-parametric quantile of `raw_x`
    train_cdf: pandas.Series
        Series of training values that represents the CDF based on
        the same process as `raw_cdf`, but overlapping in time with `ref_cdf`
    ref_cdf: pandas.Series
         Series of ref values that represents the ref CDF and that
         overlaps in time with `train_cdf`

    Returns
    -------
    multiplier : pandas.Series
        Multipliers for `raw_x`. The pandas.Series has the same index as
        `raw_x`
    """
    # Type checking (note that the checking is more strict here then it
    # probably needs to be)
    assert isinstance(raw_x, pd.Series)
    assert isinstance(raw_cdf, pd.Series)
    assert isinstance(train_cdf, pd.Series)
    assert isinstance(ref_cdf, pd.Series)

    # Given raw_x and raw_cdf determine the quantiles u_t
    # This method is slightly more efficient than using
    # scipy.percentileofscore, especially on large arrays
    cdf_idx = np.argsort(raw_cdf)
    cdf_sort = raw_cdf.iloc[cdf_idx]
    cdf_rank = 100 * scipy.stats.rankdata(cdf_sort) / len(cdf_idx)
    offset = 100 / len(cdf_idx)
    u_t = [cdf_rank[np.searchsorted(cdf_sort, x, side='left')]-offset for x in raw_x]

    # Given u_t and train_cdf determine train_x
    train_x = np.percentile(train_cdf, u_t)
    train_x[train_x < train_cdf_min] = train_cdf_min

    # Given u_t and ref_cdf determine ref_x
    ref_x = np.percentile(ref_cdf, u_t)

    # Calculate multiplier
    multiplier = ref_x / train_x
    return pd.Series(multiplier, index=raw_x.index, name='multiplier')

def train_bmorph_correction(observed, simulated, n_quantiles=25, min_pts=2000):
    """
    Train B-Morph static quantile-based correction multipliers using edcdfm.
    """
    # Clean data
    mask = np.isfinite(observed) & np.isfinite(simulated)
    obs_clean = observed[mask]
    sim_clean = simulated[mask]
    
    if len(obs_clean) < min_pts:
        return None
    
    try:
        # Convert to pandas Series for edcdfm
        obs_series = pd.Series(obs_clean)
        sim_series = pd.Series(sim_clean)
        
        # Compute empirical CDFs at quantile levels
        quantile_levels = np.linspace(0.01, 0.99, n_quantiles)
        
        # Get quantile values for both series
        obs_quantiles = obs_series.quantile(quantile_levels)
        sim_quantiles = sim_series.quantile(quantile_levels)
        
        # Create CDF series with INTEGER INDEX and quantile values as data
        # edcdfm expects CDF to be indexed by positions 0, 1, 2, ... not by values
        obs_cdf = pd.Series(obs_quantiles.values, index=range(n_quantiles))
        sim_cdf = pd.Series(sim_quantiles.values, index=range(n_quantiles))
        
        # raw_x should also be indexed by integers
        raw_x = pd.Series(sim_quantiles.values, index=range(n_quantiles))
        
        # Now call edcdfm with proper CDF inputs
        multipliers = edcdfm(
            raw_x=raw_x,              # quantile values with integer index
            raw_cdf=sim_cdf,          # simulation CDF with integer index
            train_cdf=sim_cdf,        # training CDF (same as raw for bias correction)
            ref_cdf=obs_cdf,          # reference (observed) CDF with integer index
            train_cdf_min=1e-6
        )
        
        return {
            'multipliers': np.asarray(multipliers),
            'quantile_levels': quantile_levels,
            'obs_quantiles': np.asarray(obs_quantiles),
            'sim_quantiles': np.asarray(sim_quantiles),
            'obs_mean': float(obs_clean.mean()),
            'sim_mean': float(sim_clean.mean()),
            'n_training': len(obs_clean),
            'n_quantiles': n_quantiles,
            'method': 'bmorph_edcdfm'
        }
        
    except Exception as e:
        print(f"edcdfm failed: {e}")
        return None

def interpolate_bmorph_factors_to_grid(station_coords, station_corrections, grid_coords, method='rbf', smooth_sigma=2.0):
    """
    Interpolate edcdfm-based correction factors to full model grid
    """
    from scipy.ndimage import gaussian_filter
    from scipy.interpolate import Rbf, griddata
    import numpy as np

    valid_stations = list(station_corrections.keys())
    station_xy = np.array([station_coords[s] for s in valid_stations])

    first = station_corrections[valid_stations[0]]
    
    # Handle edcdfm-based corrections
    quantile_levels = first['quantile_levels']
    n_quantiles = len(quantile_levels)
    
    # Stack multipliers and quantiles from all stations
    multipliers_mat = np.stack([station_corrections[s]['multipliers'] for s in valid_stations])
    sim_quantiles_mat = np.stack([station_corrections[s]['sim_quantiles'] for s in valid_stations])
    obs_quantiles_mat = np.stack([station_corrections[s]['obs_quantiles'] for s in valid_stations])

    npts = len(grid_coords)
    grid_multipliers = np.full((npts, n_quantiles), np.nan)
    grid_sim_quantiles = np.full((npts, n_quantiles), np.nan)
    grid_obs_quantiles = np.full((npts, n_quantiles), np.nan)

    for i in range(n_quantiles):
        mult_vals = multipliers_mat[:, i]
        sim_vals = sim_quantiles_mat[:, i]
        obs_vals = obs_quantiles_mat[:, i]
        
        # exclude stations with NaN or inf
        ok = np.isfinite(mult_vals) & np.isfinite(sim_vals) & np.isfinite(obs_vals)
        if ok.sum() < 3:
            # not enough points for RBF: use nearest
            grid_multipliers[:, i] = griddata(station_xy[ok], mult_vals[ok], grid_coords, method='nearest')
            grid_sim_quantiles[:, i] = griddata(station_xy[ok], sim_vals[ok], grid_coords, method='nearest')
            grid_obs_quantiles[:, i] = griddata(station_xy[ok], obs_vals[ok], grid_coords, method='nearest')
            continue

        if method == 'rbf':
            rbf_mult = Rbf(station_xy[ok,0], station_xy[ok,1], mult_vals[ok],
                           function='multiquadric', smooth=0.1)
            rbf_sim = Rbf(station_xy[ok,0], station_xy[ok,1], sim_vals[ok],
                          function='multiquadric', smooth=0.1)
            rbf_obs = Rbf(station_xy[ok,0], station_xy[ok,1], obs_vals[ok],
                          function='multiquadric', smooth=0.1)
            grid_multipliers[:,i] = rbf_mult(grid_coords[:,0], grid_coords[:,1])
            grid_sim_quantiles[:,i] = rbf_sim(grid_coords[:,0], grid_coords[:,1])
            grid_obs_quantiles[:,i] = rbf_obs(grid_coords[:,0], grid_coords[:,1])
        else:
            grid_multipliers[:,i] = griddata(station_xy[ok], mult_vals[ok], grid_coords,
                                           method='nearest', fill_value=np.nan)
            grid_sim_quantiles[:,i] = griddata(station_xy[ok], sim_vals[ok], grid_coords,
                                             method='nearest', fill_value=np.nan)
            grid_obs_quantiles[:,i] = griddata(station_xy[ok], obs_vals[ok], grid_coords,
                                             method='nearest', fill_value=np.nan)

    # interpolate means
    obs_means = np.array([station_corrections[s].get('obs_mean') for s in valid_stations])
    sim_means = np.array([station_corrections[s].get('sim_mean') for s in valid_stations])
    okm = np.isfinite(obs_means) & np.isfinite(sim_means)
    
    if method == 'rbf' and okm.sum() >= 3:
        rbf_om = Rbf(station_xy[okm,0], station_xy[okm,1], obs_means[okm],
                      function='multiquadric', smooth=0.1)
        rbf_sm = Rbf(station_xy[okm,0], station_xy[okm,1], sim_means[okm],
                      function='multiquadric', smooth=0.1)
        grid_obs_mean = rbf_om(grid_coords[:,0], grid_coords[:,1])
        grid_sim_mean = rbf_sm(grid_coords[:,0], grid_coords[:,1])
    else:
        grid_obs_mean = griddata(station_xy[okm], obs_means[okm], grid_coords, method='nearest')
        grid_sim_mean = griddata(station_xy[okm], sim_means[okm], grid_coords, method='nearest')

    # spatial smoothing if regular grid
    if smooth_sigma > 0:
        ux, uy = np.unique(grid_coords[:,0]), np.unique(grid_coords[:,1])
        if len(grid_coords) == ux.size * uy.size:
            shape = (uy.size, ux.size)
            for i in range(n_quantiles):
                grid_multipliers[:,i] = gaussian_filter(grid_multipliers[:,i].reshape(shape), smooth_sigma).ravel()
                grid_sim_quantiles[:,i] = gaussian_filter(grid_sim_quantiles[:,i].reshape(shape), smooth_sigma).ravel()
                grid_obs_quantiles[:,i] = gaussian_filter(grid_obs_quantiles[:,i].reshape(shape), smooth_sigma).ravel()

    return {
        'quantile_levels': quantile_levels,
        'multipliers': grid_multipliers,
        'sim_quantiles': grid_sim_quantiles,
        'obs_quantiles': grid_obs_quantiles,
        'obs_means': grid_obs_mean,
        'sim_means': grid_sim_mean,
        'n_quantiles': n_quantiles
    }

def apply_bmorph_correction(forecast_values, correction_factors):
    """
    Apply B-Morph correction to forecast values using edcdfm-based correction factors
    
    Parameters:
    -----------
    forecast_values : array-like or scalar
        Raw forecast values to be corrected
    correction_factors : dict
        Contains 'multipliers', 'quantile_levels', 'obs_quantiles', 'sim_quantiles'
        
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
    
    # Get the multipliers and quantile info
    multipliers = correction_factors['multipliers']
    sim_quantiles = correction_factors['sim_quantiles']
    
    # Apply multipliers based on quantile position
    try:
        # Find which quantile bin each forecast value falls into
        # and apply the corresponding multiplier
        corrected_clean = np.interp(
            clean_forecast, 
            sim_quantiles,  # x-coordinates (simulation quantiles)
            multipliers,    # y-coordinates (multipliers)
            left=multipliers[0],   # extrapolate with first multiplier
            right=multipliers[-1]  # extrapolate with last multiplier
        ) * clean_forecast  # multiply by the forecast value
        
        # Ensure no negative values for streamflow
        corrected_clean = np.maximum(corrected_clean, 0.0)
        
    except Exception as e:
        print(f"Warning: B-Morph correction failed, returning original values. Error: {e}")
        corrected_clean = clean_forecast
    
    # Put corrected values back
    corrected_forecast = forecast_values.copy()
    corrected_forecast[mask] = corrected_clean
    
    return corrected_forecast[0] if return_scalar else corrected_forecast
def apply_bmorph_corrections_to_forecast(forecast_ds, correction_data):
    """
    Apply B-Morph corrections to forecast using edcdfm-based grid correction factors
    """
    print("Applying B-Morph bias corrections to forecast...")
    
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
    
    # PRE-COMPUTE NEAREST NEIGHBORS
    print("Pre-computing nearest neighbors...")
    forecast_x, forecast_y = np.meshgrid(x_vals, y_vals, indexing='xy')
    forecast_coords = np.column_stack([forecast_x.ravel(), forecast_y.ravel()])
    
    from scipy.spatial import cKDTree
    grid_tree = cKDTree(grid_coords)
    _, nearest_indices = grid_tree.query(forecast_coords, k=1)
    nearest_indices = nearest_indices.reshape(len(y_vals), len(x_vals))
    
    # PRE-EXTRACT B-MORPH CORRECTION FACTORS (edcdfm-based)
    # Fix: Use the correct field names from your data structure
    multipliers_grid = grid_correction_factors['multipliers']      # shape: (n_grid_points, n_quantiles)
    sim_quantiles_grid = grid_correction_factors['sim_quantiles']  # shape: (n_grid_points, n_quantiles)
    obs_quantiles_grid = grid_correction_factors['obs_quantiles']  # shape: (n_grid_points, n_quantiles)
    
    # Create corrected dataset
    corrected_ds = forecast_ds.copy(deep=True)
    
    # B-MORPH CORRECTION FUNCTION
    def apply_bmorph_vectorized(forecast_2d):
        """Apply B-Morph corrections to entire 2D array using quantile mapping"""
        flat_forecast = forecast_2d.ravel()
        flat_nearest = nearest_indices.ravel()
        corrected_flat = flat_forecast.copy()
        
        # Only correct non-zero, non-NaN pixels
        correction_mask = ~np.isnan(flat_forecast) & (flat_forecast > 0)
        
        if not np.any(correction_mask):
            return forecast_2d
        
        valid_forecast = flat_forecast[correction_mask]
        valid_nearest = flat_nearest[correction_mask]
        
        # Get B-Morph correction factors for valid pixels
        sim_q_matrix = sim_quantiles_grid[valid_nearest, :]         # shape: (n_valid, n_quantiles)
        obs_q_matrix = obs_quantiles_grid[valid_nearest, :]         # shape: (n_valid, n_quantiles)
        
        # Apply quantile mapping correction (NOT multipliers)
        corrected_valid = np.zeros_like(valid_forecast)
        
        for i in range(len(valid_forecast)):
            try:
                # Direct quantile mapping: interpolate from sim quantiles to obs quantiles
                corrected_valid[i] = np.interp(
                    valid_forecast[i], 
                    sim_q_matrix[i, :],      # simulation quantiles for this grid point
                    obs_q_matrix[i, :],      # corresponding observed quantiles
                    left=obs_q_matrix[i, 0],     # extrapolate with first obs quantile
                    right=obs_q_matrix[i, -1]    # extrapolate with last obs quantile
                )
            except:
                corrected_valid[i] = valid_forecast[i]
        
        # Ensure no negative values
        corrected_valid = np.maximum(corrected_valid, 0.0)
        
        # Update only the corrected pixels
        corrected_flat[correction_mask] = corrected_valid
        
        return corrected_flat.reshape(forecast_2d.shape)
    
    # Process time steps
    if 'time' in forecast_data.dims:
        n_times = len(forecast_data.time)
        
        for t_idx in range(n_times):
            print(f"Processing time step {t_idx+1}/{n_times}")
            
            forecast_2d = forecast_data.isel(time=t_idx).values
            corrected_2d = apply_bmorph_vectorized(forecast_2d)
            
            if isinstance(forecast_ds, xr.DataArray):
                corrected_ds.values[t_idx, :, :] = corrected_2d
            else:
                corrected_ds[discharge_var].values[t_idx, :, :] = corrected_2d
    else:
        # Single time slice
        print("Processing single time slice")
        forecast_2d = forecast_data.values
        corrected_2d = apply_bmorph_vectorized(forecast_2d)
        
        if isinstance(forecast_ds, xr.DataArray):
            corrected_ds.values = corrected_2d
        else:
            corrected_ds[discharge_var].values = corrected_2d
    
    print("✅ B-Morph bias correction complete!")
    return corrected_ds

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

# =============================================================================
# TRAIN B-MORPH CORRECTION FACTORS AT EACH STATION
# =============================================================================
print("\n=== TRAINING B-MORPH CORRECTION FACTORS ===")

correction_factors_stations = {}
nse_results = {}
station_coords_dict = {}

ds_common = ds_nat.loc[common_dates, common_stations]
df_common = df_pivot.loc[common_dates, common_stations]

for station in common_stations:

    observed = df_common[station].values
    simulated = ds_common[station].values
    
    # Calculate NSE for evaluation
    nse = nash_sutcliffe_efficiency(observed, simulated)
    nse_results[station] = nse
    
    # Train B-Morph correction factors
    correction_model = train_bmorph_correction(observed, simulated, 
                                               n_quantiles=25)
    
    if correction_model is not None:
        correction_factors_stations[station] = correction_model
        
        # Store station coordinates
        station_coords_dict[station] = (gdf_nat.loc[station, 'geometry'].x, 
                                      gdf_nat.loc[station, 'geometry'].y)
        
        print(f"✅ {station}: NSE = {nse:.3f}, B-Morph factors trained")
    else:
        print(f"❌ {station}: NSE = {nse:.3f}, Insufficient data for B-Morph")

# =============================================================================
# INTERPOLATE B-MORPH FACTORS TO FULL GRID
# =============================================================================
print("\n=== INTERPOLATING B-MORPH FACTORS TO FULL GRID ===")

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

grid_correction_factors = interpolate_bmorph_factors_to_grid(
    station_coords_dict, 
    correction_factors_stations, 
    grid_coords, 
    method='rbf',
    smooth_sigma=2.0
)

# =============================================================================
# SAVE HISTORICAL CORRECTION FACTORS
# =============================================================================
print("\n=== SAVING CORRECTION FACTORS ===")

correction_data = {
    'multipliers':   grid_correction_factors['multipliers'],
    'sim_quantiles': grid_correction_factors['sim_quantiles'],
    'grid_coords':   grid_coords,
    'obs_quantiles': grid_correction_factors['obs_quantiles'],
}

np.savez_compressed(
    '../output/historical_correction_factors.npz',
    multipliers=grid_correction_factors['multipliers'].astype('float32') ,
    sim_quantiles=grid_correction_factors['sim_quantiles'].astype('float32') ,
    grid_coords=grid_coords.astype('float32') ,
    obs_quantiles=grid_correction_factors['obs_quantiles'].astype('float32') ,
)

# Save NSE results
nse_df = pd.DataFrame(list(nse_results.items()), columns=['Station', 'NSE'])
nse_df.to_csv(os.path.join('..', 'output', 'nse_results.csv'), index=False)

# Save summary
summary_stats = {
    'Training Period': f"{ds_nat.index.min().date()} to {ds_nat.index.max().date()}",
    'Stations with Corrections': len(correction_factors_stations),
    'Total Stations Evaluated': len(common_stations),
    'Grid Points': len(grid_coords),
    'Mean NSE': np.nanmean(list(nse_results.values())),
    'Stations with NSE > 0.5': np.sum([nse>0.5 for nse in nse_results.values() if not np.isnan(nse)])
}
summary_df = pd.DataFrame(list(summary_stats.items()), columns=['Metric', 'Value'])
summary_df.to_csv('../output/correction_factors_summary.csv', index=False)

print("\n=== HISTORICAL CORRECTION FACTORS CREATION COMPLETE ===")

def main():
    """Main function to apply corrections to GloFAS forecast"""
    
    print("=== APPLYING BIAS CORRECTION TO GLOFAS FORECAST ===")
    
    # File paths
    forecast_path = os.path.join('..','Rst','GloFAS_2019_01_01_f.nc')
    
    # Load forecast data
    print(f"Loading forecast: {forecast_path}")
    forecast_ds = xr.open_dataset(forecast_path)
    print(f"✅ Forecast loaded successfully")
    print(f"   Variables: {list(forecast_ds.data_vars)}")
    print(f"   Grid size: {forecast_ds.dims}")

    # Load correction factors
    print(f"Loading correction factors:")

    npz = np.load('../output/historical_correction_factors.npz')

    # Build correction_data exactly as apply_bmorph expects
    correction_data = {
        'grid_correction_factors': {
            'multipliers':   npz['multipliers'],    # shape (n_grid, n_q)
            'sim_quantiles': npz['sim_quantiles'],  # shape (n_grid, n_q)
            'obs_quantiles': npz['obs_quantiles'],  # shape (n_grid, n_q)
        },
        'grid_coords': npz['grid_coords'],
    }

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
    corrected_ds = apply_bmorph_corrections_to_forecast(forecast_ds, correction_data)
    
    # Save corrected forecast
    output_path = forecast_path.replace('.nc', '_bias_corrected_3.nc')
    print(f"\nSaving corrected forecast: {output_path}")
    
    # Add metadata
    corrected_ds.attrs['bias_correction_applied'] = 'True'  # String instead of boolean

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