import numpy as np
import pandas as pd
import xarray as xr
import os
import geopandas as gpd
import scipy.stats
import shapely.geometry as sg
from scipy.spatial import cKDTree

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
    
# OPTIMIZED B-MORPH CORRECTION FUNCTION
def apply_bmorph_corrections_fast(forecast_ds, correction_data):
    
    # Get pre-computed nearest neighbors (NO KDTree building!)
    nearest_indices = correction_data['nearest_indices']
    
    # Get correction factors
    grid_correction_factors = correction_data['grid_correction_factors']
    sim_quantiles_grid = grid_correction_factors['sim_quantiles']
    obs_quantiles_grid = grid_correction_factors['obs_quantiles']
    
    # Handle dataset
    forecast_data = forecast_ds
    
    # Create corrected copy
    corrected_ds = forecast_ds.copy(deep=True)
    
    # Apply correction to the 2D slice (FAST!)
    forecast_2d = forecast_data.values
    corrected_2d = apply_bmorph_correction_to_2d_ultra_fast(
        forecast_2d, nearest_indices, sim_quantiles_grid, obs_quantiles_grid
    )
    
    # Update values
    corrected_ds.values = corrected_2d
    
    return corrected_ds

def apply_bmorph_correction_to_2d_ultra_fast(forecast_2d, nearest_indices, sim_quantiles_grid, obs_quantiles_grid):
    """
    ULTRA FAST 2D B-Morph correction - No spatial computations needed!
    Uses pre-computed nearest neighbors
    """
    flat_forecast = forecast_2d.ravel()
    flat_nearest = nearest_indices.ravel()
    corrected_flat = flat_forecast.copy()
    
    # Only correct non-zero, non-NaN pixels
    valid_mask = ~np.isnan(flat_forecast) & (flat_forecast > 0)
    
    if not np.any(valid_mask):
        return forecast_2d
    
    valid_forecast = flat_forecast[valid_mask]
    valid_nearest = flat_nearest[valid_mask]
    
    # Get correction factors for valid pixels (vectorized indexing)
    sim_q_matrix = sim_quantiles_grid[valid_nearest, :]  # Shape: (n_valid, n_quantiles)
    obs_q_matrix = obs_quantiles_grid[valid_nearest, :]  # Shape: (n_valid, n_quantiles)
    
    # VECTORIZED INTERPOLATION (can be further optimized with numba)
    corrected_valid = np.zeros_like(valid_forecast, dtype=np.float32)
    
    # Process in chunks to avoid memory issues
    chunk_size = 10000
    for i in range(0, len(valid_forecast), chunk_size):
        end_idx = min(i + chunk_size, len(valid_forecast))
        
        for j in range(i, end_idx):
            try:
                corrected_valid[j] = np.interp(
                    valid_forecast[j], 
                    sim_q_matrix[j, :], 
                    obs_q_matrix[j, :],
                    left=obs_q_matrix[j, 0],
                    right=obs_q_matrix[j, -1]
                )
            except:
                corrected_valid[j] = valid_forecast[j]
    
    # Ensure no negative values
    corrected_valid = np.maximum(corrected_valid, 0.0)

    # Clamp any extreme jumps: no more than 3× the original
    max_factor = 3.0
    orig = valid_forecast
    too_big = corrected_valid / (orig + 1e-6) > max_factor
    corrected_valid[too_big] = orig[too_big] * max_factor

    # Update only corrected pixels
    corrected_flat[valid_mask] = corrected_valid
    
    return corrected_flat.reshape(forecast_2d.shape)

def save_factors():

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

    # load forecast netcdf
    forecast_ds=xr.open_dataset('tempfile.nc')
    forecast_ds_0=forecast_ds.sel(forecast_period=forecast_ds['forecast_period'].min(),number=1)
    forecast_ds_0.rio.write_crs("EPSG:32719", inplace=True)

    #  Get bounds
    bounds = forecast_ds_0.rio.bounds()
    forecast_bbox = sg.box(*bounds)

    # Create GeoDataFrame
    forecast_gdf = gpd.GeoDataFrame(
        [1], 
        geometry=[forecast_bbox], 
        crs=forecast_ds_0.rio.crs
    )

    # Clip historical data to forecast extent
    ds = ds.rio.clip(forecast_gdf.geometry, 
                     forecast_gdf.crs, 
                     drop=True)

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
        # 'multipliers':   grid_correction_factors['multipliers'],
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

def reproject(tmpfile_interp):
    # Reproject each 2D slice separately
    print("Reprojecting to EPSG:4326...")
    import gc
    
    # Create a list to store reprojected slices
    reprojected_slices = []
    
    total_slices = len(tmpfile_interp.forecast_period) * len(tmpfile_interp.number)
    processed = 0

    for period in tmpfile_interp.forecast_period:
        period_slices = []
        
        for number in tmpfile_interp.number:
            processed += 1
            if processed % 10 == 0:
                print(f"Processing slice {processed}/{total_slices}: period={period.values}, number={number.values}")
            
            # Get 2D slice
            slice_2d = tmpfile_interp.sel(forecast_period=period, number=number)
            
            # Reproject the 2D slice
            reprojected_slice = slice_2d.rio.reproject("EPSG:4326")

            fill_value = reprojected_slice.attrs.get('_FillValue', None)

            reprojected_slice = reprojected_slice.where(reprojected_slice != fill_value, np.nan)
            reprojected_slice.attrs['_FillValue'] = np.nan

            period_slices.append(reprojected_slice)

            # Clear temporary slice from memory
            del slice_2d
            del reprojected_slice

            
            # Force garbage collection every 5 slices
            if processed % 5 == 0:
                gc.collect()
        
        # Combine slices for this forecast period
        period_combined = xr.concat(period_slices, dim='number')
        period_combined = period_combined.assign_coords(number=tmpfile_interp.number.values)
        reprojected_slices.append(period_combined)
        
        # Clear period slices from memory
        del period_slices, period_combined
        gc.collect()

    # Combine all periods
    print("Combining all reprojected periods...")
    tmpfile_interp_reprojected = xr.concat(reprojected_slices, dim='forecast_period')
    tmpfile_interp_reprojected = tmpfile_interp_reprojected.assign_coords(
        forecast_period=tmpfile_interp.forecast_period.values
    )
    
    # Clear intermediate data
    del reprojected_slices
    gc.collect()

    tmpfile_interp_reprojected['dis24'] = tmpfile_interp_reprojected['dis24'].where(
    tmpfile_interp_reprojected['dis24'] < 1e30, np.nan)
    # update the _FillValue attribute so downstream writers know
    tmpfile_interp_reprojected['dis24'].attrs['_FillValue'] = np.nan

    # rename x and y to longitude and latitude
    tmpfile_interp_reprojected = tmpfile_interp_reprojected.rename({'x': 'longitude', 'y': 'latitude'})

    print("✅ Reprojection complete!")
    return tmpfile_interp_reprojected 

def bias_correct(ds):
    
    print("=== APPLYING BIAS CORRECTION TO GLOFAS FORECAST ===")
    
    # File paths
    # forecast_path = 'tempfile.nc'
    # Load forecast data
    # print(f"Loading forecast: {forecast_path}")
    forecast_ds = ds
    print(f"✅ Forecast loaded successfully")
    print(f"   Grid size: {forecast_ds.dims}")

    # Load correction factors
    print(f"Loading correction factors:")

    npz = np.load('../output/historical_correction_factors.npz')

    # Build correction_data exactly as apply_bmorph expects
    correction_data = {
        'grid_correction_factors': {
            'sim_quantiles': npz['sim_quantiles'],  # shape (n_grid, n_q)
            'obs_quantiles': npz['obs_quantiles'],  # shape (n_grid, n_q)
        },
        'grid_coords': npz['grid_coords'],
    }
    
    # Clear npz from memory immediately
    del npz
    import gc
    gc.collect()

    # Ensure forecast has CRS information
    if not hasattr(forecast_ds, 'rio') or forecast_ds.rio.crs is None:
        print("Adding CRS information to forecast...")
        forecast_ds.rio.write_crs("EPSG:32719", inplace=True)
    
    # =========================================================================
    # PRE-COMPUTE NEAREST NEIGHBORS ONCE FOR ENTIRE FORECAST GRID
    # =========================================================================
    print("Pre-computing nearest neighbors for entire forecast grid...")
    
    forecast_data = forecast_ds
    
    # Get forecast coordinates
    y_vals = forecast_data.y.values
    x_vals = forecast_data.x.values
    
    print(f"Forecast grid: {len(y_vals)} x {len(x_vals)} = {len(y_vals)*len(x_vals)} points")
    
    # Create forecast coordinate grid ONCE
    forecast_x, forecast_y = np.meshgrid(x_vals, y_vals, indexing='xy')
    forecast_coords = np.column_stack([forecast_x.ravel(), forecast_y.ravel()])

    # Clear temporary coordinate arrays
    del forecast_x, forecast_y
    gc.collect()
    
    # Build KDTree and query ONCE
    print("Building KDTree and finding nearest neighbors...")
    grid_tree = cKDTree(correction_data['grid_coords'])
    _, nearest_indices = grid_tree.query(forecast_coords, k=1)
    nearest_indices = nearest_indices.reshape(len(y_vals), len(x_vals))
    
    # Clear temporary variables
    del forecast_coords, grid_tree
    gc.collect()

    # Add pre-computed indices to correction_data
    correction_data['nearest_indices'] = nearest_indices
    correction_data['forecast_shape'] = (len(y_vals), len(x_vals))
    
    print(f"✅ Nearest neighbors computed ONCE for entire grid!")
    
    # =========================================================================
    # APPLY CORRECTIONS USING PRE-COMPUTED NEAREST NEIGHBORS
    # =========================================================================    
    # Apply corrections
    print("\nApplying bias corrections...")
    total_slices = len(forecast_ds['forecast_period']) * len(forecast_ds['number'])
    processed = 0
    
    for time in forecast_ds['forecast_period']:
        for number in forecast_ds['number']:
            processed += 1
            if processed % 10 == 0:
                print(f"Processing slice {processed}/{total_slices}: forecast_period={time.values}, number={number.values}")
            
            ds_uncorrected = forecast_ds.sel(forecast_period=time, number=number)
            corrected_ds = apply_bmorph_corrections_fast(ds_uncorrected, correction_data)
            forecast_ds.loc[dict(forecast_period=time, number=number)] = corrected_ds

                        # Clear temporary variables from this iteration
            del ds_uncorrected, corrected_ds
            
            # Force garbage collection every 20 slices to free memory
            if processed % 20 == 0:
                gc.collect()
    
    del correction_data
    gc.collect()
    print("✅ BIAS CORRECTION COMPLETE!")

    forecast_ds = forecast_ds.to_dataset(name='dis24')

    return reproject(forecast_ds)

if __name__ == "__main__":
    corrected_forecast = bias_correct()