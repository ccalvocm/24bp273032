import pandas as pd
import geopandas as gpd
import os
import xarray as xr
import rioxarray as rxr

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
print("Data and metadata saved successfully.")

path_glofas='/Users/carlos/Documents/GitHub/24bp273032/Rst/dis_1980_2018_clip.nc'
glofas_hist=xr.open_dataset(path_glofas, chunks={'time': 1, 'lat': 500, 'lon': 500})
glofas_hist.rio.write_crs("EPSG:4326", inplace=True)
glofas_h_utm=glofas_hist.rio.reproject("EPSG:32719")

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

print("Sampling complete!")
ds_nat.to_csv(os.path.join('..',
                          'output',
                          'dis_3d_idw_optimized_1980_2018.csv'),
                          encoding='utf-8',
                          index_label='FECHA',
                          date_format='%Y-%m-%dT%H:%M:%S.%fZ',
                          float_format='%.2f')
print("Data sampled and saved successfully.")

# FIX: Align timezone awareness
import numpy as np

# Convert ds_nat index to timezone-aware (UTC) to match df_pivot
ds_nat.index = pd.to_datetime(ds_nat.index).tz_localize('UTC')

# OR alternatively, remove timezone from df_pivot:
# df_pivot.index = df_pivot.index.tz_localize(None)

print(f"ds_nat dates: {ds_nat.index.min()} to {ds_nat.index.max()}")
print(f"df_pivot dates: {df_pivot.index.min()} to {df_pivot.index.max()}")

# Find common dates
common_dates = ds_nat.index.intersection(df_pivot.index)
print(f"Common dates: {len(common_dates)}")

# Find common stations (columns)
common_stations = ds_nat.columns.intersection(df_pivot.columns)
print(f"Common stations: {len(common_stations)}")

if len(common_dates) > 0 and len(common_stations) > 0:
    # Subset to common dates and stations
    ds_common = ds_nat.loc[common_dates, common_stations]
    df_common = df_pivot.loc[common_dates, common_stations]
    
    # Calculate NSE for each station
    nse_results = {}
    for station in common_stations:
        observed = df_common[station].values
        simulated = ds_common[station].values
        
        nse = nash_sutcliffe_efficiency(observed, simulated)
        nse_results[station] = nse
        
        # Count valid pairs
        valid_pairs = np.sum(~(np.isnan(observed) | np.isnan(simulated)))
        print(f"{station}: NSE = {nse:.3f} (n={valid_pairs})")
    
    # Overall statistics
    valid_nse = [nse for nse in nse_results.values() if not np.isnan(nse)]
    if valid_nse:
        print(f"\nOverall NSE statistics:")
        print(f"Mean NSE: {np.mean(valid_nse):.3f}")
        print(f"Median NSE: {np.median(valid_nse):.3f}")
        print(f"Min NSE: {np.min(valid_nse):.3f}")
        print(f"Max NSE: {np.max(valid_nse):.3f}")
        print(f"Stations with NSE > 0.5: {np.sum(np.array(valid_nse) > 0.5)}/{len(valid_nse)}")
    
    # Save NSE results
    nse_df = pd.DataFrame(list(nse_results.items()), columns=['Station', 'NSE'])
    nse_df.to_csv(os.path.join('..', 'output', 'nse_results.csv'), index=False)
    
else:
    print("No common dates or stations found!")
    print("Check date ranges overlap:")
    print(f"ds_nat: {ds_nat.index.min()} to {ds_nat.index.max()}")
    print(f"df_pivot: {df_pivot.index.min()} to {df_pivot.index.max()}")

print("NSE analysis complete!")