import pandas as pd
import geopandas as gpd
import os
import xarray as xr
import rioxarray as rxr
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
# glofas_h_utm.sel(time='2018-12-31').dis.plot()
# sample glofas_h_utm at points gdf_nat over time
