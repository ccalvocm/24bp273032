from ecmwf.opendata import Client
import xarray as xr
from datetime import date, timedelta

target_date = date.today()
day = target_date.strftime("%Y%m%d")
step=[12,24]
client = Client(source="ecmwf")

request = {
"date": day,
"time": 0,
"type": "fc",
"step": step,
"param": ["tp"],
}

ds=client.retrieve(request, "surface1.grib2")
xr_ds = xr.open_dataset("surface1.grib2", 
                        engine="cfgrib")
xr_ds.tp.isel(step=0).plot()