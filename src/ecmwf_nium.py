from ecmwf.opendata import Client
import xarray as xr
from datetime import date, timedelta

days_back = 20
target_date = date.today() - timedelta(days=days_back)
day = target_date.strftime("%Y%m%d")
run = "00"  # or "12", "06", "18" depending on what's available

client = Client(source="ecmwf")

request = {
"date": -3,
"time": 12,
"type": "fc",
"step": [3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60],
"param": [ "tp"],
}

ds=client.retrieve(request, "surface1.grib2")
xr_ds = xr.open_dataset("surface1.grib2", 
                        engine="cfgrib")
xr_ds.tp.isel(step=0).plot()