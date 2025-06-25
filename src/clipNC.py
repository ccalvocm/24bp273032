import xarray as xr

# Input / output
src = "/Users/carlos/Downloads/dis_1980_2018.nc"
dst = "/Users/carlos/Downloads/dis_1980_2018_clipped_region.nc"
# Your bbox
min_lon, min_lat, max_lon, max_lat = -71.71782, -32.28247, -69.809361, -29.0366

# Open with dask-chunks so you don't blow memory
ds = xr.open_dataset(src, chunks={"time":1, "lat":500, "lon":500})

# Ensure lat is monotonic increasing
if float(ds.lat[0]) > float(ds.lat[-1]):
    ds = ds.sortby("lat")

# Slice out the region
sub = ds.sel(
    lon=slice(min_lon, max_lon),
    lat=slice(min_lat, max_lat),
)

# Optional: set compression on all variables
enc = {v:{"zlib":True, "complevel":5} for v in sub.data_vars}

# Write to disk
# … after your sel(…) call …
# materialize in memory/dask cache
sub = sub.persist()

# now write
sub.to_netcdf(
    dst,
    encoding=enc,
    engine="h5netcdf",
    format="NETCDF4",
)
print("Written clipped file to:", dst)