# install.packages("devtools")
# install.packages("lwgeom")
# install.packages("ncdf4")
devtools::install_github("mkkallio/hydrostreamer", dependencies = TRUE)

library(hydrostreamer)
library(raster)
library(lubridate)
library(dplyr)
library(sf)
library(ncdf4) 
library(lwgeom)    # for st_make_valid()

setwd('/Users/carlos/Documents/GitHub/24bp273032/src/')
# basins <- st_read("../geodata/basinQ.gpkg")
basins <- st_read("../geodata/basins_with_elevation.gpkg")

dem <- brick("../Rst/dem90.tif", package = "hydrostreamer")

nc <- nc_open("../Rst/GloFAS_2025_06_13_f.nc")

lon <- ncvar_get(nc, "longitude")  # or "lon" depending on your file
lat <- ncvar_get(nc, "latitude")   # or "lat" depending on your file
# Get the extent
lon_min <- min(lon)
lon_max <- max(lon)
lat_min <- min(lat)
lat_max <- max(lat)

dis24_array <- ncvar_get(nc, "dis24", 
                         start = c(1, 1, 1, 1, 1), 
                         count = c(-1, -1, 1, 1, 1))

# Transpose and flip if needed
r <- raster(t(dis24_array))
# Close NetCDF
nc_close(nc)

# Set the correct extent
extent(r) <- extent(lon_min, lon_max, lat_min, lat_max)
crs(r) <- "+proj=longlat +datum=WGS84 +no_defs"

# Shift longitude if needed (GloFAS 0-360 to -180 to 180)
if (lon_min > 180) {
  extent(r) <- extent(extent(r)@xmin - 360, extent(r)@xmax - 360,
                      extent(r)@ymin, extent(r)@ymax)
}
# Now reproject to EPSG:32719
runoff_32719 <- projectRaster(
  r,
  crs    = st_crs(basins)$proj4string,
  method = "bilinear"
)

#plot(runoff_32719)
# plot(basins)
# plot(st_union(basins), add=TRUE)
# plot(river, add=TRUE)

# DASYMETRIC MAPPING WITH LINES
# Try method 1: force simple geometries
tryCatch({
  basins <- st_cast(basins, "MULTIPOLYGON")
  basins <- st_buffer(basins, 0)
  basins <- sf::st_make_valid(basins)
}, error=function(e) message("Method 1 failed:", e$message))

# Method 2: Remove problematic basins one by one
valid_basins <- st_is_valid(basins, reason=TRUE)
problem_basins <- which(valid_basins != "Valid Geometry")
if(length(problem_basins) > 0) {
  message(paste0("Removing ", length(problem_basins), " problematic basins."))
  basins <- basins[-problem_basins,]
}

# Manually bypass st_union which is causing problems
message("Creating simplified AOI...")
aoi_geom <- st_as_sfc(st_bbox(basins), crs = st_crs(basins))
aoi      <- st_sf(geometry = aoi_geom)
# Use bounding box instead of union

# 2) Read the forecast reference time (global attribute)
ref_secs   <- ncatt_get(nc, 0, "forecast_reference_time")$value
start_date <- as.Date(as.POSIXct(ref_secs, origin="1970-01-01", tz="UTC"))
# ...then use hs_unit as before...
source_runoff <- raster_to_HS(
  runoff_32719, 
  unit       = "m3/s",
  date       = ymd(start_date),
  timestep   = "day",
  aoi        = aoi,
  names      = "LORA"
)

# Descriptive statistics
summary(source_runoff$runoff_ts)

# If source_runoff is a hydrostreamer object with a raster slot:
if (!is.null(source_runoff$raster)) {
  plot(source_runoff$raster, main = "source_runoff (m3/s)")
}

###ww
river <- st_read("../geodata/riverQ.gpkg")

# basins$elevation <- sapply(raster::extract(dem, basins), mean)

river2 <- dplyr::left_join(
  river,
  as.data.frame(basins)[, c("DN", "elevation")],
  by = "DN"
)
river2 <- river2[!is.na(river2$elevation), ]
river2$river_uid <- seq_len(nrow(river2))

library(sf)
river2_buffered <- st_buffer(river2, dist = 500)  # 500 meters buffer
river2 <- st_transform(river2, crs(runoff_32719))

# Use buffered lines to get better spatial coverage
river2_buffered <- st_buffer(river2, dist = 200)  # Smaller buffer
river2$discharge_m3s <- raster::extract(runoff_32719, river2_buffered, fun = mean, na.rm = TRUE)

summary(river2$discharge_m3s)
# Check the results
sum(river2$discharge_m3s)
sum(values(runoff_32719), na.rm = TRUE)

plot(river2["discharge_m3s"], main = "River Discharge (m3/s)")

A2LDM <- interpolate_runoff(
  source_runoff,
  river2_buffered,
  riverID = "river_uid",
  intensive = FALSE
)

A2LDM$mean_runoff <- sapply(A2LDM$runoff_ts, function(x) mean(x$LORA))
plot(A2LDM["mean_runoff"])

# Check the results
sum(A2LDM$mean_runoff)
sum(values(runoff_32719), na.rm = TRUE)

summary(A2LDM$mean_runoff)
summary(values(runoff_32719))


A2LDM <- interpolate_runoff(source_runoff, river2,
                            dasymetric = "elevation",
                          riverID = "river_uid")
A2LDM$mean_runoff <- sapply(A2LDM$runoff_ts, function(x) mean(x$LORA))
plot(A2LDM[,"mean_runoff"])
# Interpolate runoff using river network
library(raster)
library(sf)

# Create a template raster with desired extent and resolution
template_raster <- raster(extent(A2LDM), 
                         res = 100, # Set your desired resolution (meters)
                         crs = st_crs(A2LDM)$proj4string)

# Rasterize the mean_runoff values from the river lines
mean_runoff_raster <- rasterize(
  A2LDM, 
  template_raster, 
  field = "mean_runoff", 
  fun = mean,
  background = NA
)

# Export as GeoTIFF
writeRaster(mean_runoff_raster, 
           filename = "../output/mean_runoff2.tif", 
           format = "GTiff", 
           overwrite = TRUE)


##other test
high_runoff <- runoff_32719 > 0.1
river_high <- river2[!is.na(raster::extract(high_runoff, river2)) & 
                     raster::extract(high_runoff, river2) == 1, ]

# Get runoff values for rivers in high-runoff areas
river_high$runoff_values <- raster::extract(runoff_32719, river_high)
summary(river_high$runoff_values)

#compare
cat("High-runoff rivers mean:", mean(river_high$runoff_values, na.rm = TRUE), "\n")
cat("All rivers mean:", mean(raster::extract(runoff_32719, river2), na.rm = TRUE), "\n")

# Plot to see where your important rivers are
plot(runoff_32719, main = "Rivers in High-Runoff Areas")
plot(st_geometry(river2), add = TRUE, col = "lightgray", lwd = 0.5)  # All rivers
plot(st_geometry(river_high), add = TRUE, col = "red", lwd = 2)      # High-runoff rivers
legend("topright", c("All rivers", "High-runoff rivers"), 
       col = c("lightgray", "red"), lwd = c(0.5, 2))

# Method 1: Use accumulate_runoff() if available
if (exists("accumulate_runoff", mode = "function")) {
  A2LDM_routed <- accumulate_runoff(A2LDM, riverID = "river_uid")
} else {
  # Method 2: Use route_runoff() or similar function
  A2LDM_routed <- route_runoff(A2LDM, riverID = "river_uid")
}

# Check the routed values
summary(A2LDM_routed$discharge)  # or whatever the output column is named