# install.packages("devtools")
# install.packages("lwgeom")
# install.packages("ncdf4")
# devtools::install_github("mkkallio/hydrostreamer")

library(hydrostreamer)
library(raster)
library(lubridate)
library(dplyr)
library(sf)
library(ncdf4) 
library(lwgeom)    # for st_make_valid()

setwd('/Users/carlos/Documents/GitHub/24bp273032/src/')
basins <- st_read("../geodata/basinQ.gpkg")

dem <- brick("../Rst/dem90.tif", package = "hydrostreamer")

runoff <- brick("../Rst/GloFAS_2025_06_13_f.nc", varname="dis24")

dis24_array <- ncvar_get(nc, "dis24", 
                         start = c(1, 1, 1, 1, 1), 
                         count = c(-1, -1, 1, 1, 1))

# Transpose and flip if needed
r <- raster(t(dis24_array))
extent(r) <- extent(288.3, 290.2, -32.3, -29.05)  # adjust as needed
crs(r)    <- "+proj=longlat +datum=WGS84 +no_defs"

# Shift longitude if needed (GloFAS 0-360 to -72 to -70 for Chile)
extent(r) <- extent(extent(r)@xmin - 360, extent(r)@xmax - 360,
                    extent(r)@ymin, extent(r)@ymax)

# Now reproject to EPSG:32719
runoff_32719 <- projectRaster(
  r,
  crs    = st_crs(basins)$proj4string,
  method = "bilinear"
)

plot(runoff_32719)
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

###ww
river <- st_read("../geodata/riverQ.gpkg")

elevation_values <- raster::extract(dem, basins)
basins$elevation <- sapply(raster::extract(dem, basins), mean)

river2 <- dplyr::left_join(
  river,
  as.data.frame(basins)[, c("DN", "elevation")],
  by = "DN"
)

A2LDM <- interpolate_runoff(source_runoff, river2,
                            dasymetric = "elevation",
                          riverID = "DN")

A2LDM$mean_runoff <- sapply(A2LDM$runoff_ts, function(x) mean(x$LORA))

library(RColorBrewer)

# Define breaks and color palette
breaks <- quantile(A2LDM$mean_runoff, probs = seq(0, 1, length.out = 8), na.rm = TRUE)
cols <- brewer.pal(7, "YlGnBu")

plot(
  A2LDM["mean_runoff"],
  breaks = breaks,
  col = cols,
  main = "Mean Runoff (m3/s)",
  key.pos = 1
)
# Interpolate runoff using river network
RB <- interpolate_runoff(source_runoff,
                         river2,
                         basins = basins,
                         riverID = "DN")
plot(RB)

