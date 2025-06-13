install.packages("devtools")
install.packages("lwgeom")
devtools::install_github("mkkallio/hydrostreamer")

library(hydrostreamer)
library(raster)
library(lubridate)
library(dplyr)
library(sf)

library(sf)
library(lwgeom)    # for st_make_valid()

data(example_basins)
data(example_rivers)

runoff <- brick(system.file("extdata", "runoff.tif", package = "hydrostreamer"))
dem <- brick(system.file("extdata", "dem.tif", package = "hydrostreamer")) 
# plot(runoff[[3]]) 
# plot(basins)
# plot(st_union(basins), add=TRUE)
# plot(river, add=TRUE)
# plot(dem, add=TRUE)


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
aoi <- st_as_sfc(st_bbox(basins))

runoff <- brick(system.file("extdata","runoff.tif", package="hydrostreamer"))
dem    <- brick(system.file("extdata","dem.tif",    package="hydrostreamer"))

# Use bounding box instead of union
source_runoff <- raster_to_HS(runoff, 
                             unit = "mm/s",
                             date = ymd("1980-01-01"), 
                             timestep = "month", 
                             aoi = aoi,  # Using bounding box as AOI
                             names = "LORA")

elevation_values <- raster::extract(dem, basins)
basins$elevation <- sapply(elevation_values, mean)
river <- dplyr::filter(river, SEGMENT_ID %in% basins$SEGMENT_ID) %>% 
  dplyr::mutate(elevation = basins$elevation)
RB <- interpolate_runoff(source_runoff,
                         river,
                         basins = basins,
                         riverID = "SEGMENT_ID")

plot(RB)

