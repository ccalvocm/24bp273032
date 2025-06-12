install.packages("devtools")
devtools::install_github("mkkallio/hydrostreamer")

library(hydrostreamer)
library(raster)
library(lubridate)
library(dplyr)
library(sf)

data(example_rivers)
data(example_basins)
runoff <- brick(system.file("extdata", "runoff.tif", package = "hydrostreamer"))
dem <- brick(system.file("extdata", "dem.tif", package = "hydrostreamer")) 
plot(runoff[[1]]) 
plot(basins)
plot(st_union(basins), add=TRUE)
plot(river, add=TRUE)
source_runoff <- raster_to_HS(runoff, 
                              unit = "mm/s",
                              date = ymd("1980-01-01"), 
                              timestep = "month", 
                              aoi = st_union(basins),
                              names = "LORA")

elevation_values <- raster::extract(dem, basins)
basins$elevation <- sapply(elevation_values, mean)
river <- dplyr::filter(river, SEGMENT_ID %in% basins$SEGMENT_ID) %>% 
  dplyr::mutate(elevation = basins$elevation)
A2LDM <- interpolate_runoff(source_runoff, river,
                            dasymetric = "elevation",
                            riverID = "SEGMENT_ID")
#> although coordinates are longitude/latitude, st_intersects assumes that they are planar

A2LDM$mean_runoff <- sapply(A2LDM$runoff_ts, function(x) mean(x$LORA))
plot(A2LDM[,"mean_runoff"])


