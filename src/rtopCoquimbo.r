# Top-Kriging interpolation of GloFAS discharge forecast onto a stream network
# Dependencies: sf, ncdf4, rtop

install.packages('rtop')
# 1. Load libraries
library(sf)
library(ncdf4)
library(rtop)

# 2. Read the stream network (riverQ.gpkg)
network <- st_read("riverQ.gpkg")
# Ensure network is in geographic CRS (lon/lat)
network <- st_transform(network, crs = 4326)

# 3. Read GloFAS NetCDF forecast
nc <- nc_open("GloFAS_forecast.nc")
# Inspect variable names
# print(nc)
# Assume variable "discharge" with dimensions [lon, lat, time]
lon <- ncvar_get(nc, "longitude")
lat <- ncvar_get(nc, "latitude")
time <- ncvar_get(nc, "time")  # units e.g. hours since...
# Select a forecast time slice (e.g., first time step)
discharge_slice <- ncvar_get(nc, "discharge", start = c(1,1,1), count = c(-1,-1,1))
nc_close(nc)

# 4. Convert gridded forecast to point observations
# Create a grid of coordinates
grid_pts <- expand.grid(lon = lon, lat = lat)
# Extract values and drop NAs
vals <- as.vector(discharge_slice)
obs_df <- data.frame(grid_pts, discharge = vals)
obs_df <- obs_df[!is.na(obs_df$discharge), ]
# Convert to sf points
oobservations <- st_as_sf(obs_df, coords = c("lon","lat"), crs = 4326)

# 5. Generate prediction locations along network
# Extract nodes (endpoints and junctions)
# Convert lines to nodes
tnodes <- st_union(st_geometry(network)) |> st_node()  # requires lwgeom
# Alternatively sample points every fixed distance
tpred_pts <- st_line_sample(network, density = 1/10000)  # one point per 10 km
predictions <- st_cast(tpred_pts, "POINT")
pred_sf <- st_as_sf(predictions)
colnames(pred_sf)[1] <- "geometry"
st_crs(pred_sf) <- 4326

# 6. Set rtop parameters: use network distance (gDist) and no area support
params <- list(gDist = TRUE, cloud = FALSE)

# 7. Build rtop object
rtop_obj <- createRtopObject(observations = oobservations,
                             predictionLocations = pred_sf,
                             params = params)
# 8. Fit variogram model
top_variogram <- rtopFitVariogram(rtop_obj)

# 9. Perform top-kriging
top_krige <- rtopKrige(rtop_obj)

# 10. Join predictions back to network for visualization
pred_values <- data.frame(pred_sf, estimated = top_krige$pred)
# interpolate back onto lines via nearest point
network_pts <- st_nearest_feature(network, pred_sf)
pred_on_network <- cbind(network, estimated = pred_values$estimated[network_pts])

# 11. Plot results
plot(st_geometry(network), col = "lightblue", lwd = 2)
plot(pred_on_network["estimated"], add = TRUE, lwd = 4)

# 12. Save output
st_write(pred_on_network, "network_discharge_topkrige.gpkg", delete_dsn = TRUE)
