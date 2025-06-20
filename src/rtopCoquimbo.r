# Top-Kriging interpolation of GloFAS discharge forecast onto a stream network
# Dependencies: sf, ncdf4, rtop

# 1. Load libraries
library(sf)
library(ncdf4)
library(rtop)
library(raster)
library(lwgeom)
library(rtop)
# 2. Read the stream network (riverQ.gpkg)
network <- st_read("../geodata/riverQ.gpkg")
# Ensure network is in geographic CRS (lon/lat)
network_proj <- st_transform(network, 32719)
network_lines <- st_cast(st_geometry(network_proj), "LINESTRING")

# 3. Read GloFAS NetCDF forecast

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
discharge_slice <- projectRaster(
  r,
  crs    = st_crs(network)$proj4string,
  method = "bilinear"
)

# 4. Convert gridded forecast to point observations
# Create a grid of coordinates
grid_pts <- expand.grid(lon = lon, lat = lat)
# Extract values and drop NAs
vals <- as.vector(discharge_slice)
# Suppose you have a mask for valid points
mask <- !is.na(vals)
obs_df <- data.frame(grid_pts[mask, ], discharge = vals[mask])
obs_df <- obs_df[!is.na(obs_df$lon) & !is.na(obs_df$lat), ]
oobservations <- st_as_sf(obs_df, coords = c("lon", "lat"), crs = 4326)
oobservations <- st_transform(oobservations, 32719)
oobservations_sp <- as(oobservations, "Spatial")
oobs_sf <- st_as_sf(oobservations_sp)
oobs_buf_sf <- st_buffer(oobs_sf, dist = 1000)
oobservations_poly_sp <- as(oobs_buf_sf, "Spatial")

# 5. Generate prediction locations along network
# Extract nodes (endpoints and junctions)
# Convert lines to nodes
# Union and extract only linear features
geom_union <- st_union(st_geometry(network))
linear <- geom_union[st_geometry_type(geom_union) %in% c("LINESTRING", "MULTILINESTRING")]

# If linear is empty, try casting:
if (length(linear) == 0) {
  linear <- st_cast(geom_union, "MULTILINESTRING")
}

# Now node
tnodes <- st_node(linear)

# Alternatively sample points every fixed distance
# Cast to LINESTRING specifically
network_linestring <- st_cast(st_geometry(network), "LINESTRING")

# Now sample
# Now convert to sp
# 6. Set rtop parameters: use network distance (gDist) and no area support
params <- list(gDist = TRUE, cloud = FALSE)
# Use centroids of network segments as prediction points
pred_sf <- st_sf(
  id = 1:nrow(network_proj),
  geometry = st_centroid(st_geometry(network_proj))
)
pred_sf_sp <- as(pred_sf, "Spatial")

# 7. Build rtop object
rtop_obj <- createRtopObject(
  observations = oobservations_poly_sp,
  predictionLocations = pred_sf_sp,
  formulaString = "discharge ~ 1",
  variogramModel = list(model = "Exp", psill = 1.0, nugget = 0.1, range = 50000),
  params = list(
    gDist = TRUE,       # Geometric sampling
    debug.level = 1,    # Enable debugging
    nmin = 3,
    nmax = 12,
    maxDist = 15000,
    rresol = 5000       # Reduce resolution if many polygons
  ),
  pdfObs = "norm",
  pdfPred = "norm"
)


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
