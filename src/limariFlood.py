# import generic packages
import numpy as np
from matplotlib import pyplot as plt

# import geospatial packages
import rasterio
from rasterio.plot import show
from shapely.geometry import LineString
import geopandas as gpd
from shapely.geometry import Point
import time

# import landlab components
from landlab import RasterModelGrid, imshow_grid
from landlab.components.overland_flow import OverlandFlow
from rasterio.features import rasterize

# Open raster image 
rasterObj = rasterio.open('../Rst/DEMLimari200Fill.tif')
limari = gpd.read_file('/Users/carlos/Documents/GitHub/24bp273032/geodata/Limari.gpkg')
mask_polygon = limari.unary_union

show(rasterObj)

elevArray = rasterObj.read(1)

#create grid from raster attributes
nrows = rasterObj.height  # number of raster cells on each side, initially 150
ncols = rasterObj.width
dxy = (rasterObj.transform.a,-rasterObj.transform.e)  # side length of a raster model cell, or resolution [m], initially 50

# define a landlab raster
mg = RasterModelGrid(shape=(nrows, ncols), 
                     xy_spacing=dxy,
                     # Use this line instead:
                     xy_of_lower_left=(rasterObj.bounds[0], rasterObj.bounds[1]))
mg.crs = rasterObj.crs
# Build a mask for active nodes: True if the node falls inside the polygon
transform = rasterObj.transform
mask_shape = (rasterObj.height, rasterObj.width)

# Create a binary mask (1 where inside polygon, 0 where outside)
binary_mask = rasterize(
    [(mask_polygon, 1)],  # List of (geometry, value) pairs
    out_shape=mask_shape,
    transform=transform,
    fill=0,               # Value for areas outside the polygon
    dtype='uint8'
)

# Ravel the mask to match node ordering, flip as needed to match your elevation
binary_mask_raveled = binary_mask[::-1,:].ravel()

# Mark nodes as inactive where mask is 0
inactive_nodes = np.where(binary_mask_raveled == 0)[0]
mg.status_at_node[inactive_nodes] = 4
# show number of rows, cols and resolution
print(nrows, ncols, dxy)

# create a dataset of zero values
zr = mg.add_zeros("topographic__elevation", at="node")

# apply cell elevation to defined arrray
zr += elevArray[::-1,:].ravel()
print(zr)
imshow_grid(mg, "topographic__elevation", shrink=0.5)

#set and apply and initial height
initialHeight = 0.25
depthArray = np.ones(elevArray.shape)*initialHeight
mg.at_node["surface_water__depth"] = depthArray

#define the flood objeds
of = OverlandFlow(mg, steep_slopes=True, 
                  rainfall_intensity=0,
                  mannings_n=0.01,)
                  # Manning's n value)

#list to store times
dtList = []

#Run once and store elapsed time
of.run_one_step()
dtList.append(of.dt)
print(of.dt)

# explore the output data and location

# model outputs
print(of.output_var_names)

# where this nodes are locates
print(of.var_loc("surface_water__depth"))

# show the water depth array
print(mg.at_node["surface_water__depth"])

# plot the resulting water depth for the first run
imshow_grid(mg, "surface_water__depth", shrink=0.5)

total_time = 24*3600  # Total simulation time in seconds
elapsed_time = 0
dtList = []
time_steps = 0

# Define start and end rainfall intensity values:
# Here using 0 to 0.020/3600/24 as before
rain_start = 0
rain_end = 0.020/3600/24

start_time = time.time()
while elapsed_time < total_time:
    # Compute the current rainfall intensity as a linear function of elapsed time.
    # (Alternatively, adjust this calculation to your desired temporal evolution.)
    fraction = elapsed_time / total_time
    current_intensity = rain_start + fraction * (rain_end - rain_start)
    of.rainfall_intensity = current_intensity
    
    # Run the model for one time step
    of.run_one_step()
    
    # Update elapsed time and store the dt
    dtList.append(of.dt)
    elapsed_time += of.dt
    time_steps += 1
    
    print("Time step: %.2f seconds. Elapsed time %.2f seconds" % (of.dt, elapsed_time))

end_time = time.time()
print("Total elapsed time: %.2f seconds" % (end_time - start_time))
    # plot the resulting water depth for the 101th run
# imshow_grid(mg, "surface_water__depth", shrink=0.5,  vmax=10)

# #convert the resulting data to a numpy array
# zArray = mg.at_node["surface_water__depth"].reshape((nrows,ncols))[::-1,:]
# print(zArray)

def export_surface_water_discharge():
    discharge_link = mg.at_link["surface_water__discharge"]
    node_discharge = np.zeros(mg.number_of_nodes, dtype=float)
    node_count = np.zeros(mg.number_of_nodes, dtype=float)

    head = mg.node_at_link_head
    tail = mg.node_at_link_tail

    for i in range(mg.number_of_links):
        node_discharge[ head[i] ] += discharge_link[i]
        node_count[ head[i] ] += 1
        node_discharge[ tail[i] ] += discharge_link[i]
        node_count[ tail[i] ] += 1

    # Average the accumulated discharge values for nodes that are connected to more than one link
    mask = node_count > 0
    node_discharge[mask] /= node_count[mask]

    # --- Reshape the node array to a 2D grid ---
    # For a RasterModelGrid the nodes are arranged in an [n_rows, n_cols] array.
    # When the original elevation was assigned, the raster array was reversed vertically.
    # To match the original projection, flip vertically.
    mapped_discharge = node_discharge.reshape(mg.number_of_node_rows, mg.number_of_node_columns)
    mapped_discharge = np.flipud(mapped_discharge)

    mapped_discharge*=250
    # --- Export to GeoTIFF using the original raster's metadata ---
    import rasterio

    # Use a copy of the original raster's metadata; update the dtype and number of bands
    out_meta = rasterObj.meta.copy()
    out_meta.update({
        "driver": "GTiff",
        "height": mapped_discharge.shape[0],
        "width": mapped_discharge.shape[1],
        "count": 1,
        "dtype": "float32",
        "compress": "lzw"
    })

    # Write the mapped discharge to a GeoTIFF file
    with rasterio.open("discharge.tif", "w", **out_meta) as dst:
        dst.write(mapped_discharge.astype("float32"), 1)

# export_surface_water_discharge()