# import generic packages
import numpy as np
from matplotlib import pyplot as plt

# import geospatial packages
import rasterio
from rasterio.plot import show
from shapely.geometry import LineString
import geopandas as gpd

# import landlab components
from landlab import RasterModelGrid, imshow_grid
from landlab.components.overland_flow import OverlandFlow
from landlab.components.overland_flow import KinwaveImplicitOverlandFlow

# Open raster image 
rasterObj = rasterio.open('../Rst/DEMLimari500.tif')
show(rasterObj)

elevArray = rasterObj.read(1)

#create grid from raster attributes
nrows = rasterObj.height  # number of raster cells on each side, initially 150
ncols = rasterObj.width
dxy = (rasterObj.transform.a,-rasterObj.transform.e)  # side length of a raster model cell, or resolution [m], initially 50

# define a landlab raster
mg = RasterModelGrid(shape=(nrows, ncols), 
                     xy_spacing=dxy,
                     #xy_of_lower_left=(rasterObj.bounds[0],rasterObj.bounds[1]))
                     xy_of_lower_left=(0,0))

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
                  rainfall_intensity=0)

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

total_time = 3600  # Total simulation time in seconds
elapsed_time = 0
dtList = []
time_steps = 0

# Define start and end rainfall intensity values:
# Here using 0 to 0.020/3600/24 as before
rain_start = 0
rain_end = 0.020/3600/24

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

    # plot the resulting water depth for the 101th run
imshow_grid(mg, "surface_water__depth", shrink=0.5,vmax=5)

#convert the resulting data to a numpy array
zArray = mg.at_node["surface_water__depth"].reshape((nrows,ncols))[::-1,:]
print(zArray)
