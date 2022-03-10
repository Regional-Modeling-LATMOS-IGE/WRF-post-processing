# Example plotting script for wrf output data. This script plots 2m temperature.

# Shaddy Ahmed
# 10/03/2022

#Load packages
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from mpl_toolkits.basemap import Basemap
import numpy as np
from wrf import to_np, getvar, smooth2d, get_basemap, latlon_coords, ll_to_xy

# Open the NetCDF file
ncfile = Dataset("/path/to/WRF/output/file/wrfout_d0X_XXXX-XX-XX_XX:XX:XX")

# Get the 2m temperature
T2 = getvar(ncfile, "T2")

# Get the latitude and longitude points
lats, lons = latlon_coords(T2)

# Get the basemap object
bm = get_basemap(T2, resolution='h')

# Create a figure
fig = plt.figure(figsize=(12,9))

# Add geographic outlines
bm.drawcoastlines(linewidth=0.5)
bm.drawstates(linewidth=0.25)
bm.drawcountries(linewidth=0.5)

# Add meridians and parallels
latcorners=[40.,90.]
loncorners=[0,360]
parallels = np.arange(int(np.min(latcorners)),int(np.max(latcorners)),4)
meridians = np.arange(int(np.min(loncorners)),int(np.max(loncorners)),20)

bm.drawmeridians(meridians,color='#6D5F47',labels=[True,True,True,True],linewidth=0.45, fontsize=14)
bm.drawparallels(parallels,color='#6D5F47',labels=[False,False,False,False],linewidth=0.45,fontsize=14)

# Convert the lats and lons to x and y.  Make sure you convert the lats and
# lons to numpy arrays via to_np, or basemap crashes with an undefined
# RuntimeError.
x, y = bm(to_np(lons), to_np(lats))

# Draw the contours and filled contours (x, y, var, number of color levels, color map)
bm.contourf(x, y, to_np(T2), 30, cmap=get_cmap("coolwarm"))

# Add a color bar
plt.colorbar(pad=0.1)

# Add figure title
plt.title("2m Temperature (K)",pad=30,fontsize=14)

# Show plot
fig.tight_layout()
plt.show()

