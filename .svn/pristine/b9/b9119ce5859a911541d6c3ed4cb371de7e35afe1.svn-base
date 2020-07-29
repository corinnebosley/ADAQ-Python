# -*- coding: iso-8859-1 -*-
"""
Rings.py
Authors: Frankie Beckett, Susan Leadbetter, Peter Killick
Last Modified: 21 January 2014

.. note::
   Imports distance.py from adaqcode/cube_statistics.py

Program to plot rings showing the maximum distance at which a threshold concentration is predicted
Can plot a number of different rings (e.g. for different fields)

.. todo::
    * Generalise - e.g. different NAME fields, different thresholds, different source locations
    * Add ability to change projection of plot
    * Ability to plot different thresholds for one field

"""

from six.moves.builtins import range
import os.path
import sys
#Find adaq directory -
# should be in directory above the one containing this file
adaq_path = os.path.dirname(os.path.realpath(__file__))+'/../'
sys.path.append(adaq_path)
from adaqcode import cube_statistics

import iris
import numpy as np
import matplotlib.pyplot as plt
from pyproj import Geod
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader

filename='/project/NAME/apnm/NAMEIIIout/PythonHackathon/DataFiles/Fields_grid23_C1_T10_201004161600.txt'

# Define source latitude and longitude
reflong = -19.62
reflat  =  63.63

# Define threshold concentration
threshold = 1E-25

# Load cube
cubes=iris.load(filename)

# Set up plot axes
fig = plt.figure(figsize=[12,4.8])
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([-179, 179, -20, 88])#xmin xmax ymin ymax

# Read in cubes, each represents a different field
for bin1 in cubes:

   # Calculate distance to each point
   bin1_distance = cube_statistics.distance(bin1, reflong, reflat)

   # Extract and flatten (reduce to a single dimension) the data
   cube_data = bin1_distance.data.flat
   cube_distance = bin1_distance.coord('great circle distance').points.flat

   # Remove data values less than the threshold
   dmask = np.where( cube_data >= threshold )
   new_data = cube_data[dmask]
   gcdist = cube_distance[dmask]

   # Set up projection for range ring calculation
   lons = cubes[0].coord('longitude')
   globe = lons.coord_system.as_cartopy_globe()
   #elipse = lons.coord_system.as_cartopy_crs().proj4_params['ellps']
   g = Geod(a=globe.semimajor_axis, b=globe.semiminor_axis)

   # Calculate latitudes and longitudes of range ring
   circ_lon = []
   circ_lat = []
   for azimuth in range(0, 361):
       lon,lat,_ = g.fwd(reflong,reflat,azimuth,max(gcdist),radians=False)
       circ_lon.append(lon)
       circ_lat.append(lat)

   # Plot ring
   plt.plot(circ_lon,circ_lat,label=bin1.attributes['Sources'],
            transform=ccrs.Geodetic())

# Plot country outlines
countries = cfeature.NaturalEarthFeature(
    category='cultural',
    name='admin_0_countries',
    scale='50m',
    facecolor='none')
ax.add_feature(countries, edgecolor='black',zorder=2)

# Set-up the gridlines
gl = ax.gridlines(draw_labels=True,
                  linewidth=0.8,
                  alpha=0.9)
plt.legend(loc=0)

plt.tight_layout(pad=2)

plt.show()





