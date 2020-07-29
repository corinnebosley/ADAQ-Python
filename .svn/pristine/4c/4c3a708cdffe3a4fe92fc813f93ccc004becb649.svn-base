'''
Python code which demonstrates two methods for plotting a vertical
cross section through a plume.

  1. In the first method the line of the cross section is specified
     by listing all the latitudes and longitudes along the
     cross-section
  2. In the second method the line of the cross section is
     specified by providing the start and the end of the cross-
     section as well as the number of sample points along the
     cross-section.

To use this code you only need to copy the method you wish to use
'''

from six.moves.builtins import str
import iris
import iris.plot as iplt
import matplotlib.pyplot as plt
import iris.analysis.trajectory as trajectory

# Directory containing data and name of file to load into IRIS
workdir = '/data/users/apdg/python_sample_data/name/'
filename = 'SO2_grid1_C1_T113_201409061700.txt'

# Name of column to load into IRIS
name = 'Unnamed Field Req 1'

# Create attribute constraint based on name of column
attConstraint = iris.AttributeConstraint(**{'Name':name})

# Load data into a cube
cube = iris.load_cube(workdir+'/'+filename, attConstraint)

#----------------------
# Option 1: Specify the trajectory by listing all the latitudes and
# longitudes along it
sample_points = [('latitude', [52, 54, 56, 58]),
                 ('longitude', [-20, -15, -10, -5])]

# Extract the data at the trajectory points by linearly interpolating
# from the nearest points
interpolated_cube = trajectory.interpolate(cube, sample_points)

# Plot
# Choose the coordinates to be shown on the x and y-axis, these are passed to
# contourf as coords = [x-coord, y-coord]
altitude = interpolated_cube.coord('altitude')
latitude = interpolated_cube.coord('latitude')

cf = iplt.contourf(interpolated_cube, coords=[latitude, altitude])
plt.xlabel(latitude.name().title() + ' (' + str(latitude.units) + ')')
plt.ylabel(altitude.name().title() + ' (' + str(altitude.units) + ')')
plt.title(interpolated_cube.name() + ' (' + str(interpolated_cube.units) + ')')
cb = plt.colorbar(cf, orientation='vertical', shrink=0.7, format='%.1e')

plt.show()


#-------------------------
# Option 2: Specify the trajectory (straight-line) by giving the
# start and end of the line and specifying how many points there
# should be on the line (in this case 30)
waypoints = [{'latitude': 52, 'longitude': -20},
             {'latitude': 58, 'longitude': -5}]
sample_points2 = trajectory.Trajectory(waypoints, sample_count=30)

# This bit is slightly odd - it is necessary to convert the sample points from a
# dictionary to a tuple (found this out here:
# https://ocefpaf.github.io/python4oceanographers/blog/2015/04/06/bathymetry/)
lon = [d['longitude'] for d in sample_points2.sampled_points]
lat = [d['latitude'] for d in sample_points2.sampled_points]

sampled_points = [('longitude', lon),
                  ('latitude', lat)]

# Again extract the data from the cube at the sample points
section = trajectory.interpolate(cube, sampled_points)

# Choose the coordinates to be shown on the x and y-axis, these are
# passed to contourf as coords = [x-coord, y-coord]
altitude = interpolated_cube.coord('altitude')
longitude = interpolated_cube.coord('longitude')

# And plot
cf = iplt.contourf(section, coords=[longitude, altitude])
plt.xlabel(longitude.name().title() + ' (' + str(longitude.units) + ')')
plt.ylabel(altitude.name().title() + ' (' + str(altitude.units) + ')')
plt.title(section.name() + ' (' + str(section.units) + ')')
cb = plt.colorbar(cf, orientation='vertical', shrink=0.7, format='%.1e')
plt.show()
