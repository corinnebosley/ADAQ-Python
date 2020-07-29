'''
Test code to visualise a plume as the sum of its individual
particles. An example is shown in the gallery

This codes is designed to be used on large numbers of particle
trajectories stored by time rather than by particle. An example
is stored in the python sample data directory

.. todo ::
   * Make the tick labels identical on all three plots
   * Add code with works out the ratio of the x and y dimensions
     and scales the figure size accordingly
   * Add code to automatically determine the z-coordinate

Last updated: 2 July 2020 (Susan)

'''

import matplotlib.pyplot as plt
import matplotlib.colors as mplcolors
import iris.plot as iplt
import iris
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import glob

workdir = '/data/users/apdg/python_sample_data/name_trajectory/'
filenames = sorted(glob.glob(workdir + 'Trajectories_C1_T*_*.txt'))

# Plot limits
hlevels = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
hlim = [0, 10]
extent = [13.5, 24.5, 46.3, 53.3]
dateline = 1
clon = 180

# Set up a colormap and normalisation
cmap = 'autumn_r'
norm = mplcolors.BoundaryNorm(hlevels, 256)

for filename in filenames:
    # Load in data
    cube = iris.load_cube(filename, 'U Turb')

    # Convert height coord to km
    cube.coord('height').convert_units('km')

    if dateline == 1:
        # Move longitude points to range 0:360
        lon_points = cube.coord('longitude').points
        lon_points = [lp + 360 if lp < 0 else lp for lp in lon_points]
        cube.coord('longitude').points = lon_points
        clon = 180

    # Set up figure
    fig = plt.figure(figsize=[11, 7])

    # Scatter plot the plume showing altitude by colour
    ax1 = plt.subplot2grid((3, 3), (0, 0), colspan=2, rowspan=2,
                           projection=ccrs.PlateCarree(central_longitude=clon))
    cf = iplt.scatter(cube.coord('longitude'), cube.coord('latitude'),
                      s=20, c=cube.coord('height').points,
                      edgecolor='', cmap=cmap, norm=norm)
    ax1.coastlines('10m')
    ax1.gridlines()
    ax1.set_extent(extent)
    gl = ax1.gridlines(draw_labels=True,
                       linewidth=0.8,
                       alpha=0.9, zorder=9)

    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    # Scatter plot a cross-section through the plume
    ax2 = plt.subplot2grid((3, 3), (2, 0), colspan=2)
    iplt.scatter(cube.coord('longitude'), cube.coord('height'),
                 s=20, c=cube.coord('height').points, edgecolor='',
                 cmap=cmap, norm=norm)
    plt.ylim(hlim)
    plt.xlim(extent[:2])
    plt.colorbar(cf, orientation='horizontal')

    # Scatter plot a cross-section through the plume
    ax3 = plt.subplot2grid((3, 3), (0, 2), rowspan=2)
    iplt.scatter(cube.coord('height'), cube.coord('latitude'),
                 s=20, c=cube.coord('height').points, edgecolor='',
                 cmap=cmap, norm=norm)
    plt.xlim(hlim)
    plt.ylim(extent[2:])
    ax3.yaxis.tick_right()

    # Compute time bounds
    t_coord = cube.coord('time')
    t_stamp1 = t_coord.units.num2date(t_coord.points[0])
    t_fmt = '%H:%M %d/%m/%Y'
    title = 'Valid at {}'.format(t_stamp1.strftime(t_fmt))
    plt.suptitle(title)

    datestring = filename.strip('.txt').split('_')[-1]

    plt.show()

    #plotdir = workdir.replace('Output', 'Plots')
    #pltname = '{}Trajectory_plume_{}.png'.format(plotdir, datestring)
    #print('Plotting {}'.format(pltname))
    #plt.tight_layout(rect=[0.02, 0.03, 1, 0.95])
    #plt.savefig(pltname)
    #plt.close()
