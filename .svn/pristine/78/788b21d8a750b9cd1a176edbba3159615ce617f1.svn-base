"""
Code for producing time-series plots.
Contains class for producing the plot, as well as
generic functions and colours variable which
may be useful for other plotting.
"""
from __future__ import division

import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import iris.plot as iplt
import numpy as np
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.feature as cfeature

import line_plot
from field_plot import compute_grid_line_locs
from plotting_functions import insert_logo

if not os.getenv('DISPLAY'):
    #Enable DISPLAY if running under cron
    mpl.use('Agg')

Z_COORD_YLABEL = {'flight_level': 'FL',
                  'altitude': 'm asl',
                  'height': 'm agl',
                  'Z (Pa)': 'Pa'
                 }

def compute_traj_direction(cube_list):
    '''
    Code to compute the direction of the trajectories: forwards or
    backwards for the plot title. If files contain a mixture of forwards
    and backwards trajectories the title will be set to "Trajectories"
    '''
    forward = 0
    backward = 0
    for cube in cube_list:
        time = cube.coord('time')
        timediff = time.points[1] - time.points[0]

        if timediff >= 0:
            forward += 1
        else:
            backward += 1

    if forward > 0 and backward > 0:
        direct_title = 'Trajectories'
    elif forward > 0:
        direct_title = 'Forward Trajectories'
    else:
        direct_title = 'Backward Trajectories'

    return direct_title

def compute_trajectory_extent(cube_list, dateline, meridian,
                              minheight=2.2):
    '''
    Code to compute the extent of a trajectory plot

    :param cube_list: List of iris cubes to calculate extent for.
    :param dateline: Set to 1 if points cross dateline (lon=180deg), otherwise
                     set to 0.
    :param meridian: Set to 1 if points cross meridian (lon=0deg), otherwise
                     set to 0.
    :param minheight: Minimum height of the rectangle used in the plot.
    :returns: extent: [xmin, xmax, ymin, ymax]

    '''

    # First set min and max using first cube
    xmin = min(cube_list[0].coord('longitude').points)
    xmax = max(cube_list[0].coord('longitude').points)
    ymin = min(cube_list[0].coord('latitude').points)
    ymax = max(cube_list[0].coord('latitude').points)

    for cube in cube_list:
        lon_points = cube.coord('longitude').points
        lat_points = cube.coord('latitude').points

        xmin = min(xmin, min(lon_points))
        xmax = max(xmax, max(lon_points))
        ymin = min(ymin, min(lat_points))
        ymax = max(ymax, max(lat_points))

    dx = np.max(np.diff(lon_points))
    x_range = xmax - xmin
    y_range = ymax - ymin
    xmid = xmin + x_range/2
    ymid = ymin + y_range/2

    xfactor = 1.6

    xy_range = max([25.0*dx, x_range/xfactor, y_range])*1.1
    xy_range = max([xy_range, minheight])

    ymin = min([ymin, ymid-xy_range*0.5])
    ymax = max([ymax, ymid+xy_range*0.5])
    xmin = min([xmin, xmid-xy_range*xfactor*0.5])
    xmax = max([xmax, xmid+xy_range*xfactor*0.5])

    # shift whole area to be within global limits

    if dateline == 1 and meridian == 0:
        lonmax = 360
        lonmin = 0
    else:
        lonmax = 180
        lonmin = -180

    if xmax > lonmax:
        xshift = xmax - lonmax
        xmax = xmax - xshift
        xmin = xmin - xshift

    if xmin < lonmin:
        xshift = xmin - lonmin
        xmax = xmax - xshift
        xmin = xmin - xshift

    xmax = min([xmax, lonmax])
    xmin = max([xmin, lonmin])

    if ymax > 90:
        yshift = ymax - 90.0
        ymax = ymax - yshift
        ymin = ymin - yshift

    if ymin < -90.0:
        yshift = ymin + 90
        ymax = ymax - yshift
        ymin = ymin - yshift

    ymax = min([ymax, 90.0])
    ymin = max([ymin, -90.0])

    extent = [xmin, xmax, ymin, ymax]
    return extent

#--- Class ---
class TrajectoryPlot(line_plot.LinePlot):
    """
    Class for plotting trajectories, a subclass of LinePlot.

    """

    def __init__(self):
        """
        Initiates class as a subset of line_plot.LinePlot
        """

        line_plot.LinePlot.__init__(self)
        self.extent = []
        self.clon = 0
        self.mapping = 'coastlines' #Mapping to add
        self.mobrand = False
        self.release_info = None
        self.rsmc = False
        self.annote = False

    def plot(self):
        """
        Produce trajectory plot.

        Returns fig object for further plotting if needed.
        """

        if not self.lines:
            raise ValueError("TrajectoryPlot: no lines have been added")

        if self.fig is None:
            if self.rsmc:
                self.fig = plt.figure(figsize=[12, 6])
                ax = plt.subplot2grid(
                    (3, 3), (0, 0), rowspan=2, colspan=2,
                    projection=ccrs.PlateCarree(central_longitude=self.clon))
            elif self.annote:
                self.fig = plt.figure(figsize=[12, 6])
                ax = plt.subplot2grid(
                    (3, 3), (0, 0), rowspan=2, colspan=2,
                    projection=ccrs.PlateCarree(central_longitude=self.clon))
            else:
                self.fig = plt.figure(figsize=[7, 9])
                ax = plt.subplot2grid(
                    (3, 1), (0, 0), rowspan=2,
                    projection=ccrs.PlateCarree(central_longitude=self.clon))
        ax = plt.gca()

        for line in self.lines:

            add_settings = {}
            if 'add_settings' in line:
                add_settings = line['add_settings']

            style = {'label':line['label'], 'color':line['colour'],
                     'linestyle':line['linestyle'],
                     'linewidth':line['linewidth'], 'marker':line['marker']}
            style2 = style.copy()
            style2.update(add_settings)

            iplt.plot(line['x'], line['y'], **style2)

            # Add a black square at the trajectory start point
            iplt.scatter(line['x'][0], line['y'][0], color='k', marker='s')

        #Add title and axis labels
        if self.title is not None:
            ax.set_title(self.title)

        # Set the extent
        # Bug in Cartopy Dec'17 - Global extent will not be plotted with
        # extent[0] = 0, extent[1] = 360 So the longitudinal extents are
        # deliberately taken in by 0.1
        if abs(self.extent[1] - self.extent[0]) > 330:
            self.extent[0] = -179.9
            self.extent[1] = 179.9
        ax.set_extent(self.extent, crs=ccrs.PlateCarree())

        # Determine extent of plotting region and use this to
        # select an appropriate mapping zoom
        if abs(self.extent[1] - self.extent[0]) < 15.0:
            res = '10m'
        elif abs(self.extent[1] - self.extent[0]) < 50.0:
            res = '50m'
        else:
            res = '110m'

        if self.mapping == 'countries' or self.mapping == 'states':
            countries = cfeature.NaturalEarthFeature(
                category='cultural',
                name='admin_0_countries_lakes',
                scale=res,
                facecolor='none')

        if self.mapping == 'states':
            states = cfeature.NaturalEarthFeature(
                category='cultural',
                name='admin_1_states_provinces_shp',
                scale=res,
                facecolor='none')

        if self.mapping == 'coastlines':

            ax.coastlines(res)


        elif self.mapping == 'countries':

            ax.coastlines(res, zorder=3)
            ax.add_feature(countries, edgecolor='gray',
                           zorder=2, linewidth=0.5)

        elif self.mapping == 'states':

            ax.coastlines(res, zorder=3)
            ax.add_feature(countries, edgecolor='gray',
                           zorder=2, linewidth=1)
            ax.add_feature(states, edgecolor='lightgray',
                           zorder=2, linewidth=0.5)

        elif self.mapping == 'wms':
            # NOTE WMS mapping does not appear to work for extents
            # greater than 130 degrees in either direction for a
            # typical 6x6 sized map. For smaller maps, the useable
            # WMS extents are smaller.
            # It should also be noted that if the WMS map
            # crosses 180E/W, if the northern or southern
            # edge of the map is on the equator, this will
            # result in the size of page and the placing of
            # the map on the page being altered.
            num_layers = np.linspace(0, 40, 41)
            layers = ['{:.0f}'.format(x) for x in num_layers]
            ax.add_wms(wms='http://exxdmmsprd01:6080/arcgis/services/DMMS/' +
                       'Global_NE_HC_Hybrid_Greyscale/MapServer/WMSServer',
                       layers=layers)

        #Add gridlines
        if self.gridlines:
            try:
                if self.extent[0] < 180 and self.extent[1] > 180:
                    xlocs, xlocs_extend = compute_grid_line_locs(self.extent)
                    ax.gridlines(xlocs=xlocs_extend)
                    gl = ax.gridlines(draw_labels=True,
                                      xlocs=xlocs,
                                      linewidth=0.001)
                else:
                    gl = ax.gridlines(draw_labels=True,
                                      linewidth=0.8,
                                      alpha=0.9, zorder=9)

                gl.xlabels_top = False
                gl.ylabels_right = False
                gl.xformatter = LONGITUDE_FORMATTER
                gl.yformatter = LATITUDE_FORMATTER
            except:
                gl = ax.gridlines()

        # Add information about the release location and time if provided
        if self.release_info is not None:
            release_text = 'Release location: {}, {},\n'.format(
                self.release_info[0], self.release_info[1])
            release_text += 'Release time: {}'.format(
                self.release_info[2])
            ax.annotate(release_text,
                        xy=(0.5, 0.34),
                        xycoords=('axes fraction', 'figure fraction'),
                        xytext=(0, 10),
                        textcoords='offset points',
                        size=12, ha='center', va='bottom')

        # Apply branding
        if self.mobrand:
            insert_logo()

        return self.fig

if __name__ == '__main__':

    import doctest
    doctest.testmod()
