"""
field_plot.py

Contains code required for plotting layers, adding maps, titles and
colorbars
"""
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import os
import unicodedata

from six.moves.builtins import str
from six.moves.builtins import object

import iris.plot as iplt
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.io.img_tiles import GoogleTiles

# adaqcode
from plotting_functions import insert_logo, resize_fig
from name_functions import height_text_fmt


def compute_grid_line_locs(extent, lat_or_lon='lon', maxbins=8):
    """
    Function to compute the locations of gridlines for plots which
    extend over the dateline. Requires the extent as an input
    and lat_or_lon = 'lat' for latitude lines

    >>> extent = [133, 203, 20, 50]
    >>> xlocs, xlocs_extend = compute_grid_line_locs(extent)
    >>> print(xlocs)
    [140.0, 150.0, 160.0, 170.0, 180.0, -170.0, -160.0]
    >>> ylocs, ylocs_extend = compute_grid_line_locs(extent, lat_or_lon='lat')
    >>> print(ylocs)
    [24.0, 30.0, 36.0, 42.0, 48.0]
    """

    steps = [0.0001, 0.0002, 0.0005,
             0.001, 0.002, 0.005,
             0.01, 0.02, 0.05,
             0.1, 0.2, 0.5,
             1, 2, 3, 6, 10, 15, 30, 45, 60]

    if lat_or_lon == 'lon':
        range_min = extent[0]
        range_max = extent[1]
    if lat_or_lon == 'lat':
        range_min = extent[2]
        range_max = extent[3]
    range1 = range_max - range_min

    nstep = 0
    nbins = range1 / steps[nstep]

    while nbins > maxbins:
        nstep = nstep + 1
        nbins = range1 / steps[nstep]

    step = steps[nstep]
    gg1 = np.ceil(range_min / float(step))
    locs = np.arange(gg1 * step, range_max, step)
    locs_extend = np.arange(gg1 * step - step, range_max + step, step)

    # Ensure all gridline locations lie in range -180:180 (lon) or -90:90 (lat)
    if lat_or_lon == 'lon':
        locs = [x if x <= 180 else x - 360 for x in locs]
        locs = [x + 360 if x < -180 else x for x in locs]
    if lat_or_lon == 'lat':
        locs = locs.clip(-90, 90)
        locs = [y for y in locs]

    return locs, locs_extend


class GreyscaleOSMMapTiles(GoogleTiles):
    """
    Class for using greyscale Ordnance Survey web mapping service.

    """

    def _image_url(self, tile):
        x, y, z = tile
        url = 'http://a.tiles.wmflabs.org/bw-mapnik/{}/{}/{}.png'
        return url.format(z, x, y)


# -----------------------------------------------------------
# Field Plot Object
# -----------------------------------------------------------


class FieldPlot:
    """
    Class for plotting field layer objects

    Read in data

    >>> import name_data
    >>> import config
    >>> from adaqcode import inifile
    >>> ini_dict = inifile.get_inidict(
    ... defaultfilename='adaqcode/inifile.ini') # doctest: +ELLIPSIS
    Reading inifile ...adaqcode/inifile.ini
    >>> sample_data_path = config.SAMPLE_DATADIR+'name/'
    >>> name = name_data.NAMEData()
    >>> name.readdata(sample_data_path + 'Fields_grid3_201511291800.txt',
    ... field_attributes = {'Quantity': 'Total deposition'})
    [<iris 'Cube' of CAESIUM-137_TOTAL_DEPOSITION / (Bq/m2) \
(latitude: 205; longitude: 245)>]

    Set up field_layer class and add cube to it

    >>> from field_layer import FieldLayer
    >>> fl = FieldLayer (name.gridded_cube_list[0])

    Set the layer style (note that pcolormesh is the python equivalent
    of 'pixel' in IDL)

    >>> fl.set_layerstyle(plottype='pcolormesh', colorscale='log',
    ... mask=True, autozoom=True)

    Set up plot and add layer

    >>> fp = FieldPlot(ini_dict)
    >>> fp.add_layer(fl)
    >>> fp2 = FieldPlot(ini_dict)
    >>> fp2.add_layer(fl)
    >>> fp3 = FieldPlot(ini_dict)
    >>> fp3.add_layer(fl)

    Setup some background mapping

    >>> fp.setup_mapping(extent=[-12, 4.0, 48., 60.], mapping='coastlines',
    ... gridlines=True)
    >>> fp2.setup_mapping(extent=[-4, 2, 51., 55.], mapping='wms',
    ... gridlines=True)
    >>> fp3.setup_mapping(extent=[-3., 2., 51., 55.], mapping='grey_os',
    ... gridlines=True)

    Plot

    >>> fp.plot()  # doctest: +ELLIPSIS
    <...>
    >>> fp2.plot()  # doctest: +ELLIPSIS
    <...>
    >>> fp3.plot()  # doctest: +ELLIPSIS
    <...>

    Save to file

    >>> import config
    >>> fp.save_fig(plotdir=config.CODE_DIR + "/adaqdocs/figures")
    ... # doctest: +ELLIPSIS
    Saved figure  .../adaqdocs/figures/Fieldplot_NAME_CAESIUM-137_\
TOTAL_DEPOSITION_Boundarylayer_201511291800.png

    .. image:: ../adaqdocs/figures/Fieldplot_NAME_CAESIUM-137_\
TOTAL_DEPOSITION_Boundarylayer_201511291800.png

    >>> fp2.save_fig(plotdir=config.CODE_DIR + "/adaqdocs/figures/wms")
    ... # doctest: +ELLIPSIS
    Saved figure  .../adaqdocs/figures/wms/Fieldplot_NAME_CAESIUM-137_\
TOTAL_DEPOSITION_Boundarylayer_201511291800.png

    .. image:: ../adaqdocs/figures/wms/Fieldplot_NAME_CAESIUM-137_\
TOTAL_DEPOSITION_Boundarylayer_201511291800.png

    >>> fp3.save_fig(plotdir=config.CODE_DIR + "/adaqdocs/figures/grey_os")
    ... # doctest: +ELLIPSIS
    Saved figure  .../adaqdocs/figures/grey_os/Fieldplot_NAME_CAESIUM-137_\
TOTAL_DEPOSITION_Boundarylayer_201511291800.png

    .. image:: ../adaqdocs/figures/grey_os/Fieldplot_NAME_CAESIUM-137_\
TOTAL_DEPOSITION_Boundarylayer_201511291800.png

    """

    def __init__(self, ini_dict=None):
        """
        Initialises class
        """
        # Figure
        self.fig = None

        # Layers
        self.layer_list = []

        # x-coord and y-coord will only be determined when the
        # field layer is plotted.
        self.x_coord = None
        self.y_coord = None

        # Mapping and extent
        self.projection = None
        self.extent = None
        self.central_longitude = 0
        self.mapping = None
        self.gridlines = False
        self.gridlines_xmaxbins = None
        self.gridlines_ymaxbins = None
        self.gridlines_xlabel_size = None
        self.gridlines_ylabel_size = None

        # Title
        self.title = None
        self.titleprefix = None
        self.suptitle = None

        # Pull some values from ini_dict, otherwise set default values:
        if ini_dict is not None:
            self.title = ini_dict.get('title', None)
            self.longtitle = True if self.title == 'name_verbose' else False
            self.mobrand = ini_dict.get('mobrand', False)
            self.rsmc = ini_dict.get('rsmc', False)
            self.default_figsize = ini_dict.get('figsize', 'Auto')
            self.annote_location = ini_dict.get('annote_location', None)
            self.cbar = ini_dict.get('cbar', None)
            self.cbar_orientation = ini_dict.get('cbar_orientation', None)
        else:
            self.longtitle = False
            self.mobrand = False
            self.rsmc = False
            self.default_figsize = 'Auto'
            self.annote_location = None
            self.cbar = None
            self.cbar_orientation = None

        self.mappable = None

        # Output
        self.plotdir = None
        self.figformat = "png"

    def __str__(self):
        """
        Output from 'print'
        """
        string = str(type(self)) + '\n'
        keys = sorted(self.__dict__.keys())
        for key in keys:
            string += key + ': ' + str(self.__dict__[key]) + '\n'

        return string

    def add_layer(self, layer):
        """
        Add in a new layer

        .. note::
           May also add here the option to add the importance of the layer
        """
        self.layer_list.append(layer)

    def auto_figsize(self, height=None, print_size=False):
        """
        Automatically resize figure based on ratio of x-y coordinates and
        extra components to be added to the plot space.  All components are
        extracted from the FieldPlot object except for height.

        Sets up figure size such that the xsize is scaled according to required
        size of y, based on the number of points in x and y coordinates,
        but allowing some space for colorbar (based on colorbar fraction
        default = 0.4 of axis), and the location of annotation if requested.

        If a colorbar is to be constructed but its orientation is not
        specified, then matplotlib will choose where it goes, so the colorbar
        fraction is increased to 0.6 to make sure that elements of the plot
        are not pushed off the figure.

        :param height: required height of figure.

        :param print_size: boolean value representing whether to print out the
                           new figure size as standard output.  Defaults to
                           False.


        Firstly set up an arbitrary field plot object with a figure and check
        what the default figure size would be:

        >>> import iris
        >>> import matplotlib.pyplot as plt
        >>> from adaqcode import inifile
        >>> ini_dict = inifile.get_inidict(defaultfilename=
        ... 'adaqcode/inifile.ini')  # doctest: +ELLIPSIS
        Reading inifile .../adaqcode/inifile.ini
        >>> x_coord = iris.coords.AuxCoord(np.arange(85),long_name='x_coord')
        >>> y_coord = iris.coords.AuxCoord(np.arange(100),long_name='y_coord')

        Note less points in xcoord (85) than y_coord (100).

        Now initialise a FieldPlot object and some of its properties:

        >>> fplt = FieldPlot(ini_dict)
        >>> fplt.x_coord = x_coord
        >>> fplt.y_coord = y_coord
        >>> fplt.fig = plt.Figure()

        Set the FieldPlot property 'default_figsize' to 'True' to keep
        original (default) sizing:

        >>> fplt.default_figsize = True
        >>> fplt.auto_figsize(print_size=True)
        Resizing figure to  [8.0, 6.0]

        # Set FieldPlot property 'default_figsize' to 'False' to resize
        automatically:

        >>> fplt.default_figsize = False
        >>> fplt.auto_figsize(print_size=True)
        Resizing figure to  [9.18, 10.8]

        Now resize on the basis that a horizonal colour bar will be included.
        This will mean reducing the x coordinate relative to the y coordinate
        in the FieldPlot object:

        >>> fplt.cbar_orientation = 'horizontal'
        >>> fplt.auto_figsize(print_size=True)
        Resizing figure to  [5.1, 8.399999999999999]

        Resize if a vertical colour bar is included:

        >>> fplt.cbar_orientation = 'vertical'
        >>> fplt.auto_figsize(print_size=True)
        Resizing figure to  [7.139999999999999, 6.0]

        Resize if a colorbar is to be plotted but without an orientation:

        >>> fplt.cbar_orientation = None
        >>> fplt.auto_figsize(print_size=True)
        Resizing figure to  [9.18, 10.8]

        Resize if horizontal annotation is included:

        >>> fplt.annote_location = 'below'
        >>> fplt.auto_figsize(print_size=True)
        Resizing figure to  [9.18, 14.96]

        Also try setting height of figure and then resizing (remember to set
        annotation and colorbars back to None to prevent additional figure
        stretching):

        >>> fplt.cbar_orientation = None
        >>> fplt.annote_location = None
        >>> fplt.auto_figsize(height=10.0, print_size=True)
        Resizing figure to  [8.499999999999998, 10.0]

        """
        # Set some default or assumed values to work with.
        figsize = [8.0, 6.0]

        if self.default_figsize in [None, 'Auto', False]:
            cbar_fraction = 0.4

            # y-axis size is set as a standard 6.0 or specified height.
            # x-axis size will be the proportion of x points to y points
            # multiplied by the standard size of the y-axis:
            figsize[0] = len(self.x_coord.points) / len(self.y_coord.points) \
                * figsize[1]
            if self.cbar is False:
                # If there isn't a colorbar, no need to change the figsize.
                figsize = [figsize[0], figsize[1]]
            else:
                if self.cbar_orientation == 'vertical':
                    # If colorbar is next to plot, make figure wider.
                    figsize = [figsize[0] * (1 + cbar_fraction), figsize[1]]
                elif self.cbar_orientation == 'horizontal':
                    # If colorbar is underneath plot, make figure taller.
                    figsize = [figsize[0], figsize[1] * (1 + cbar_fraction)]
                elif self.cbar_orientation is None:
                    # In this case, a colorbar will be plotted anyway, but
                    # matplotlib will select the position and orientation, so
                    # just make the figure bigger in both planes and by a
                    # slightly larger amount to cover matplotlib's decisions.
                    cbar_fraction = 0.8
                    figsize = [figsize[0] * (1 + cbar_fraction),
                               figsize[1] * (1 + cbar_fraction)]

            # Space for branding (logo)
            if self.mobrand is True:
                figsize = [figsize[0] * (6. / 5.), figsize[1] * (3. / 2.)]

            # Space for annotation
            if self.annote_location == 'right':
                figsize = [figsize[0] * (6. / 5.) + 2., figsize[1]]
            elif self.annote_location == 'below':
                figsize = [figsize[0], figsize[1] * (6. / 5.) + 2.]

            # Space for a long title
            if self.longtitle:
                figsize = [figsize[0], figsize[1] + 2.]

        # Ignore all of the above for figsize[1] if height is set manually.
        if height is not None:
            figsize[0] *= (height / figsize[1])
            figsize[1] = height

        if print_size is True:
            print('Resizing figure to ', figsize)

        resize_fig(self.fig, figsize)

    def setup_mapping(self, projection=None, extent=None,
                      mapping=None, gridlines=False,
                      gridlines_xmaxbins=None, gridlines_ymaxbins=None,
                      gridlines_xlabel_size=None, gridlines_ylabel_size=None,
                      central_longitude=0):
        """
        Method to setup various mapping parameters

        :param projection: projection for the plot
        :param extent: plotting extent in WGS84 reference frame in format
                       [xmin,xmax,ymin,ymax]
        :param mapping: choice of background map
                        (currently 'None', 'coastlines', 'countries',
                        or 'states')
        :param gridlines: add gridlines (True or False)
        :param gridlines_xmaxbins: max number of bins for the x gridlines
        :param gridlines_ymaxbins: max number of bins for the y gridlines
        :param gridlines_xlabel_size: text size of labels for x gridlines
        :param gridlines_ylabel_size: text size of labels for y gridlines
        :param central_longitude: centre longitude of plot

        Load an ini dict and use it to set up a plot:

        >>> from adaqcode import inifile
        >>> ini_dict = inifile.get_inidict(
        ... defaultfilename='adaqcode/inifile.ini') # doctest: +ELLIPSIS
        Reading inifile .../adaqcode/inifile.ini
        >>> fplt = FieldPlot(ini_dict)

        Add some mapping

        >>> fplt.setup_mapping(extent=[-12, 8.8, 45, 62.8],
        ... mapping='coastlines')

        >>> for key in sorted(fplt.__dict__.keys()):
        ...    print(key, ':', fplt.__dict__[key])
        annote_location : None
        cbar : None
        cbar_orientation : None
        central_longitude : 0
        default_figsize : Auto
        extent : [-12, 8.8, 45, 62.8]
        fig : None
        figformat : png
        gridlines : False
        gridlines_xlabel_size : None
        gridlines_xmaxbins : None
        gridlines_ylabel_size : None
        gridlines_ymaxbins : None
        layer_list : []
        longtitle : False
        mappable : None
        mapping : coastlines
        mobrand : False
        plotdir : None
        projection : None
        rsmc : False
        suptitle : None
        title : None
        titleprefix : None
        x_coord : None
        y_coord : None

        .. note::
           The projection parameter can be one of 'Mercator', 'PlateCarree',
           'NorthPolarStereo' or 'SouthPolarStereo' or a cartopy coordinate
           reference system
           (e.g. ccrs.LambertConformal(central_longitude=-20.0))

        .. note::
           The central longitude is particularly important for layers which
           extend across the dateline

        """

        if extent is not None:
            self.extent = extent
        if mapping is not None:
            self.mapping = mapping
        if gridlines is not None:
            self.gridlines = gridlines
            if gridlines:
                if gridlines_xmaxbins is not None:
                    self.gridlines_xmaxbins = gridlines_xmaxbins
                if gridlines_ymaxbins is not None:
                    self.gridlines_ymaxbins = gridlines_ymaxbins
                if gridlines_xlabel_size is not None:
                    self.gridlines_xlabel_size = gridlines_xlabel_size
                if gridlines_ylabel_size is not None:
                    self.gridlines_ylabel_size = gridlines_ylabel_size
        if central_longitude is not None:
            self.central_longitude = central_longitude
        if projection is not None:
            projection = projection.lower()

        if projection == 'mercator':
            self.projection = ccrs.Mercator(
                central_longitude=self.central_longitude)
        elif projection == 'osgb':
            self.projection = ccrs.OSGB()
        elif projection == 'platecarree':
            self.projection = ccrs.PlateCarree(
                central_longitude=self.central_longitude)
        elif projection == 'southpolarstereo':
            self.projection = ccrs.SouthPolarStereo(
                central_longitude=self.central_longitude)
        elif projection == 'northpolarstereo':
            self.projection = ccrs.NorthPolarStereo(
                central_longitude=self.central_longitude)
        elif projection == 'lamberticeland':
            self.projection = ccrs.LambertConformal(central_longitude=-18.0)
        elif projection is not None:
            self.projection = projection

    def plot(self, fig=None, figsize=None, field_layers=None):
        """
        Contains all the main plotting commands or calls to them

        Returns the figure object so the user can add to the plot
        if required

        :param fig: matplotlib figure object on which to create the plot.

        :param figsize: Selected size of figure to create if fig is None.

        :param field_layers: field layer with which to create plot.

        """
        if fig is None:
            self.fig = plt.figure()
        else:
            self.fig = fig

        # Examine first layer to extract some properties of the plot:
        if field_layers is not None:
            layer_list = [layer for layer in field_layers]
        else:
            layer_list = [layer for layer in self.layer_list]

        layer = layer_list[0]
        if layer.plottype in ['pcolormesh', 'contour', 'contourf',
                              'contourf-edge']:
            self.x_coord = layer.cube.dim_coords[-2]  # lat-type coord
            self.y_coord = layer.cube.dim_coords[-1]  # lon-type coord
        elif layer.plottype == 'scatter_latlon':
            self.x_coord = layer.cube.coord('longitude')
            self.y_coord = layer.cube.coord('latitude')

        # Amend figsize to fit completed plot unless figsize is set manually:
        if figsize in [None, 'Auto', 'Default']:
            # Automatically determine figure size if it isn't set by user.
            self.auto_figsize()
        else:
            # If figure size is set by user, just make the figure that size.
            resize_fig(self.fig, figsize)

        # Make sure to specify fig to draw on:
        fig = self.fig

        # Set the projection if requested
        if self.projection is not None:
            if self.annote_location == 'right':
                plt.subplot2grid((1, 3), (0, 0), colspan=2,
                                 projection=self.projection)
            elif self.annote_location == 'below':
                plt.subplot2grid((3, 1), (0, 0), rowspan=2,
                                 projection=self.projection)
            else:
                plt.axes(projection=self.projection)
        elif self.annote_location == 'right':
            plt.subplot2grid((1, 3), (0, 0), colspan=2)
        elif self.annote_location == 'below':
            plt.subplot2grid((3, 1), (0, 0), rowspan=2)

        # Find field layer(s) and plot it/them.

        for layer in layer_list:
            self.plot_layer(layer, self.mapping)

        if self.mapping is not None:
            self.plot_map()

        # Add source location if it exists
        coord_names = [coord.name()
                       for coord in self.layer_list[0].cube.coords()]
        if 'source_longitude' in coord_names:
            lyr_cube = self.layer_list[0].cube
            sourcelon = lyr_cube.coord('source_longitude').points[0]
            sourcelat = lyr_cube.coord('source_latitude').points[0]
            if self.rsmc:
                rsmc_cube = self.layer_list[0].cube
                plt.plot(sourcelon, sourcelat, 'ks',
                         markersize=6, label='Source Location',
                         transform=ccrs.PlateCarree())
                maxval = np.argmax(rsmc_cube.data)
                xyvals = np.unravel_index(maxval, rsmc_cube.data.shape)
                lon_max = rsmc_cube.coord('longitude').points[xyvals[0]]
                lat_max = rsmc_cube.coord('latitude').points[xyvals[1]]
                if rsmc_cube.attributes['short_name'] != 'TimeOfArrival':
                    plt.plot(lon_max, lat_max, 'k*',
                             transform=ccrs.PlateCarree())
            else:
                plt.plot(sourcelon, sourcelat, 'r*',
                         markersize=10, label='Source Location',
                         transform=ccrs.PlateCarree())

        # Add a title
        if self.title is None:
            self.title = self.gen_title()
        elif self.title == 'name_verbose':
            self.title = self.gen_name_verbose_title()

        # Add prefix to title if required
        if self.titleprefix is not None:
            self.title = self.titleprefix + self.title
        plt.title(self.title, linespacing=1.5)

        # Add a supertitle
        if self.suptitle is not None:
            plt.suptitle(self.suptitle, fontsize=16)

        # Add branding
        if self.mobrand:
            insert_logo()
            copyright_sign = unicodedata.lookup('COPYRIGHT SIGN')
            copyright_txt = copyright_sign + " Met Office Crown Copyright"
            plt.annotate(copyright_txt,
                         (0, 0), (-50, -50),
                         xycoords='axes fraction',
                         textcoords='offset points',
                         va='top')
        # Ensure that the field_plot object has all the changes made here
        # included in self.fig:
        self.fig = fig
        return self.fig

    def multiplot(self, fig=None):
        """
        Generates a single postage stamp plot over a set of ensemble members.

        Contains all the main plotting commands or calls to them

        Returns the figure object so the user can add to the plot
        if required
        """

        if fig is None:
            self.fig = plt.figure()
        else:
            self.fig = fig

        fig.subplots_adjust(left=0.075, right=0.95, bottom=0.05, top=0.88,
                            hspace=0.3, wspace=0.3)

        # Identify number of ensemble members defined within the cube
        nmembers = len(self.layer_list[0].cube.coord('realization').points)

        # Save state of cbar variable for later use
        plot_colourbar = self.layer_list[0].cbar

        # Add the individual postage stamps as subplots
        for ens_layer in self.layer_list:
            for mem, layer in enumerate(
                    ens_layer.layer_slice_over('realization')):

                # Do not plot colour bar with each postage stamp
                # A single colour bar is added to the plot later on
                layer.cbar = False

                # Specify position of each postage stamp
                # and set the projection if requested
                if self.projection is not None:
                    plt.subplot(4, 5, mem + 3, projection=self.projection)
                else:
                    plt.subplot(4, 5, mem + 3)

                self.plot_layer(layer, self.mapping)

                # Set an axes extent (assuming extent is given in WGS84)
                # Needs to follow plot in order that axes have become geoaxes
                if self.extent is not None:
                    plt.gca().set_extent(self.extent, crs=ccrs.PlateCarree())

                if self.mapping is not None:
                    self.plot_map()

                # Add subtitle for each postage stamp
                mem_title = layer.cube.coord('ensemble_member').points[0]
                plt.title(mem_title, fontsize=10)

                # Add source location if it exists
                coord_names = [coord.name() for coord in layer.cube.coords()]
                if 'source_longitude' in coord_names:
                    cube = layer.cube
                    sourcelon = cube.coord('source_longitude').points[0]
                    sourcelat = cube.coord('source_latitude').points[0]
                    if self.rsmc:
                        plt.plot(sourcelon, sourcelat, 'ks',
                                 markersize=6, label='Source Location',
                                 transform=ccrs.PlateCarree())
                        maxval = np.argmax(cube.data)
                        xyvals = np.unravel_index(maxval, cube.data.shape)
                        lon_max = cube.coord('longitude').points[xyvals[0]]
                        lat_max = cube.coord('latitude').points[xyvals[1]]
                        if cube.attributes['short_name'] != 'TimeOfArrival':
                            plt.plot(lon_max, lat_max, 'k*',
                                     transform=ccrs.PlateCarree())
                    else:
                        plt.plot(sourcelon, sourcelat, 'r*',
                                 markersize=10, label='Source Location',
                                 transform=ccrs.PlateCarree())

        ens_layer = self.layer_list[0]

        # Add colour bar, titles, etc.
        # Set title and suptitle text
        if self.title is not None:
            title = self.title
        else:
            title = ens_layer.cube.attributes['Species'] + ' '
            title += ens_layer.cube.attributes['Quantity']
            cbar_label = ens_layer.cube.attributes['Quantity']

        # Add prefix to title if required
        if self.titleprefix is not None:
            title = self.titleprefix + title

        if self.suptitle is not None:
            suptitle = self.suptitle
        else:
            shortname_text = ens_layer.cube.attributes['short_name']
            if 'Time Av or Int' in ens_layer.cube.attributes:
                tavint_text = ens_layer.cube.attributes['Time Av or Int']
            else:
                tavint_text = ens_layer.cube.attributes['Time av/int info']
            height_text = height_text_fmt(ens_layer.cube)
            t_coord = ens_layer.cube.coord('time')
            t1 = t_coord.units.num2date(t_coord.bounds[0][0])
            t2 = t_coord.units.num2date(t_coord.bounds[0][1])
            t_fmt = '%H:%M %d/%m/%Y'
            suptitle = 'NAME ensemble with {} members\n'.format(nmembers)
            suptitle += '{}, {}, {}\n'.format(shortname_text, height_text,
                                              tavint_text)
            suptitle += 'Valid from {} to {}'.format(t1.strftime(t_fmt),
                                                     t2.strftime(t_fmt))

        # Make a colourbar if required
        if plot_colourbar:
            ens_layer.cbar = plot_colourbar

            cbar_axes = plt.gcf().add_axes([0.14, 0.80, 0.24, 0.02])

            ens_layer.construct_cbar(self.mappable,
                                     position=cbar_axes,
                                     orientation='horizontal',
                                     title=title)

        plt.suptitle(suptitle, fontsize=14, linespacing=1.2)

        # Add branding
        if self.mobrand:
            insert_logo()
            copyright_sign = unicodedata.lookup('COPYRIGHT SIGN')
            copyright_txt = copyright_sign + " Met Office Crown Copyright"
            plt.annotate(copyright_txt,
                         (0, 0), (0, -40),
                         xycoords='axes fraction',
                         textcoords='offset points',
                         va='top')

        return self.fig

    def plot_layer(self, layer, mapping):
        """
        Method for plotting each individual layer using the information
        provided in the layer class.


        :param layer: field layer to plot.

        :param mapping: map type to use for plot background.

        """
        # First make sure matplotlib knows which figure to draw on:
        fig = self.fig
        alpha = 1.0
        if mapping in ['grey_os', 'wms', 'jam']:
            alpha = 0.5

        # Quick plot on a linear scale if only a colormap is given
        if layer.levels is None and layer.colorscale == 'linear':
            if layer.plottype == 'pcolormesh':
                fig.cf1 = iplt.pcolormesh(layer.cube,
                                          cmap=layer.cmap,
                                          alpha=alpha)
            elif layer.plottype == 'contour':
                fig.cf1 = iplt.contour(layer.cube,
                                       cmap=layer.cmap,
                                       alpha=alpha)
            elif layer.plottype == 'contourf':
                fig.cf1 = iplt.contourf(layer.cube,
                                        cmap=layer.cmap,
                                        alpha=alpha)
            elif layer.plottype == 'contourf_edge':
                fig.cf1 = iplt.contourf(layer.cube,
                                        cmap=layer.cmap,
                                        alpha=alpha)
                iplt.contour(layer.cube, levels=fig.cf1.levels, colors='k')
            elif layer.plottype == 'scatter_latlon':
                fig.cf1 = iplt.scatter(layer.cube.coord('longitude'),
                                       layer.cube.coord('latitude'),
                                       c=np.atleast_1d(layer.cube.data),
                                       marker=layer.marker,
                                       s=layer.markersize,
                                       edgecolors="k",
                                       cmap=layer.cmap)

        else:

            if layer.levels[0] > np.max(layer.cube.data):
                # If max data point is less than lowest level don't
                # plot and don't add a colorbar
                layer.cbar = False
            elif layer.plottype == 'pcolormesh':
                fig.cf1 = iplt.pcolormesh(layer.cube,
                                          cmap=layer.cmap,
                                          norm=layer.norm,
                                          alpha=alpha)
            elif layer.plottype == 'contour':
                fig.cf1 = iplt.contour(layer.cube,
                                       levels=layer.levels,
                                       cmap=layer.cmap,
                                       norm=layer.norm,
                                       alpha=alpha)
            elif layer.plottype == 'contourf':
                fig.cf1 = iplt.contourf(layer.cube,
                                        levels=layer.levels,
                                        cmap=layer.cmap,
                                        norm=layer.norm,
                                        alpha=alpha)
            elif layer.plottype == 'contourf_edge':
                fig.cf1 = iplt.contourf(layer.cube,
                                        levels=layer.levels,
                                        cmap=layer.cmap,
                                        norm=layer.norm,
                                        alpha=alpha)
                iplt.contour(layer.cube,
                             levels=layer.levels,
                             colors='k')
            elif layer.plottype == 'scatter_latlon':
                fig.cf1 = iplt.scatter(layer.cube.coord('longitude'),
                                       layer.cube.coord('latitude'),
                                       c=np.atleast_1d(layer.cube.data),
                                       marker=layer.marker,
                                       s=layer.markersize,
                                       edgecolors="k",
                                       label=layer.label,
                                       cmap=layer.cmap,
                                       norm=layer.norm)

        # Set an axes extent (assuming extent is given in WGS84)
        # Needs to follow plot in order that axes have become geoaxes
        if self.extent is not None:
            plt.gca().set_extent(self.extent, crs=ccrs.PlateCarree())

        # Add a colorbar if required
        # First extract or determine required orientation:
        if layer.cbar_orientation:
            layer.cbar_orientation = layer.cbar_orientation
        elif layer.colorscale == 'linear':
            layer.cbar_orientation = 'horizontal'
        else:
            layer.cbar_orientation = 'vertical'

        # Please note that at this point, unless the 'layer.cbar' variable is
        # explicitly set as False, a colorbar will be constructed anyway.
        if layer.cbar is not False:
            layer.construct_cbar(fig.cf1,
                                 position=None,
                                 orientation=layer.cbar_orientation,
                                 title=None,
                                 title_fontsize=None,
                                 label_fontsize=None,
                                 tickmark_size=None)

        # Tell the object which fig to hold on to:
        self.fig = fig

    def plot_map(self):
        """
        If requested add some mapping background to the plot
        """

        ax = plt.gca()

        # Coordinates are transformed to lat-lon to
        # compute extent
        if self.extent is not None:
            extent = self.extent
        elif self.layer_list[0].cube_extent is not None:
            extent = self.layer_list[0].cube_extent
        else:
            extent = ax.get_extent(crs=ccrs.PlateCarree())

        lat_range = abs(extent[3] - extent[2])
        if lat_range < 10.0:
            res = '10m'
        elif lat_range < 50.0:
            res = '50m'
        else:
            res = '110m'

        if self.mapping == 'grey_os':
            if 0.0 < lat_range < 0.5:
                zoom_level = 11
            elif 0.5 <= lat_range < 1.0:
                zoom_level = 10
            elif 1.0 <= lat_range < 1.5:
                zoom_level = 9
            elif 1.5 <= lat_range < 4.0:
                zoom_level = 8
            elif lat_range >= 4.0:
                zoom_level = 7
            tiles = GreyscaleOSMMapTiles()
            ax.add_image(tiles, zoom_level)

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
            if self.rsmc:
                ax.coastlines(res, color='dimgray', zorder=3,
                              linewidth=1.5)
                ax.add_feature(countries, edgecolor='darkgray',
                               zorder=2, linewidth=1)
                ax.add_feature(states, edgecolor='darkgray',
                               linestyle='-.',
                               zorder=2)
            else:
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

        if self.mapping == 'jam':
            # for larger JAM maps - use coastlines only
            if (extent[1] - extent[0] > 60.0) or (extent[3] - extent[2] > 60.0):
                ax.coastlines(res)
            # use WMS mapping for smaller maps
            else:
                num_layers = np.linspace(0, 40, 41)
                layers = ['{:.0f}'.format(x) for x in num_layers]
                ax.add_wms(wms='http://exxdmmsprd01:6080/arcgis/' +
                           'services/DMMS/Global_NE_HC_Hybrid_Greyscale/' +
                           'MapServer/WMSServer', layers=layers)

        # Set-up the gridlines
        if self.gridlines:
            try:
                # for JAM map styles: use calculated gridlines
                # as required when setting font size to be small
                if self.mapping == 'jam':
                    xlocs, xlocs_extend = compute_grid_line_locs(
                        extent, maxbins=6)
                    ylocs, ylocs_extend = compute_grid_line_locs(
                        extent, lat_or_lon='lat', maxbins=6)
                    ax.gridlines(xlocs=xlocs_extend, ylocs=ylocs_extend)
                    gl = ax.gridlines(draw_labels=True,
                                      xlocs=xlocs, ylocs=ylocs,
                                      linewidth=0.001)
                    gl.xlabel_style = {'size': 8, 'color': 'black'}
                    gl.ylabel_style = {'size': 8, 'color': 'black'}
                elif extent[0] < 180 and extent[1] > 180:
                    xlocs, xlocs_extend = compute_grid_line_locs(extent)
                    ax.gridlines(xlocs=xlocs_extend)
                    gl = ax.gridlines(draw_labels=True,
                                      xlocs=xlocs,
                                      linewidth=0.001)
                elif (self.gridlines_xmaxbins is not None and
                      self.gridlines_ymaxbins is not None):
                    xlocs, xlocs_extend = compute_grid_line_locs(
                        extent, lat_or_lon='lon',
                        maxbins=self.gridlines_xmaxbins)
                    ylocs, ylocs_extend = compute_grid_line_locs(
                        extent, lat_or_lon='lat',
                        maxbins=self.gridlines_ymaxbins)
                    ax.gridlines(xlocs=xlocs_extend, ylocs=ylocs_extend)
                    gl = ax.gridlines(draw_labels=True,
                                      xlocs=xlocs, ylocs=ylocs,
                                      linewidth=0.8,
                                      alpha=0.9, zorder=9)
                else:
                    gl = ax.gridlines(draw_labels=True,
                                      linewidth=0.8,
                                      alpha=0.9, zorder=9)

                if self.gridlines_xlabel_size is not None:
                    gl.xlabel_style = {'size': self.gridlines_xlabel_size,
                                       'color': 'black'}

                if self.gridlines_ylabel_size is not None:
                    gl.ylabel_style = {'size': self.gridlines_ylabel_size,
                                       'color': 'black'}

                gl.xlabels_top = False
                gl.ylabels_right = False
                gl.xformatter = LONGITUDE_FORMATTER
                gl.yformatter = LATITUDE_FORMATTER
            except:
                gl = ax.gridlines()

    def gen_title(self):
        """
        Generates an automatic title based on validity time
        of the first layer in the list.

        If layer is a DAQI quantity, title only gives date.
        """
        t_coord = self.layer_list[0].cube.coord('time')
        t_stamp = t_coord.units.num2date(t_coord.points[0])
        short_name = self.layer_list[0].cube.attributes['short_name']
        # If plotting a DAQI quantity only give valid date (and not time)
        if short_name.find('DAQI') >= 0:
            title = 'Valid on {}'.format(t_stamp.strftime('%d/%m/%Y'))
        else:
            title = 'Valid at {}'.format(t_stamp.strftime('%H:%M %d/%m/%Y'))

        return title

    def gen_name_verbose_title(self):
        """
        Adds a more verbose title
        For NAME plots
        """

        ncube = self.layer_list[0].cube

        quantity = ncube.attributes['Quantity']
        if 'Time Av or Int' in ncube.attributes:
            ave_int = ncube.attributes['Time Av or Int']
        else:
            ave_int = ncube.attributes['Time av/int info']

        if quantity == 'Dosage':
            quantity = 'Air concentration'

        height_text = height_text_fmt(self.layer_list[0].cube)

        # Compute time bounds
        t_coord = ncube.coord('time')
        t_stamp1 = t_coord.units.num2date(t_coord.bounds[0][0])
        t_stamp2 = t_coord.units.num2date(t_coord.bounds[0][1])

        title = '{}\n'.format(quantity)
        title += '{}\n'.format(ave_int.title())
        if quantity not in ['Deposition', 'Total deposition',
                            'Wet deposition', 'Dry deposition',
                            'Dry Deposition', 'Wet Deposition']:
            title += '{}\n'.format(height_text)
        t_fmt = '%H:%M %d/%m/%Y'
        # The NAME output files are in UTC therefore the title can have UTC
        # hard coded in the title.
        title += 'Valid from {} to {} UTC'.format(t_stamp1.strftime(t_fmt),
                                                  t_stamp2.strftime(t_fmt))

        return title

    def gen_filename(self, filesuffix='', fileprefix='Fieldplot_'):
        """
        Generate automatic filename (doesn't include path)
        Extracts short name from layer 1 and date stamp
        if available.

        :param filesuffix: String to append to filename, prior to file
                           extension.
        """

        filename = fileprefix
        label = self.layer_list[0].cube.attributes['label']
        filename += label + '_'
        short_name = self.layer_list[0].cube.attributes['short_name']
        short_name = short_name.replace('.',
                                        '_')  # short name should not contain a decimal point.
        filename += short_name

        if 'NAME Version' in self.layer_list[0].cube.attributes:
            height = height_text_fmt(self.layer_list[0].cube,
                                     text_type='filename')
            filename += '_' + height

        try:
            t_coord = self.layer_list[0].cube.coord('time')
            t_stamp = t_coord.units.num2date(t_coord.points[0])
            filename += t_stamp.strftime('_%Y%m%d%H%M')
        except:
            filename += ''

        filename += filesuffix
        filename += '.' + self.figformat

        return filename

    def save_fig(self, plotdir=None, filename=None, filesuffix='',
                 fileprefix='Fieldplot_', tight=False, verbose=1):
        """
        Save figure to file.

        :param plotdir: Directory to save plot in. If not set,
                        will use plotdir set previously,
                        or default to current directory.
        :param filename: Filename to save plot as.
                         If filename not given, will generate automatically.
        :param filesuffix: String to append to default filename, prior to file
                       extension. Ignored if filename is not None
        :param tight: If set to True, will adjust figure size to
                      be tight to bounding box of figure and minimise
                      whitespace. Note this will change size of output file
                      and may not be consistent amongst different plots
                      (particularly if levels on colorbar are different).
        :param verbose: 0 (no extra print output) or
                        1 (standard level of print output)
                        2 (extra level of print output for debugging)

        """

        if plotdir is None:
            if self.plotdir is None:
                self.plotdir = './'
            plotdir = self.plotdir
        # Check finishes with /
        if plotdir[-1] != '/':
            plotdir += '/'

        if filename is None:
            filename = self.gen_filename(filesuffix=filesuffix,
                                         fileprefix=fileprefix)

        if not os.path.isdir(plotdir):
            if verbose > 1:
                print('Creating output directory:', plotdir)
            # When we do not need Python 2 anymore this one line
            # can replace the two below
            # os.makedirs(plotdir, exist_ok=True)
            if not os.path.exists(plotdir):
                os.makedirs(plotdir)

        if tight:
            self.fig.savefig(plotdir + filename, bbox_inches='tight')
        else:
            self.fig.savefig(plotdir + filename)
        if verbose >= 1:
            print('Saved figure ', plotdir + filename)
        plt.close()  # Clear current figure to avoid opening too many at once


if __name__ == '__main__':
    import doctest

    doctest.testmod()
