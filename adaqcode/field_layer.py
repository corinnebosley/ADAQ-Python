"""
field_layer.py

Contains all the code for constructing, querying and defining a field layer
which can then be plotted using the field_plot class
"""
from __future__ import division

from matplotlib import pyplot as plt
from six.moves.builtins import zip
from six.moves.builtins import str
from six.moves.builtins import range
from six.moves.builtins import object

import warnings

import iris
import numpy as np
import matplotlib.colors as mplcolors
import matplotlib.pyplot as plt

#------------------------------------------------------------------
# Some additional modules which are not part of the class
#--------------------------------------------------------
from adaqcode.plotting_functions import units_str


def _everything_but_coords(cube, everything_but_coord_name):
    """
    Quick script which returns the names of all coordinates in a cube
    except those in the list provided
    """

    keep_dims = []
    for item in everything_but_coord_name:
        keep_dims += cube.coord_dims(cube.coord(item))
    # Return the dimension coordinates which are not on the given
    # coordinate's dimensions.
    return [coord for coord in cube.coords(dim_coords=True)
            if not any(dim in keep_dims for dim in cube.coord_dims(coord))]

def _expand_mask(array):
    """
    Used in mask_data function for 'contour' plots to expand
    mask by one to ensure smooth outer contour
    """

    hgt, wdt = array.shape
    mi = np.where(not array.mask)
    for j, i in zip(mi[0], mi[1]):
        i0 = max(0, i-1)
        j0 = max(0, j-1)
        i1 = min(wdt-1, i+1)
        j1 = min(hgt-1, j+1)
        array.mask[j0:j1+1, i0:i1+1] = False


def _get_extent_ge(cube, thresh, border=2.0):
    """
    Used in setup_extent function to return the (x0,x1,y0,y1) range
    where the cube exceeds the lowest contour (thresh).
    """

    if len(cube.shape) > 2:
        cube = cube.collapsed(_everything_but_coords(cube,
                                                     ['latitude', 'longitude']),
                              iris.analysis.MAX)

    x_dim = cube.shape[-1]
    y_dim = cube.shape[-2]

    # find out which columns and rows contain useful values
    x_ge = np.zeros((x_dim,))
    y_ge = np.zeros((y_dim,))

    for x in range(x_dim):
        if np.any(cube.data[:, x] >= thresh):
            x_ge[x] = True
    for y in range(y_dim):
        if np.any(cube.data[y, :] >= thresh):
            y_ge[y] = True

    cols = x_ge.nonzero()[0]
    rows = y_ge.nonzero()[0]

    if cols.size == 0 or rows.size == 0:
        warnings.warn("No data points >= threshold ("+str(thresh)+"), " +
                      "a cube extent cannot be calculated, " +
                      "defaulting to global extent", RuntimeWarning)
        crs = cube.coord(axis="X").coord_system.as_cartopy_crs()
        extent = [-179.9, 179.9, -89.9, 89.9] # Taken in by 0.1 deg due to
        return (extent, crs)                  # Cartopy bug - Dec 2017

    # get the tight range of interest
    x_min = cols.min()
    x_max = cols.max()
    x_centre = (x_min + x_max) / 2
    x_range = x_max - x_min

    y_min = rows.min()
    y_max = rows.max()
    y_centre = (y_min + y_max) / 2
    y_range = y_max - y_min

    # Assuming roughly square cell sizes, add a border (limited to cube range).
    dist_from_cen = max(x_range, y_range) * border
    dist_from_cen = np.ceil(dist_from_cen)
    x_min = int(max(0, x_centre - dist_from_cen))
    x_max = int(min(x_dim - 1, x_centre + dist_from_cen))
    y_min = int(max(0, y_centre - dist_from_cen))
    y_max = int(min(y_dim - 1, y_centre + dist_from_cen))

    # turn data indices into native coord values
    x_pts = cube.coord(axis="X").points
    y_pts = cube.coord(axis="Y").points
    extent = [x_pts[x_min], x_pts[x_max], y_pts[y_min], y_pts[y_max]]

    crs = cube.coord(axis="X").coord_system.as_cartopy_crs()
    return (extent, crs)


#---------------------------------------------------------------------------
# Field Layer Object
#---------------------------------------------------------------------------

class FieldLayer(object):

    """
    Object for holding information about a field layer

    As an example, load in some model data and compute suitable levels for
    displaying the layer.

    First import name_data and config

    >>> import name_data
    >>> import config

    Use the sample data path to locate data

    >>> sample_data_path = config.SAMPLE_DATADIR+'name/'

    Read in the data

    >>> name = name_data.NAMEData()
    >>> name.readdata(sample_data_path + 'Fields_grid2_201304172000.txt',
    ... field_attributes = {'Quantity': 'Wet deposition'})
    [<iris 'Cube' of TRACER_WET_DEPOSITION / (g/m2) \
(latitude: 200; longitude: 200)>]
    >>> mod_cube = name.gridded_cube_list[0]
    >>> print(mod_cube)
    TRACER_WET_DEPOSITION / (g/m2)      (latitude: 200; longitude: 200)
         Dimension coordinates:
              latitude                           x               -
              longitude                          -               x
         Scalar coordinates:
              source_latitude: 52.0394 degrees
              source_longitude: -2.6868 degrees
              time: 2013-04-17 20:00:00, bound=\
(2013-04-17 14:00:00, 2013-04-17 20:00:00)
              z: Boundary layer
         Attributes:
              End of release: 2000UTC 17/04/2013
              Forecast duration: 6 hours
              Met data: NWP Flow.UKV_PT1_flow; NWP Flow.UKV_PT2_flow
              NAME Version: NAME III (version 6.1)
              Quantity: Wet deposition
              Release height: 0.000 to 1.000m agl
              Release location: 2.6868W   52.0394N
              Release rate: 1.000000g/s
              Run time: 1245UTC 17/04/2013
              Species: TRACER
              Species Category: CHEMISTRY-SPECIES
              Start of release: 1400UTC 17/04/2013
              Time Av or Int: 006 hr time integrated
              Title: USER_DEFINED_NG_17042013_1245Z
              label: NAME
              short_name: TRACER_WET_DEPOSITION
         Cell methods:
              sum: time

    Set up field_layer class and add cube to it

    >>> fl = FieldLayer ( mod_cube )

    Set the layer style

    >>> fl.set_layerstyle(plottype='pcolormesh', colorscale='log',
    ... mask=True, autozoom=True)

    Print out the levels to check they are sensible

    >>> np.set_printoptions(formatter={'float':lambda x: '{:5.2e}'.format(x)})
    >>> print(fl.levels)
    [3.16e-09 1.00e-08 3.16e-08 1.00e-07 3.16e-07 1.00e-06 3.16e-06 1.00e-05
     3.16e-05 1.00e-04]

    Print out the field layer object to see everything set up:

    >>> print(fl) # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    <class '...FieldLayer'>
    autozoom: True
    cbar: None
    cbar_label: default
    cbar_num_fmt: None
    cbar_orientation: None
    cmap: YlGnBu (<matplotlib.colors.LinearSegmentedColormap object at ...>)
    colors: None
    colorscale: log
    cube: TRACER_WET_DEPOSITION / (g/m2)      (latitude: 200; longitude: 200)
         Dimension coordinates:
              latitude                           x               -
              longitude                          -               x
         Scalar coordinates:
              source_latitude: 52.0394 degrees
              source_longitude: -2.6868 degrees
              time: 2013-04-17 20:00:00, bound=\
(2013-04-17 14:00:00, 2013-04-17 20:00:00)
              z: Boundary layer
         Attributes:
              End of release: 2000UTC 17/04/2013
              Forecast duration: 6 hours
              Met data: NWP Flow.UKV_PT1_flow; NWP Flow.UKV_PT2_flow
              NAME Version: NAME III (version 6.1)
              Quantity: Wet deposition
              Release height: 0.000 to 1.000m agl
              Release location: 2.6868W   52.0394N
              Release rate: 1.000000g/s
              Run time: 1245UTC 17/04/2013
              Species: TRACER
              Species Category: CHEMISTRY-SPECIES
              Start of release: 1400UTC 17/04/2013
              Time Av or Int: 006 hr time integrated
              Title: USER_DEFINED_NG_17042013_1245Z
              label: NAME
              short_name: TRACER_WET_DEPOSITION
         Cell methods:
              sum: time
    cube_crs: <cartopy._crs.Geodetic object at ...>
    cube_extent: [-2.752..., -2.608..., 52.029..., 52.118...]
    label: None
    levels: [3.16e-09 1.00e-08 3.16e-08 1.00e-07 3.16e-07 1.00e-06 3.16e-06 1.00e-05
     3.16e-05 1.00e-04]
    marker: o
    markersize: 20
    mask: True
    nlevels: 10
    norm: <matplotlib.colors.BoundaryNorm object at ...>
    plottype: pcolormesh
    step: 0.5
    <BLANKLINE>
    >>> np.set_printoptions()

    """

    def __init__(self, cube=None):

        """
        Initialisation of the FieldLayer Class
        """

        # The cube
        self.cube = cube

        # Plottype
        self.plottype = 'pcolormesh'
        self.mask = False

        #For scatter plots
        self.marker = 'o'
        self.markersize = 20
        self.label = None

        # Levels for contours or colormapping
        self.nlevels = 10
        self.step = 0.5
        self.levels = None
        self.colorscale = 'linear'

        # Colormap and colormap normalisation
        self.cmap = 'YlGnBu'
        self.colors = None
        self.norm = None

        # Colorbar
        self.cbar = None
        self.cbar_orientation = None
        self.cbar_label = 'default'
        self.cbar_num_fmt = None

        # Extent information
        self.autozoom = False
        self.cube_extent = None
        self.cube_crs = None

        # If cube is not None check it is a cube
        if self.cube is not None:
            if not isinstance(self.cube, iris.cube.Cube):
                raise ValueError("cube is not an Iris Cube")


    def __str__(self):
        """
        Output from 'print'
        """
        string = str(type(self))+'\n'
        keys = sorted(self.__dict__.keys())
        for key in keys:
            value = self.__dict__[key]
            if key == 'cmap':
                value = value.name + ' ('+str(value)+')'
            string += key + ': ' + str(value) + '\n'

        return string


    def set_layerstyle(self, nlevels=10, levels=None, step=0.5,
                       plottype='pcolormesh', mask=False, autozoom=False,
                       colorscale='linear', cmap='YlGnBu', colors=None):

        """
        Method to allow various plotting styles to be set in a single
        command.

        :param nlevels: number of contour levels to use in plotting
        :param levels: a list of contour levels for plotting
        :param plottype: contour, contourf (filled contour) or pcolormesh
                         (pixel) or scatter_latlon (scatter plot - cube
                         must have latitude and longitude coordinates)
        :param mask: mask out values below lowest contour (True or False)
        :param autozoom: zoom in to plume extent (True or False)
        :param colorscale: log or linear colorscale
        :param cmap: colormap to use
        :param colors: a list of colors to use in place of a colormap

        .. note::
            colors should either be a hex array or an array of RGB values in
            range 0-1

        .. note::
            self.nlevels will be ignored if self.levels is not None

        """

        if nlevels is not None:
            self.nlevels = nlevels
        if step is not None:
            self.step = step
        if colorscale is not None:
            self.colorscale = colorscale
        if levels is not None:
            self.levels = levels
        if cmap is not None:
            self.cmap = cmap
        if colors is not None:
            self.colors = colors
        if plottype is not None:
            self.plottype = plottype
        if mask is not None:
            self.mask = mask
        if autozoom is not None:
            self.autozoom = autozoom

        # Set levels
        if self.levels is None and self.colorscale == 'log':
            self.get_loglevels()
        if self.levels is None and self.colorscale == 'linear':
            self.get_linearlevels()

        # Then obtain a colormap and normalisation
        if self.levels is not None:
            self.get_cmap_and_norm()

        # Apply mask if required
        if self.mask:
            self.mask_data()

        # Get plume extent if autozoom is True and coordinates
        # include latitude and longitude
        if self.autozoom:
            coord_names = [coord.name() for coord in self.cube.coords()]
            if 'latitude' not in coord_names and 'longitude' not in coord_names:
                warnings.warn('Autozoom currently only available for ' +
                              'latitude and longitude coordinates')
            else:
                self.setup_extent()

    def get_linearlevels(self):
        """
        Determine levels for a linear scale

        Uses:
          * self.cube - single IRIS cube (may be multi-dimensional)
          * self.nlevels: number of levels (optional defaults to 10)

        Produces:
          * self.levels - levels for plotting

        As an example, load in some model data and compute linear levels
        Use the sample data path to locate data

        >>> import config
        >>> sample_data_path = config.SAMPLE_DATADIR+'name/'
        >>> cube = iris.load_cube(sample_data_path +
        ... 'Fields_grid99_201105230000.txt',
        ...  'VOLCANIC_ASH_TOTAL_DEPOSITION')

        Set up a field_layer and add the cube

        >>> fl = FieldLayer(cube=cube)
        >>> fl.nlevels = 4
        >>> np.set_printoptions(formatter={'float':lambda x: '{:5.2f}'.format(x)})
        >>> fl.get_linearlevels()
        array([ 0.00, 11.42, 22.84, 34.26])
         >>> np.set_printoptions()
        """

        data_min = self.cube.data.min()
        data_max = self.cube.data.max()

        self.levels = np.linspace(data_min, data_max, self.nlevels)

        return self.levels


    def get_loglevels(self):
        """
        Determine levels for a log scale

        Uses:
          * self.cube - single IRIS cube (may be multi-dimensional)
          * self.nlevels: number of levels (optional defaults to 10)
          * self.step: step between contours(defaults to 0.5 (or half-log))

        Produces:
          * self.levels - levels for plotting

        As an example, load in some model data and compute log levels
        Use the sample data path to locate data

        >>> import config
        >>> sample_data_path = config.SAMPLE_DATADIR+'name/'
        >>> cube = iris.load_cube(sample_data_path +
        ... 'Fields_grid99_201105230000.txt',
        ...  'VOLCANIC_ASH_TOTAL_DEPOSITION')

        Set up a field_layer and add the cube

        >>> fl = FieldLayer(cube=cube)
        >>> fl.nlevels = 4
        >>> np.set_printoptions(formatter={'float':lambda x: '{:5.2f}'.format(x)})
        >>> fl.get_loglevels()
        array([ 3.16, 10.00, 31.62, 100.00])
        >>> fl.step = 1
        >>> fl.get_loglevels()
        array([ 0.10,  1.00, 10.00, 100.00])
        >>> np.set_printoptions()
        """

        data_max = self.cube.data.max()

        maxval = np.ceil(np.log10(data_max))

        if np.isfinite(maxval):
            self.levels = np.logspace(maxval - self.nlevels*self.step
                                      + self.step, maxval, num=self.nlevels)
        else:
            self.levels = [1.0, 10.0, 100.0]

        return self.levels


    def get_cmap_and_norm(self):
        """
        Given a set of colors build a colormap and a normaliser

        Uses:
          * self.levels - levels for color boundaries
          * self.cmap - colormap name
          * self.colors (optional) - array of colors for constructing own
            colormap

        Produces:
          * self.cmap - colormap
          * self.norm - normaliser

        >>> fl = FieldLayer()
        >>> fl.levels = [1., 2., 5., 10., 20., 50.]
        >>> cmap, norm = fl.get_cmap_and_norm()
        >>> print(norm)  # doctest: +ELLIPSIS
        <matplotlib.colors.BoundaryNorm object at ...>
        """

        if self.colors is not None:
            self.cmap = mplcolors.ListedColormap(self.colors, name=self.cmap)
            self.norm = mplcolors.BoundaryNorm(self.levels, len(self.colors))
        else:
            if isinstance(self.cmap, str):
                #Convert to valid colour map object
                self.cmap = plt.get_cmap(self.cmap)
            self.norm = mplcolors.BoundaryNorm(self.levels, self.cmap.N)

        return self.cmap, self.norm

    def mask_data(self):
        """
        Mask data below the lowest contour level.
        Also masks any invalid (nan) data.
        If contour plot specified then the mask is expanded by one to
        ensure smoother outer contour.
        """

        data_field1 = self.cube.data
        data_field_finite = np.ma.masked_invalid(data_field1)

        if self.levels is not None:
            data_field = np.ma.masked_less(data_field_finite, self.levels[0])
        else:
            warnings.warn("Data cannot be masked below lowest contour level \
                                as field_layer levels not set")
            data_field = data_field_finite

        if self.plottype == 'contour':
            _expand_mask(data_field)
            self.cube.data = data_field
        else:
            self.cube.data = data_field

    def setup_extent(self):
        """
        Extra method for use in determining an extent which can then be passed
        to FieldPlot. Returns an extent and a reference projection

        .. note::
           This still requires some work to cope with the wide range of
           files and plume extents
        """

        if self.levels is None:
            raise ValueError("Cannot determine extent as no levels set")

        cube_extent, cube_crs = _get_extent_ge(self.cube,
                                               self.levels[0],
                                               border=0.6)

        # This bit of code checks to see whether a smaller extent could be used
        # if a -180:180 longitude range is used instead of a 0:360 longitude
        # range or vice versa
        if abs(cube_extent[1] - cube_extent[0]) > 330:
            lon_points = self.cube.coord('longitude').points
            if min(lon_points) < -10:
                subset = self.cube.intersection(longitude=(0, 359))
            else:
                subset = self.cube.intersection(longitude=(-179, 179))

            cube_extent_2, cube_crs_2 = _get_extent_ge(subset,
                                                       self.levels[0],
                                                       border=0.6)

            if abs(cube_extent_2[1] - cube_extent_2[0]) <= 330:
                self.cube = subset
                cube_extent = cube_extent_2
                cube_crs = cube_crs_2

        # If longitude extent is close to global then extend to global
        # and centre on the Greenwich Meridian due to British interest
        # even if the extent was smaller when centred on the Date Line.
        #
        # Bug in Cartopy Dec'17 - Global extent will not be plotted with
        # extent[0] = -180, extent[1] = 180 So the extents are deliberately
        # taken in by 0.1

        if cube_extent[1] - cube_extent[0] > 330:
            cube_extent[0] = -179.9
            cube_extent[1] = 179.9

        self.cube_extent = cube_extent
        self.cube_crs = cube_crs

        return self.cube_extent, self.cube_crs

    def layer_slice(self, dim_slices):
        """
        Method for slicing the layer into 2D arrays which can
        be plotted. In reality this means slicing the cube bit and
        keeping all the rest of the information with the slice.
        """

        if not isinstance(dim_slices, list):
            dim_slices = [dim_slices]

        for cube_slice in self.cube.slices(dim_slices):
            layer_slice = self
            layer_slice.cube = cube_slice
            yield layer_slice

    def layer_slice_over(self, dim_slices):
        """
        Method for slicing the layer over the specified dimensions.
        In reality this means slicing the cube bit and keeping
        all the rest of the information with the slice.
        """

        if not isinstance(dim_slices, list):
            dim_slices = [dim_slices]

        for cube_slice in self.cube.slices_over(dim_slices):
            layer_slice = self
            layer_slice.cube = cube_slice
            yield layer_slice

    def construct_cbar(self, mappable,
                       position=None,
                       orientation=None,
                       title=None,
                       title_fontsize=12,
                       label_fontsize=10,
                       tickmark_size=8):
        """
        Creates a standard colourbar for use within FieldPlots when plotting
        multi-layer datasets i.e. FieldPlot.plot_layer() or
        FieldPlot.multiplot().


        :param mappable: plot for which to construct colorbar.

        :param position: position of colorbar axes within figure.

        :param orientation: orientation of colorbar, i.e. 'vertical' or
                            'horizontal'.

        :param title: title of phenomenon to plot.

        :param label_fontsize: fontsize for colorbar label (defaults to 10).

        :param title_fontsize: fontsize for colorbar title (defaults to 12).

        :param tickmark_size: fontsize for colorbar ticks (defaults to 8).

        """
        # Determine required properties for colorbar construction.

        # Orientation:
        if orientation is not None:
            self.cbar_orientation = orientation

        # Number format (look for a requested format, otherwise try and choose a
        # sensible one):
        if self.cbar_num_fmt is not None:
            cbar_num_fmt = self.cbar_num_fmt
        elif self.colorscale != 'linear':
            cbar_num_fmt = '%.1e'
        else:
            cbar_num_fmt = None

        # Labelling:
        if self.cbar_label == 'default':
            cbar_label = self.cube.attributes['short_name']
        else:
            cbar_label = self.cbar_label
        # Append units to label
        units = units_str(self.cube.units)
        if units != '1':
            cbar_label = '{} [{}]'.format(cbar_label, units)
        else:
            cbar_label = '{}'.format(cbar_label)

        # Make colorbar using properties defined or determined and in the
        # position either specified or calculated:
        if position is not None:
            cbar = plt.colorbar(mappable,
                                cmap=self.cmap,
                                norm=self.norm,
                                orientation=self.cbar_orientation,
                                format=cbar_num_fmt,
                                ticks=self.levels,
                                cax=position)
        else:
            if self.cbar_orientation == 'vertical':
                fraction = 0.042
            elif self.cbar_orientation == 'horizontal':
                fraction = 0.08
            cbar = plt.colorbar(mappable,
                                cmap=self.cmap,
                                norm=self.norm,
                                orientation=self.cbar_orientation,
                                format=cbar_num_fmt,
                                ticks=self.levels,
                                fraction=fraction,
                                pad=0.1)
        cbar.ax.tick_params(labelsize=tickmark_size)
        cbar.set_label(cbar_label, fontsize=label_fontsize)
        if title is not None:
            cbar.ax.set_title(title, fontsize=title_fontsize)


if __name__ == '__main__':

    import doctest
    doctest.testmod()

