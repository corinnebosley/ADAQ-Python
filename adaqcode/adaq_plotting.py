"""
Plotting of ADAQData objects.
These routines all generally take the following arguments:

 * **ini_dict -** Dictionary of a :class:`inifile` object.
   Must contain 'plot_dir'.
 * **od -** observation data in the form of an
   :class:`adaq_data.ADAQData` object.
 * **md_list -** list of model data objects in the form of an
   :class:`adaq_data.ADAQData` objects.
 * **short_name -** string to match to short_name attribute in cubes

Plus some other routine-specific arguments.
"""
from __future__ import division
from __future__ import print_function

import os
import warnings

import iris
import iris.iterate
import matplotlib.pyplot as plt
import numpy as np

# adaqcode
import cube_statistics
import cube_time
import plotting_functions
from field_layer import FieldLayer
from field_plot import FieldPlot
from statistics_plotting import Histogram, QQPlot, SoccerPlot
from timeseries_plot import TimeSeriesPlot, tsp_statistic

import timeseries_stats
from trajectory_plot import (compute_trajectory_extent,
                             compute_traj_direction, TrajectoryPlot)

Z_COORD_YLABEL = {'flight_level': 'FL',
                  'altitude': 'm asl',
                  'height': 'm agl',
                  'Z (Pa)': 'hPa'
                  }


def plot_diurnal(ini_dict, od, md_list, short_name, filesuffix=''):
    """
    Plot diurnal mean plots for a given short_name.

    :param ini_dict: Dictionary of an :class:`inifile` object.
                     Must contain 'plot_dir'.
                     May contain 'line_colours_list'
                     (default: :any:`plotting_functions.COLOURS`).
    :param od: observation data in the form of an :class:`adaq_data.ADAQData`
               object
    :param md_list: list of model data objects in the form of an
                    :class:`adaq_data.ADAQData` objects.
    :param short_name: string to match to short_name attribute in cubes.
    :param filesuffix: String to add to end of filename for specific naming
                       purposes. By default adds nothing.

    :return: :class:`timeseries_plot.TimeSeriesPlot` object.

    >>> import adaq_functions
    >>> ini_dict, sites_data, od, md_list = adaq_functions.get_exampledata(
    ... exampletype="full") # doctest: +ELLIPSIS
    Reading inifile ...example_data_5days.ini
    Number of sites:  5

    For this example only, manually change plot directory to save to gallery:

    >>> import config
    >>> ini_dict['plot_dir'] = (config.CODE_DIR +
    ...     "/adaqdocs/figures/adaq_plotting")
    >>> dp = plot_diurnal(ini_dict, od, md_list, 'O3') # doctest: +ELLIPSIS
    Plotting diurnal
    Saved figure  .../figures/adaq_plotting/Diurnal_O3.png

    .. image:: ../adaqdocs/figures/adaq_plotting/Diurnal_O3.png
       :scale: 50%

    Can also check the 24 mean values from the observations (the first cube in
    the lines dictionary):

    >>> np.set_printoptions(formatter={'float':lambda x: '{:5.2f}'.format(x)})
    >>> print(dp.lines[0]['cube'].data)
    [60.20 55.44 55.56 55.52 54.16 52.20 51.80 50.16 49.36 53.44 58.48 64.88
     69.84 74.84 78.80 81.60 82.60 82.56 81.36 76.24 73.00 71.88 70.76 67.64]
    >>> np.set_printoptions()
    """

    print('Plotting diurnal')

    colours = ini_dict.get('line_colours_list', plotting_functions.COLOURS)

    dp = TimeSeriesPlot()
    dp.xlim = (0, 23)  # Manually set xlimit
    dp.xticks = [0, 6, 12, 18]  # Force markers to be every 6 hours
    dp.xlabel = 'Hour (Z)'

    startdt, enddt = None, None

    obs_cube = od.extract(short_name=short_name, singlecube=True)

    # Plot obs
    if obs_cube is not None:
        startdt, enddt = cube_time.get_startenddt(obs_cube)
        obs_cube_diurnal = cube_statistics.diurnal(obs_cube, collapsed=True)
        if np.sum(np.isfinite(obs_cube_diurnal.data)) > 0:
            dp.add_line(obs_cube_diurnal, x=obs_cube_diurnal.coord('hour'),
                        colour=colours[0], linestyle='--')

    # Plot models
    for imd, md in enumerate(md_list):

        mod_cube = md.extract(short_name=short_name, singlecube=True)

        if mod_cube is not None:

            if startdt is None:
                # Start date/time not set by obs cube, so set from first
                # available model cube
                startdt, enddt = cube_time.get_startenddt(mod_cube)

            mod_cube_diurnal = cube_statistics.diurnal(mod_cube,
                                                       collapsed=True)

            if np.sum(np.isfinite(mod_cube_diurnal.data)) > 0:
                dp.add_line(mod_cube_diurnal, x=mod_cube_diurnal.coord('hour'),
                            colour=colours[imd + 1])

    if not dp.lines:
        warnings.warn("No diurnal plot for " + short_name + "(No obs/models)")
        return dp

    dp.title = 'Diurnal Mean Cycle '
    startstr = startdt.strftime("%d/%m/%Y %H:%M")
    endstr = enddt.strftime("%d/%m/%Y %H:%M")
    dp.title += '(' + startstr + ' to ' + endstr + ')'
    dp.title += '\n ' + dp.lines[0]['cube'].name().replace('_', ' ')
    dp.plot()
    dp.save_fig(plotdir=ini_dict['plot_dir'],
                filename='Diurnal_' + short_name + filesuffix + '.png')

    return dp


def plot_histogram(ini_dict, od, md_list, short_name, filesuffix='',
                   binsize=None, maxperc=None):
    """
    Plot histogram for given short_name.

    :param ini_dict: Dictionary of an :class:`inifile` object.
                     Must contain 'plot_dir'.
                     May contain 'line_colours_list'
                     (default: :any:`plotting_functions.COLOURS`).
    :param od: observation data in the form of
               an :class:`adaq_data.ADAQData` object
    :param md_list: list of model data objects in the form of an
                    :class:`adaq_data.ADAQData` objects.
    :param short_name: string to match to short_name attribute in cubes.
    :param filesuffix: String to add to end of filename for specific naming
                       purposes. By default adds nothing.
    :param binsize: Binsize for histogram (uses default value if not set).
    :param maxperc: Maximum percentile from all data to display on x-axis.
                    For example setting to 99.5 will cut off the extreme
                    maximum values and therefore allow plot to zoom in
                    slightly.
    :return: :class:`statistics_plotting.Histogram` object.

    >>> import adaq_functions
    >>> ini_dict, sites_data, od, md_list = adaq_functions.get_exampledata()
    ... # doctest: +ELLIPSIS
    Reading inifile ...example_data_1days.ini
    Number of sites:  5

    For this example only, manually change plot directory to save to gallery:

    >>> import config
    >>> ini_dict['plot_dir'] = (config.CODE_DIR +
    ...     "/adaqdocs/figures/adaq_plotting")

    Produce histogram:

    >>> hist = plot_histogram(ini_dict,od,md_list,'O3') # doctest: +ELLIPSIS
    Plotting histogram
    Saved figure  .../figures/adaq_plotting/Histogram_O3.png
    >>> print(hist.alldata.min(), hist.alldata.max())
    0.0 96.0

    Example output file:

    .. image:: ../adaqdocs/figures/adaq_plotting/Histogram_O3.png
       :scale: 50%
    """

    print('Plotting histogram')

    colours = ini_dict.get('line_colours_list', plotting_functions.COLOURS)

    hist = Histogram()
    if binsize is not None:
        hist.binsize = binsize
    if maxperc is not None:
        hist.maxperc = maxperc

    # Add obs in black
    obs_cube = od.extract(short_name=short_name, singlecube=True)
    if obs_cube is not None:
        if np.sum(np.isfinite(obs_cube.data)) > 0:
            hist.add_line(obs_cube, colour=colours[0])

    # Add models
    for imd, md in enumerate(md_list):
        mod_cube = md.extract(short_name=short_name, singlecube=True)
        if mod_cube is not None:
            if np.sum(np.isfinite(mod_cube.data)) > 0:
                hist.add_line(mod_cube,
                              colour=colours[imd + 1])

    if not hist.lines:
        warnings.warn("No histogram plot for " + short_name +
                      "(No obs/models)")
        return hist

    hist.plot()
    hist.save_fig(plotdir=ini_dict['plot_dir'],
                  filename='Histogram_' + short_name + filesuffix + '.png')

    return hist


def plot_md_gridded_fields(ini_dict, md_list, short_name,
                           griddedcube=None, defaults=None,
                           filesuffix='', plot_subdir='',
                           titleprefix=None, tight=False,
                           verbose=1):
    """
    Plot pcolormesh fields automatically for all ADAQData objects in md_list.
    This assumes various default settings.
    Only does X (eg longitude) versus Y (eg latitude) axis plots.
    For more detailed adjustments and plotting, do not use this routine.

    :param ini_dict: Dictionary of an :class:`inifile` object. Could contain:

                      * extent_list - list of [xmin,xmax,ymin,ymax] in lat/lon
                        coords. Defaults to None if not in ini_dict.
                      * levels_list - list of contouring levels.
                        Defaults to None if not in ini_dict.
                      * levels_dict - dictionary whose keys are short_names
                        and whose values are the list of contouring levels
                        for that short_name. If required short_name is available
                        in dictionary then values from this dictionary are
                        used instead of those from levels_list. Defaults
                        to empty dictionary if in not in ini_dict.
                      * nlevels - Number of contour levels if levels_list
                        not set. Defaults to 10 if not in ini_dict.
                      * cmap - Colour map. Defaults to 'YlGnBu' if not
                        in ini_dict.
                      * cbar_orientation - Orientation of color bar
                        ('horizontal' or 'vertical').
                        Defaults to None if not in ini_dict which will add
                        color bar in default position according to levels type
                        (log or linear).
                      * cbar_label - label to add to color bar. If 'default'
                        will use the short name.
                      * cbar_num_fmt - Colourbar number format. Default depends
                        on if an AQ or NAME plot
                      * cbar - Set to True (default) to plot colour bar, set to
                        False to not include a colour bar
                      * annote_location - location of annotation in
                        relation to plot ('right' - next to plot,
                        'below' - below plot).
                      * annote - text to add to plot. Text will added to
                        right hand third of page or bottom third of page
                        depending on annote_location
                      * plot_dir - Top level output directory for plots.
                      * title - String to specify title
                      * projection - string containing map projection. Default
                        is 'PlateCarree'. If keyword defaults='AQ' then the
                        default projection is the model native grid (eg rotated
                        pole).

    :param md_list: list of model data objects in the form of an
                    :class:`adaq_data.ADAQData` objects.
    :param short_name: string to match to short_name attribute in cubes.
    :param griddedcube: iris Cube. Overrides md_list and plots gridded cube
                         instead of plotting cubes extracted from
                         md_list.gridded_cube_list.
    :param defaults: Use predefined default setups. Options available are:

                      * None - no default setups.
                      * 'NAME' - default setups for NAME plotting.
                        Adds coastlines, gridlines, central longitude and
                        PlateCarree projection.
                      * 'AQ' - default setups for Air Quality (AQUM etc)
                        plotting. Just adds coastlines and gridlines.
                        Uses native cube projections.

    :param filesuffix: String to add to end of filename for specific naming
                       purposes. By default adds nothing.
    :param plot_subdir: String to add to end of plotdir defined within ini_dict
                        to give add a subdirectory to the output plot location.
    :param titleprefix: String to prefix main title with.
                        If prefix ends with ``\\n`` this will put it
                        on the line above the automatic title.
    :param tight: If set to True, will adjust figure size to
                  be tight to bounding box of figure and minimise
                  whitespace. Note this will change size of output file
                  and may not be consistent amongst different plots
                  (particularly if levels on colorbar are different).
    :param verbose: Level of print output:

                    * 0 = No extra printing
                    * 1 = Standard level (recommended)
                    * 2 = Extra level for debugging


    >>> import adaq_functions
    >>> import config
    >>> ini_dict, sites_data, od, md_list = adaq_functions.get_exampledata(
    ... gridded_cube_list=True, sites_cube_list=False) # doctest: +ELLIPSIS
    Reading inifile .../example_data_1days.ini

    For this example only, ensure no defaults are set from ini_dict,
    other than plot_dir. We will also shorten the total fields used.
    In general, this routine will loop through all times in the
    gridded_cube_list for all models in md_list.

    >>> ini_dict = {}
    >>> ini_dict['plot_dir'] = (config.CODE_DIR +
    ...     "/adaqdocs/figures/adaq_plotting/")

    >>> md_list = [md_list[0]] #Shorten to a single model for this example
    >>> for i, cube in enumerate(md_list[0].gridded_cube_list):
    ...     md_list[0].gridded_cube_list[i] = cube[0]

    Do a basic plot using the defaults option of 'AQ'.
    This will use the native coordinate system of the cube for its
    projection and use a linear colour scale. Request output to be in
    subdirectory 'gridded_fields'

    >>> plot_md_gridded_fields(ini_dict, md_list, 'O3',
    ... defaults='AQ', filesuffix='_basicAQ', plot_subdir='gridded_fields')
    ... # doctest: +ELLIPSIS
    Plotting gridded fields
    Saved figure  ...fields/Fieldplot_aqum_oper_O3_201404020000_basicAQ.png

    .. image:: ../adaqdocs/figures/adaq_plotting/gridded_fields/
               Fieldplot_aqum_oper_O3_201404020000_basicAQ.png
       :scale: 75%

    Try adding some extra options in, such as a list of contour levels etc:

    >>> import aq_indices
    >>> import numpy as np
    >>> ini_dict['figsize'] = [6.0, 6.0]
    >>> ini_dict['cbar_orientation'] = 'vertical'
    >>> ini_dict['levels_list'] = np.linspace(0, 220, 11)
    >>> ini_dict['cmap'] = ini_dict.get('cmap', aq_indices.daqi_cmap())
    >>> plot_md_gridded_fields(ini_dict, md_list, 'O3',
    ... defaults='AQ', filesuffix='_fullAQ',
    ... plot_subdir='gridded_fields') # doctest: +ELLIPSIS
    Plotting gridded fields
    Saved figure ...fields/Fieldplot_aqum_oper_O3_201404020000_fullAQ.png

    .. image:: ../adaqdocs/figures/adaq_plotting/gridded_fields/
               Fieldplot_aqum_oper_O3_201404020000_fullAQ.png
       :scale: 75%

    Now also try a basic plot using 'NAME' as defaults option.
    This will instead use a log-based colour scale and use the PlateCarree
    projection. To make room for the longer numbers on the colorbar, we can
    also set the figure size using 'figsize'.

    >>> ini_dict = {}
    >>> ini_dict['plot_dir'] = (config.CODE_DIR +
    ...     "/adaqdocs/figures/adaq_plotting/")
    >>> ini_dict['figsize'] = [9.0, 6.0]
    >>> plot_md_gridded_fields(ini_dict, md_list, 'O3',
    ... defaults='NAME', filesuffix='_basicNAME', plot_subdir='gridded_fields')
    ... # doctest: +ELLIPSIS
    Plotting gridded fields
    Saved figure  ...fields/Fieldplot_aqum_oper_O3_201404020000_basicNAME.png

    .. image:: ../adaqdocs/figures/adaq_plotting/gridded_fields/
               Fieldplot_aqum_oper_O3_201404020000_basicNAME.png
        :scale: 75%
    """

    print('Plotting gridded fields')

    # ---
    # Set up variables which should be used for all models
    # Get some variables from ini_dict if possible.

    # Get extent from ini_dict if possible
    extent = ini_dict.get('extent_list', None)
    if extent:
        # Ensure we have floats
        extent = [float(v) for v in extent]

    # If levels not already set, get levels from ini_dict if possible
    # (and corresponding number of levels)
    levels_dict = ini_dict.get('levels_dict', {})
    if short_name in levels_dict:
        levels = levels_dict[short_name]
    else:
        levels = ini_dict.get('levels_list', None)
    if levels is not None:
        levels = [float(v) for v in levels]
        nlevels = len(levels)
    else:
        nlevels = ini_dict.get('nlevels', 10)

    cmap = ini_dict.get('cmap', 'YlGnBu')
    cbar_orientation = ini_dict.get('cbar_orientation', None)
    cbar_label = ini_dict.get('cbar_label', 'default')
    cbar = ini_dict.get('cbar', True)
    cbar_num_fmt = ini_dict.get('cbar_num_fmt', None)
    rsmc = ini_dict.get('rsmc', False)
    lstep = 1.0 if rsmc else 0.5

    plotdir = ini_dict.get('plot_dir', './')
    if plot_subdir:
        plotdir = os.path.join(plotdir, plot_subdir)
    # Ensure location is a directory by adding '/' on end:
    if plotdir[-1] != os.sep:
        plotdir += os.sep

    mapping = ini_dict.get('mapping', 'countries')
    projection = ini_dict.get('projection', 'PlateCarree')
    plottype = ini_dict.get('plottype', 'pcolormesh')
    mask_name = (True if plottype == 'pcolormesh' else False)

    # If titles are not already set, get these from ini_dict if possible
    title = ini_dict.get('title', None)
    suptitle = ini_dict.get('suptitle', None)

    # If annotation is requested set location to 'right'
    # if annote_location not set
    annote = ini_dict.get('annote', None)
    a_loc = (None if annote is None else 'right')
    annote_location = ini_dict.get('annote_location', a_loc)
    ini_dict['annote_location'] = annote_location

    # Loop over all models and extract cube to plot for required short_name
    if griddedcube is not None:
        md_list = [griddedcube]
    for md in md_list:

        if griddedcube is None:
            if not md.gridded_cube_list:
                continue

            cube = md.extract(short_name=short_name, gridded=True,
                              singlecube=True)
            if cube is None:
                continue
        else:
            cube = griddedcube

        # Now can start setting up plotting for this cube...

        # Deterime which coordinates to use
        # Currently only finds to plot X versus Y
        for coord in cube.dim_coords:
            axis = iris.util.guess_coord_axis(coord)
            if axis == 'X':
                x_coord = coord
            if axis == 'Y':
                y_coord = coord

        # ---
        # Set up field layer
        flr = FieldLayer(cube)

        # Set up layer style, plus any other default changes.
        if defaults == 'NAME':
            # Determine if autozoom is required
            if extent is None:
                autozoom = True
            else:
                autozoom = False
            flr.set_layerstyle(plottype=plottype,
                               colorscale='log',
                               nlevels=int(nlevels),
                               step=lstep,
                               levels=levels,
                               mask=mask_name,
                               autozoom=autozoom,
                               cmap=cmap)
            flr.cbar_orientation = cbar_orientation or 'vertical'
            flr.cbar_label = cbar_label
            flr.cbar = cbar
            flr.cbar_num_fmt = cbar_num_fmt
        elif defaults == 'AQ':
            flr.set_layerstyle(plottype='pcolormesh',
                               colorscale='linear',
                               levels=levels,
                               nlevels=nlevels,
                               mask=True,
                               cmap=cmap)
            flr.cbar_orientation = cbar_orientation
            if flr.levels is None:
                # Setup levels to be same for all fields in cube
                flr.get_linearlevels()
        else:
            flr.set_layerstyle()

        # ---

        # Set the figsize
        figsize = ini_dict.get('figsize', None)

        # ---
        # Loop over fields within layer
        for layer_slice in flr.layer_slice([x_coord.name(), y_coord.name()]):
            fplt = FieldPlot(ini_dict)
            fplt.x_coord = x_coord
            fplt.y_coord = y_coord
            # Don't plot fields which are entirely nan data
            if isinstance(layer_slice.cube.data, np.ma.MaskedArray):
                if np.sum(np.isfinite(layer_slice.cube.data.data)) == 0:
                    # All nan data
                    warnings.warn('All nan data, no gridded field plot created')
                    continue
            else:
                # Not a masked array
                if np.sum(np.isfinite(layer_slice.cube.data)) == 0:
                    # All nan data
                    warnings.warn('All nan data, no gridded field plot created')
                    continue
            fplt.add_layer(layer_slice)

            # Retrieve the extent set from ini_dict if provided or
            # the layer if not. Also calculate central longitude.
            if extent is None:
                layer_extent = flr.cube_extent
            else:
                layer_extent = extent

            # Set up mapping
            if rsmc:
                fplt.rsmc = True
            if defaults == 'NAME':
                # Need to give a central longitude
                if layer_extent is None:
                    clon = 0.0
                else:
                    clon = layer_extent[0] + \
                           (layer_extent[1] - layer_extent[0]) / 2
                fplt.setup_mapping(projection=projection,
                                   extent=layer_extent,
                                   central_longitude=clon,
                                   mapping=mapping,
                                   gridlines=True)
            elif defaults == 'AQ':
                # Use native coordinate system unless set in ini_dict
                fplt.setup_mapping(projection=ini_dict.get('projection', None),
                                   extent=layer_extent,
                                   mapping='coastlines',
                                   gridlines=True)
            else:
                fplt.setup_mapping()

            # Produce and save plot
            fplt.titleprefix = titleprefix
            fplt.suptitle = suptitle
            fplt.title = title
            fplt.plot(figsize=figsize,
                      field_layers=fplt.layer_list)
            if annote is not None:
                plotting_functions.annotate_plot(flr.cube,
                                                 annote, annote_location)

            fplt.save_fig(plotdir=plotdir, filesuffix=filesuffix,
                          tight=tight, verbose=verbose)


def plot_md_multi_fields(ini_dict, md_list, short_name,
                         griddedcube=None, defaults=None,
                         filesuffix='', plot_subdir='',
                         titleprefix=None, tight=False,
                         figsize=None, verbose=1):
    """
    Postage stamp plot for a collection of 2-d (horizontal) fields.
    This assumes various default settings.
    Only does X (eg longitude) versus Y (eg latitude) plots
    for a collection of fields over a 'realization' coordinate.
    For more detailed adjustments and plotting, do not use this routine.

    :param ini_dict: Dictionary of an :class:`inifile` object. Could contain:

                      * extent_list - list of [xmin,xmax,ymin,ymax] in lat/lon
                        coords. Defaults to None if not in ini_dict.
                      * levels_list - list of contouring levels.
                        Defaults to None if not in ini_dict.
                      * nlevels - Number of contour levels if levels_list
                        not set. Defaults to 7 if not in ini_dict.
                      * cmap - Colour map. Defaults to 'YlGnBu' if not
                        in ini_dict.
                      * cbar - Set to True (default) to plot colour bar, set to
                        False to not include a colour bar
                      * cbar_label - label to add to color bar. If 'default'
                        will use the short name.
                      * cbar_num_fmt - Colourbar number format. Default depends
                        on if an AQ or NAME plot
                      * plot_dir - Top level output directory for plots.
                      * title - String to specify title
                      * projection - string containing map projection. Default
                        is 'PlateCarree'. If keyword defaults='AQ' then the
                        default projection is the model native grid (eg rotated
                        pole).

    :param md_list: list of model data objects in the form of an
                    :class:`adaq_data.ADAQData` objects.

    :param short_name: string to match to short_name attribute in cubes.

    :param griddedcube: iris Cube. Overrides md_list and plots gridded cube
                         instead of plotting cubes extracted from
                         md_list.gridded_cube_list.

    :param defaults: Use predefined default setups. Options available are:

                      * None - no default setups.
                      * 'NAME' - default setups for NAME plotting.
                        Adds coastlines, gridlines, central longitude and
                        PlateCarree projection.
                      * 'AQ' - default setups for Air Quality (AQUM etc)
                        plotting. Just adds coastlines and gridlines.
                        Uses native cube projections.

    :param filesuffix: String to add to end of filename for specific naming
                       purposes. By default adds nothing.

    :param plot_subdir: String to add to end of plotdir defined within ini_dict
                        to add a subdirectory to the output plot location.

    :param titleprefix: String to prefix main title with.
                        If prefix ends with ``\\n`` this will put it
                        on the line above the automatic title.

    :param tight: If set to True, will adjust figure size to
                  be tight to bounding box of figure and minimise
                  whitespace. Note this will change size of output file
                  and may not be consistent amongst different plots
                  (particularly if levels on colorbar are different).

    :param figsize: List of two values representing the required size in the
                    x-direction and y-direction of the figure.  If not set,
                    figsize will be automtaically determined using values
                    from the ini_dict.

    :param verbose: Level of print output:

                    * 0 = No extra printing
                    * 1 = Standard level (recommended)
                    * 2 = Extra level for debugging


    >>> import adaq_functions
    >>> import config
    >>> import name_data

    Initialise class:

    >>> name = name_data.NAMEData()

    Read sample data using an ensemble dictionary:

    >>> directory = config.SAMPLE_DATADIR+'name_ensemble/member*/'
    >>> filenames = [directory + 'Fields_grid4*.txt']
    >>> ensemble_dict = {
    ... config.SAMPLE_DATADIR + 'name_ensemble/member0/': ['member0', 0],
    ... config.SAMPLE_DATADIR + 'name_ensemble/member1/': ['member1', 1],
    ... config.SAMPLE_DATADIR + 'name_ensemble/member2/': ['member2', 2],
    ... config.SAMPLE_DATADIR + 'name_ensemble/member3/': ['member3', 3],
    ... config.SAMPLE_DATADIR + 'name_ensemble/member4/': ['member4', 4],
    ... config.SAMPLE_DATADIR + 'name_ensemble/member5/': ['member5', 5],
    ... config.SAMPLE_DATADIR + 'name_ensemble/member6/': ['member6', 6],
    ... config.SAMPLE_DATADIR + 'name_ensemble/member7/': ['member7', 7],
    ... config.SAMPLE_DATADIR + 'name_ensemble/member8/': ['member8', 8],
    ... config.SAMPLE_DATADIR + 'name_ensemble/member9/': ['member9', 9],
    ... config.SAMPLE_DATADIR + 'name_ensemble/member10/': ['member10', 10],
    ... config.SAMPLE_DATADIR + 'name_ensemble/member11/': ['member11', 11],
    ... config.SAMPLE_DATADIR + 'name_ensemble/member12/': ['member12', 12],
    ... config.SAMPLE_DATADIR + 'name_ensemble/member13/': ['member13', 13],
    ... config.SAMPLE_DATADIR + 'name_ensemble/member14/': ['member14', 14],
    ... config.SAMPLE_DATADIR + 'name_ensemble/member15/': ['member15', 15],
    ... config.SAMPLE_DATADIR + 'name_ensemble/member16/': ['member16', 16],
    ... config.SAMPLE_DATADIR + 'name_ensemble/member17/': ['member17', 17]
    ... }

    >>> gcl = name.readdata(filenames, ensemble_dict=ensemble_dict,
    ... field_attributes={'Field Name':'Unnamed Field Req 8'})

    No defaults are set from ini_dict except for plot_dir and nlevels
    (the optimal number of contour levels for postage stamps is around 7).

    >>> ini_dict = {}
    >>> ini_dict['nlevels'] = 7
    >>> ini_dict['plot_dir'] = (config.CODE_DIR +
    ...     "/adaqdocs/figures/adaq_plotting/")

    >>> md_list = [name]

    Plot a sample postage stamp using the 'NAME' defaults option.
    Request output graphic to be written in subdirectory 'postage_stamp'.

    >>> plot_md_multi_fields(ini_dict, md_list, 'CAESIUM-137_AIR_CONCENTRATION',
    ... defaults='NAME', plot_subdir='postage_stamp')
    ... # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    Plotting postage stamps
    Saved figure  .../postage_stamp/Fieldplot_NAME_CAESIUM-137_AIR_CONCENTRATION\
_Boundarylayeraverage_201908290200.png

    .. image:: ../adaqdocs/figures/adaq_plotting/postage_stamp/
               Fieldplot_NAME_CAESIUM-137_AIR_CONCENTRATION\
_Boundarylayeraverage_201908290200.png
       :scale: 75%
    """

    print('Plotting postage stamps')

    # ---
    # Set up variables which should be used for all models
    # Get some variables from ini_dict if possible.

    # Get extent from ini_dict if possible
    extent = ini_dict.get('extent_list', None)
    if extent:
        # Ensure we have floats
        extent = [float(v) for v in extent]

    # If levels not already set, get levels from ini_dict if possible
    levels = ini_dict.get('levels_list', None)
    if levels is not None:
        levels = [float(v) for v in levels]
        nlevels = len(levels)
    else:
        nlevels = ini_dict.get('nlevels', 7)

    cmap = ini_dict.get('cmap', 'YlGnBu')
    cbar = ini_dict.get('cbar', True)
    cbar_label = ini_dict.get('cbar_label', 'default')
    cbar_num_fmt = ini_dict.get('cbar_num_fmt', None)
    rsmc = ini_dict.get('rsmc', False)
    lstep = 1.0 if rsmc else 0.5

    plotdir = ini_dict.get('plot_dir', './')
    if plot_subdir:
        plotdir = os.path.join(plotdir, plot_subdir)
    # Ensure location is a directory by adding '/' on end:
    if plotdir[-1] != os.sep:
        plotdir += os.sep

    mapping = ini_dict.get('mapping', 'countries')
    projection = ini_dict.get('projection', 'PlateCarree')
    plottype = ini_dict.get('plottype', 'pcolormesh')
    mask_name = (True if plottype == 'pcolormesh' else False)

    # If titles are not already set, get these from ini_dict if possible
    title = ini_dict.get('title', None)
    suptitle = ini_dict.get('suptitle', None)

    # Loop over all models and extract cube to plot for required short_name
    if griddedcube is not None:
        md_list = [griddedcube]
    for md in md_list:

        if griddedcube is None:
            if not md.gridded_cube_list:
                continue

            cube = md.extract(short_name=short_name, gridded=True,
                              singlecube=True)
            if cube is None:
                continue
        else:
            cube = griddedcube

        # Now can start setting up plotting for this cube...

        # Deterime which coordinates to use
        # Currently only finds coordinates to plot X versus Y
        # for different ensemble realisations
        for coord in cube.dim_coords:
            axis = iris.util.guess_coord_axis(coord)
            if axis == 'X':
                x_coord = coord
            if axis == 'Y':
                y_coord = coord
            if coord.name() == 'realization':
                r_coord = coord

        # ---
        # Set up field layer
        flr = FieldLayer(cube)

        # Set up layer style, plus any other default changes.
        if defaults == 'NAME':
            # Determine if autozoom is required
            if extent is None:
                autozoom = True
            else:
                autozoom = False
            flr.set_layerstyle(plottype=plottype,
                               colorscale='log',
                               nlevels=int(nlevels),
                               step=lstep,
                               levels=levels,
                               mask=mask_name,
                               autozoom=autozoom,
                               cmap=cmap)
            flr.cbar = cbar
            flr.cbar_label = cbar_label
            flr.cbar_num_fmt = cbar_num_fmt
        elif defaults == 'AQ':
            flr.set_layerstyle(plottype='pcolormesh',
                               colorscale='linear',
                               levels=levels,
                               nlevels=nlevels,
                               mask=True,
                               cmap=cmap)
            if flr.levels is None:
                # Setup levels to be same for all fields in cube
                flr.get_linearlevels()
        else:
            flr.set_layerstyle()

        # ---
        # Loop over fields within layer grouped over the ensemble
        for layer_slice in flr.layer_slice([x_coord.name(), y_coord.name(),
                                            r_coord.name()]):
            fplt = FieldPlot(ini_dict)
            # Don't plot fields which are entirely nan data
            if isinstance(layer_slice.cube.data, np.ma.MaskedArray):
                if np.sum(np.isfinite(layer_slice.cube.data.data)) == 0:
                    # All nan data
                    warnings.warn('All nan data, no gridded field plot created')
                    continue
            else:
                # Not a masked array
                if np.sum(np.isfinite(layer_slice.cube.data)) == 0:
                    # All nan data
                    warnings.warn('All nan data, no gridded field plot created')
                    continue
            fplt.add_layer(layer_slice)

            # Retrieve the extent set from ini_dict if provided or
            # the layer if not. Also calculate central longitude.
            if extent is None:
                layer_extent = flr.cube_extent
            else:
                layer_extent = extent

            # Set up mapping
            if rsmc:
                fplt.rsmc = True
            if defaults == 'NAME':
                # Need to give a central longitude
                if layer_extent is None:
                    clon = 0.0
                else:
                    clon = layer_extent[0] + \
                           (layer_extent[1] - layer_extent[0]) / 2
                fplt.setup_mapping(projection=projection,
                                   extent=layer_extent,
                                   central_longitude=clon,
                                   mapping=mapping,
                                   gridlines=True,
                                   gridlines_xmaxbins=5,
                                   gridlines_ymaxbins=5,
                                   gridlines_xlabel_size=8,
                                   gridlines_ylabel_size=8)
            elif defaults == 'AQ':
                # Use native coordinate system unless set in ini_dict
                fplt.setup_mapping(projection=ini_dict.get('projection', None),
                                   extent=layer_extent,
                                   mapping='coastlines',
                                   gridlines=True)
            else:
                fplt.setup_mapping()

            # Set figsize:
            if figsize is not None:
                figsize = figsize
            elif plot_subdir == 'postage_stamp':
                figsize = [16, 11]
            else:
                figsize = None

            # Produce and save plot
            fig = plt.figure(figsize=figsize)
            fplt.titleprefix = titleprefix
            fplt.suptitle = suptitle
            fplt.title = title
            fplt.multiplot(fig)
            fplt.save_fig(plotdir=plotdir, filesuffix=filesuffix,
                          tight=tight, verbose=verbose)


def plot_qq(ini_dict, od, md_list, short_name):
    """
    Plot quantile-quantile plot for given short_name.

    :param ini_dict: Dictionary of an :class:`inifile` object.
                     Must contain 'plot_dir'.
                     May contain 'line_colours_list'
                     (default: :any:`plotting_functions.COLOURS`).
    :param od: observation data in the form of an
               :class:`adaq_data.ADAQData` object.
    :param md_list: list of model data objects in the form of
                    :class:`adaq_data.ADAQData` objects.
    :param short_name: string to match to short_name attribute in cubes.
    :return: :class:`statistics_plotting.Histogram` object.

    >>> import adaq_functions
    >>> ini_dict, sites_data, od, md_list = adaq_functions.get_exampledata()
    ... # doctest: +ELLIPSIS
    Reading inifile ...example_data_1days.ini
    Number of sites:  5

    For this example only, manually change plot directory to save to gallery:

    >>> import config
    >>> ini_dict['plot_dir'] = (config.CODE_DIR +
    ...     "/adaqdocs/figures/adaq_plotting")

    Now call plot_qq:

    >>> qq = plot_qq(ini_dict, od, md_list, 'O3')  # doctest: +ELLIPSIS
    Plotting quantile-quantile plot
    Saved figure .../figures/adaq_plotting/Quantile-Quantile_O3.png
    >>> np.set_printoptions(formatter={'float':lambda x: '{:5.2f}'.format(x)})
    >>> print(qq.xpercentiles[:5])
    [ 2.00  2.00  2.96 12.64 19.84]
    >>> print(qq.lines[0]['ypercentiles'][:5]) # doctest: +NORMALIZE_WHITESPACE
    [ 0.20  1.98  2.80  4.76  5.79]
    >>> np.set_printoptions()

    Example output file:

    .. image:: ../adaqdocs/figures/adaq_plotting/Quantile-Quantile_O3.png
       :scale: 50%
    """

    print('Plotting quantile-quantile plot')

    colours = ini_dict.get('line_colours_list', plotting_functions.COLOURS)

    obs_cube = od.extract(short_name=short_name, singlecube=True)
    if obs_cube is None:
        warnings.warn("No qq plot for " + short_name + "(No obs)")
        return None
    if np.sum(np.isfinite(obs_cube.data)) == 0:
        warnings.warn("No qq plot for " + short_name + "(obs all nan)")
        return None

    qq = QQPlot(xcube=obs_cube)
    qq.ylabel = 'Model'

    for imd, md in enumerate(md_list):
        mod_cube = md.extract(short_name=short_name, singlecube=True)
        if mod_cube is not None:
            if np.sum(np.isfinite(mod_cube.data)) > 0:
                qq.add_line(mod_cube, colour=colours[imd + 1])

    if not qq.lines:
        warnings.warn("No qq plot for " + short_name + "(No finite model data)")
        return qq

    qq.plot()
    qq.save_fig(plotdir=ini_dict['plot_dir'])

    return qq


def plot_sitescube_maps(ini_dict, od=None, md_list=None, short_name=None,
                        classifications_dict=None, sitescubelist=None,
                        filesuffix='', plot_subdir='', titleprefix=None,
                        tight=False, figsize=None, verbose=1):
    """
    Plot data in sites cubes from both od and md_list on a map for all times.
    Each time and observation and model cube will be in a different plot.
    Can also distinguish between different site types on the plot by using
    the classifications_dict to give a different marker type to the different
    site types which are then also shown in a legend.

    :param ini_dict: Dictionary of an :class:`inifile` object. Could contain:

                      * marker_size - size of marker of plots. Default is 20.
                      * extent_list - list of [xmin,xmax,ymin,ymax] in lat/lon
                        coords. Defaults to None if not in ini_dict, however if
                        obs_fmt is set to 'aurn' or 'camsaqobs' in ini_dict,
                        then a default extent is set to an appropriate domain
                        for the obs type (if all sites within this domain).
                      * levels_list - list of contouring levels.
                        Defaults to None if not in ini_dict.
                      * levels_dict - dictionary whose keys are short_names
                        and whose values are the list of contouring levels
                        for that short_name. If required short_name is available
                        in dictionary then values from this dictionary are
                        used instead of those from levels_list. Defaults
                        to empty dictionary if in not in ini_dict.
                      * nlevels - Number of contour levels if levels_list
                        not set. Defaults to 10 if not in ini_dict.
                      * cmap - Colour map. Defaults to 'YlGnBu' if not
                        in ini_dict.
                      * cbar_orientation - Orientation of color bar
                        ('horizontal' or 'vertical').
                        Defaults to None if not in ini_dict which will add
                        color bar in default position according to levels type
                        (log or linear).
                      * cbar_label - label to add to color bar. If 'default'
                        will use the short name.
                      * cbar_num_fmt - Colourbar Number Format ie. '%.1e'
                      * cbar - Set to True (default) to plot colour bar, set to
                        False to not include a colour bar
                      * colorscale - 'log' or 'linear' colorscale. Defaults to
                        'linear'
                      * annote_location - location of annotation in
                        relation to plot ('right' - next to plot,
                        'below' - below plot).
                      * annote - text to add to plot. Text will added to
                        right hand third of page or bottom third of page
                        depending on annote_location
                      * plot_dir - Top level output directory for plots.
                        Plots will then be saved in plot_dir/gridded_fields.
                        Defaults to './gridded_fields/' if not in ini_dict.
                      * title - String to specify title
                      * mapping - string containing type of mapping to apply
                        Default is 'countries'.
                      * projection - string containing map projection. Default
                        is 'PlateCarree'.

    :param od: observation data in the form of an
               :class:`adaq_data.ADAQData` object. Will plot all times in
               sites_cube_list that match short_name.

    :param md_list: list of model data objects in the form of an
                    :class:`adaq_data.ADAQData` objects. Will plot all times in
                    sites_cube_list that match short_name.

    :param short_name: string to match to short_name attribute in cubes.
                       (Compulsory)

    :param classifications_dict: list of classifications.

    :param sitescubelist: iris Cube. Overrides od and md_list and plots cubes
                          contained in sitescubelist instead.

    :param filesuffix: String to add to end of filename for specific naming
                       purposes. By default adds nothing.

    :param plot_subdir: String to add to end of plotdir defined within ini_dict
                        to give add a subdirectory to the output plot location.

    :param titleprefix: String to prefix main title with.  If prefix ends
                        with ``\\n`` this will put it on the line above the
                        automatic title.

    :param tight: If set to True, will adjust figure size to
                  be tight to bounding box of figure and minimise
                  whitespace. Note this will change size of output file
                  and may not be consistent amongst different plots
                  (particularly if levels on colorbar are different).

    :param figsize: List containing an x-value and y-value to represent the
                    required figsize.

    :param verbose: Level of print output:

                    * 0 = No extra printing
                    * 1 = Standard level (recommended)
                    * 2 = Extra level for debugging

    Example usage... first read in some example data and limit to just the
    first two times rather than plotting full 24 hours of data:

    >>> import adaq_functions
    >>> ini_dict, sites_data, od, md_list = adaq_functions.get_exampledata(
    ... exampletype='short') # doctest: +ELLIPSIS
    Reading inifile .../example_data_1days.ini
    Number of sites:  5

    >>> od.sites_cube_list = iris.cube.CubeList(
    ...     [cube[:,:2] for cube in od.sites_cube_list])
    >>> for md in md_list:
    ...     md.sites_cube_list = iris.cube.CubeList(
    ...         [cube[:,:2] for cube in md.sites_cube_list])

    For this example only, ensure no defaults are set from ini_dict already
    and then set our required options:

    >>> ini_dict = {}
    >>> ini_dict['cbar_orientation'] = 'vertical'
    >>> ini_dict['extent_list'] = [-12, 4.0, 48., 60.]
    >>> ini_dict['marker_size'] = 40
    >>> import config
    >>> ini_dict['plot_dir'] = (config.CODE_DIR +
    ...     "/adaqdocs/figures/adaq_plotting/")

    Set up the classifications dictionary to enable Rural and Urban sites to be
    plotted with different markers and then included on a legend:

    >>> classifications_dict = {
    ... 'Rural': {'site_type': ['REMOTE', 'RURAL'],  'marker': 'o'},
    ... 'Urban': {'site_type': ['SUBURBAN', 'URBAN_BACKGROUND'], 'marker': 's'}}

    Finally pick which short_name to plot:

    >>> short_name = 'O3'

    And now call the plotting routine, requesting an automatic figure size and
    to put the output in a subdirectory 'gridded_fields'

    >>> plot_sitescube_maps(ini_dict, od, md_list, short_name,
    ... classifications_dict, figsize='Auto', plot_subdir='gridded_fields')
    ... # doctest: +ELLIPSIS
    Plotting sitescube maps
    Saved figure  .../Fieldplot_sites_Obs_O3_201404020000.png
    Saved figure  .../Fieldplot_sites_Obs_O3_201404020100.png
    Saved figure  .../Fieldplot_sites_aqum_oper_O3_201404020000.png
    Saved figure  .../Fieldplot_sites_aqum_oper_O3_201404020100.png
    Saved figure  .../Fieldplot_sites_aqum_casestudy_O3_201404020000.png
    Saved figure  .../Fieldplot_sites_aqum_casestudy_O3_201404020100.png

    .. image:: ../adaqdocs/figures/adaq_plotting/gridded_fields/
               Fieldplot_sites_Obs_O3_201404020000.png
       :scale: 75%

    """

    print('Plotting sitescube maps')

    if short_name is None:
        raise ValueError('Require a short_name to enable cube extraction')

    if classifications_dict is None:
        # Set up a default dictionary, containing just one key
        classifications_dict = {'All': {}}

    # Generate list of site cubes (if not passed in) from od and md_list
    if sitescubelist is None:
        sitescubelist = iris.cube.CubeList()

        if od is not None:
            obs_cube = od.extract(short_name=short_name, singlecube=True)
            if obs_cube is not None:
                sitescubelist.append(obs_cube)

        if md_list is not None:
            for md in md_list:
                mod_cube = md.extract(short_name=short_name, singlecube=True)
                if mod_cube is not None:
                    sitescubelist.append(mod_cube)
    if not sitescubelist:
        warnings.warn('No cubes available to plot')
        return

    # Set up variables which should be used for sitecubes
    # Get some variables from ini_dict if possible.

    markersize = ini_dict.get('marker_size', 20)

    # Get extent from ini_dict if possible
    extent = ini_dict.get('extent_list', None)
    if extent:
        # Ensure we have floats
        extent = [float(v) for v in extent]
    # Set up extent (if required) to ensure same extent used across all plots,
    # on the basis of observation type (if given)
    elif 'obs_fmt' in ini_dict:
        # extent_list: [xmin,xmax,ymin,ymax]
        if ini_dict['obs_fmt'] == 'aurn':
            extent = [-11., 3.0, 49., 61.]
        elif ini_dict['obs_fmt'] == 'camsaqobs':
            extent = [-25., 45., 30., 70.]
        else:
            extent = None
        # Check all points are entirely within this domain
        if extent is not None:
            lon_min = np.min([np.min(cube.coord('longitude').points)
                              for cube in sitescubelist])
            lon_max = np.max([np.max(cube.coord('longitude').points)
                              for cube in sitescubelist])
            lat_min = np.min([np.min(cube.coord('latitude').points)
                              for cube in sitescubelist])
            lat_max = np.max([np.max(cube.coord('latitude').points)
                              for cube in sitescubelist])

            if not (lon_min > extent[0] and lon_max < extent[1] and
                    lat_min > extent[2] and lat_max < extent[3]):
                extent = None

    # If levels not already set, get levels from ini_dict if possible
    # (and corresponding number of levels)
    levels_dict = ini_dict.get('levels_dict', {})
    if short_name in levels_dict:
        levels = levels_dict[short_name]
    else:
        levels = ini_dict.get('levels_list', None)
    if levels is not None:
        levels = [float(v) for v in levels]
        nlevels = len(levels)
    else:
        nlevels = ini_dict.get('nlevels', 10)

    cmap = ini_dict.get('cmap', 'YlGnBu')
    cbar_label = ini_dict.get('cbar_label', 'default')
    cbar_num_fmt = ini_dict.get('cbar_num_fmt', None)
    colorscale = ini_dict.get('colorscale', 'linear')

    # Set up linear levels (if required) to ensure same set of levels
    # used across all plots
    if levels is None and colorscale == 'linear':
        data_min = np.nanmin([np.nanmin(cube.data) for cube in sitescubelist])
        data_max = np.nanmax([np.nanmax(cube.data) for cube in sitescubelist])
        levels = np.linspace(data_min, data_max, nlevels)

    plotdir = ini_dict.get('plot_dir', './')
    if plot_subdir:
        plotdir = os.path.join(plotdir, plot_subdir)
    # Ensure location is a directory by adding '/' on end:
    if plotdir[-1] != os.sep:
        plotdir += os.sep

    mapping = ini_dict.get('mapping', 'countries')
    projection = ini_dict.get('projection', 'PlateCarree')

    # If titles are not already set, get these from ini_dict if possible
    title = ini_dict.get('title', None)
    suptitle = ini_dict.get('suptitle', None)

    # If annotation is requested set location to 'right'
    # if annote_location not set
    annote = ini_dict.get('annote', None)
    a_loc = (None if annote is None else 'right')
    annote_location = ini_dict.get('annote_location', a_loc)
    cbar_orientation = ini_dict.get('cbar_orientation', 'vertical')
    ini_dict['annote_location'] = annote_location

    # Loop over all site cubes
    for sitescube in sitescubelist:
        # Loop over all times in cube
        for cube in sitescube.slices_over('time'):

            # Set up plot
            fplt = FieldPlot(ini_dict)

            # Loop over different classifications from the input dictionary;
            # Set up logical to determine if plotting first layer or not
            first_layer = True
            for classification, class_dict in classifications_dict.items():
                # Get cbar properties for first layer, save them to apply
                # to all layers:
                if first_layer:
                    # Get colorbar properties for first layer to be applied to
                    # all layers:
                    all_cbar_orientation = cbar_orientation
                    all_cbar_label = cbar_label
                    all_cbar_num_fmt = cbar_num_fmt
                    cbar = True
                    first_layer = False
                else:
                    # Set boolean status for colorbar plotting to False for
                    # all but first layer so that only one cbar is plotted.
                    cbar = False

                # Get marker from dictionary if available,
                # otherwise default to a circle.
                marker = class_dict.get('marker', 'o')

                # Extract only those sites which match required site types
                # This could be adapted in the future to allow different options
                # to be passed in through the classifications_dict,
                # eg 'site_name'.
                if 'site_type' in class_dict:
                    site_class_cube = cube.extract(iris.Constraint(
                        site_type=class_dict['site_type']))
                else:
                    site_class_cube = cube

                if site_class_cube is not None:
                    if np.sum(np.isfinite(site_class_cube.data)) == 0:
                        # All nan data, continue onto next classification
                        continue

                    # Set up field layer
                    flr = FieldLayer(site_class_cube)
                    flr.set_layerstyle(plottype='scatter_latlon',
                                       colorscale=colorscale,
                                       levels=levels,
                                       nlevels=nlevels,
                                       cmap=cmap)
                    flr.marker = marker
                    flr.markersize = markersize
                    flr.label = classification

                    # Use colorbar properties from first layer for all layers:
                    flr.cbar_orientation = all_cbar_orientation
                    flr.cbar_label = all_cbar_label
                    flr.cbar_num_fmt = all_cbar_num_fmt
                    flr.cbar = cbar

                    # Add layer to plot
                    fplt.add_layer(flr)

            # Check at least one layer has been added:
            if not fplt.layer_list:
                # If not, then continue onto next cube
                continue

            # ---
            # Setup background mapping

            # If extent is not set, then get it from the last layer
            # Required to calculate central longitude
            if extent is None:
                layer_extent = flr.cube_extent
            else:
                layer_extent = extent
            # Now calculate central longitude
            if layer_extent is None:
                clon = 0.0
            else:
                clon = layer_extent[0] + \
                       (layer_extent[1] - layer_extent[0]) / 2

            # Can now setup mapping
            fplt.setup_mapping(projection=projection,
                               extent=layer_extent,
                               central_longitude=clon,
                               mapping=mapping)

            # ---
            # Produce and save plot
            fplt.titleprefix = titleprefix
            fplt.suptitle = suptitle
            fplt.title = title
            fplt.plot(field_layers=fplt.layer_list, figsize=figsize)

            # Add a legend if more than one classification passed in
            if len(classifications_dict) > 1:
                plotting_functions.add_legend_belowaxes(scatterpoints=1)

            if annote is not None:
                plotting_functions.annotate_plot(flr.cube,
                                                 annote, annote_location)

            fplt.save_fig(plotdir=plotdir, filesuffix=filesuffix,
                          tight=tight, verbose=verbose,
                          fileprefix='Fieldplot_sites_')


def plot_soccer_plots(ini_dict, od, md_list, short_name, classifications_dict,
                      filesuffix=''):
    """
    Plot soccer plots for given short_name.

    :param ini_dict: Dictionary of an :class:`inifile` object.
                     May contain 'line_colours_list'
                     (default: :any:`plotting_functions.COLOURS`).
    :param od: observation data in the form of an
               :class:`adaq_data.ADAQData` object.
    :param md_list: list of model data objects in the form of an
                    :class:`adaq_data.ADAQData` objects.
    :param short_name: string to match to short_name attribute in cubes.
    :param classifications_dict: Dictionary, whose keys are classification
                                 names, and whose values are dictionaries,
                                 which should include:

                                  * 'site_type' - list of site_types to use to
                                    constrain cubes by to get this particular
                                    classification. (See example below).
                                  * 'marker' - (optional) to indicate the marker
                                    on the scatter plot - defaults to circle
                                    if not set.

    :param filesuffix: String to add to end of filename for specific naming
                    purposes. By default adds nothing.
    :return: :class:`statistics_plotting.SoccerPlot` object

    .. note:: Currently, this routine limits on the basis of site_type.
              This routine could easily be adapted to constraint on other
              coordinates.

    >>> import adaq_functions
    >>> ini_dict, sites_data, od, md_list = adaq_functions.get_exampledata()
    ... # doctest: +ELLIPSIS
    Reading inifile ...example_data_1days.ini
    Number of sites:  5

    For this example only, manually change plot directory to save to gallery:

    >>> import config
    >>> ini_dict['plot_dir'] = (config.CODE_DIR +
    ...     "/adaqdocs/figures/adaq_plotting")

    Set up a dictionary of classifications:

    >>> classifications_dict = {
    ... 'Rural': {'site_type': ['REMOTE','RURAL'], 'marker': 'o'},
    ... 'Urban': {'site_type': ['SUBURBAN','URBAN_BACKGROUND'], 'marker': 's'}
    ... }

    Then produce the soccer plots:

    >>> sp = plot_soccer_plots(ini_dict, od, md_list, 'O3',
    ... classifications_dict) # doctest: +ELLIPSIS
    Plotting Soccer Plot
    Saved figure  .../figures/adaq_plotting/Soccer_Plot_O3.png
    >>> print(sp.title)
    Soccer Plot of O3
     01/04/2014 23:00 to 03/04/2014 00:00
    >>> np.set_printoptions(formatter={'float':lambda x: '{:5.2f}'.format(x)})
    >>> print(sp.lines[0]['cube'].coord('mnmb').points)
    [ 0.24 -0.41 -0.20 -0.41]
    >>> np.set_printoptions()

    Example output file:

    .. image:: ../adaqdocs/figures/adaq_plotting/Soccer_Plot_O3.png
       :scale: 50%
    """

    print('Plotting Soccer Plot')

    colours = ini_dict.get('line_colours_list', plotting_functions.COLOURS)

    sp = SoccerPlot()
    sp.gridlines = False

    # Extract observation cube for this short_name
    obs = od.extract(short_name=short_name, singlecube=True)

    for imd, md in enumerate(md_list):

        # Extract model data cube for this short_name
        mod = md.extract(short_name=short_name, singlecube=True)

        if obs is not None and mod is not None:
            # Calculate statistics and get a cube with statistics included in
            # For optimisation, calculate stats over all sites and then limit
            # by classification later - intersection of cube times
            # is the slowest component
            stats = sp.get_stats(obs, mod)
        else:
            continue

        for classification, class_dict in classifications_dict.items():
            # Get marker from dictionary if available,
            # otherwise default to a circle.
            marker = class_dict.get('marker', 'o')

            # Extract only those sites which match required site types
            stats_classification = stats.extract(iris.Constraint(
                site_type=class_dict['site_type']))

            if stats_classification is not None:
                # Then add as a line (scatter point) to soccer plot
                sp.add_line(stats_classification,
                            marker=marker,
                            colour=colours[imd + 1],
                            label=mod.attributes['label']
                            + '-' + classification)

    # Plot
    # Alter number of columns in legend
    if len(md_list) == 1:
        legendcols = len(classifications_dict)
    else:
        legendcols = len(md_list)

    if not sp.lines:
        warnings.warn("No soccer plot for " + short_name + "(No obs/models)")
        return sp

    sp.plot(legendcols=legendcols)
    sp.save_fig(plotdir=ini_dict['plot_dir'],
                filename='Soccer_Plot_' + short_name + filesuffix + '.png')

    return sp


def plot_timeseries(ini_dict, od, md_list, short_name, site, verbose=1):
    """
    Produce times-series plot for given short_name and site.

    :param ini_dict: Dictionary of an :class:`inifile` object.
                     May contain 'line_colours_list'
                     (default: :any:`plotting_functions.COLOURS`).
    :param od: observation data in the form of an
               :class:`adaq_data.ADAQData` object.
    :param md_list: list of model data objects in the form of an
                    :class:`adaq_data.ADAQData` objects.
    :param short_name: string to match to short_name attribute in cubes
    :param site: element from numpy ndarray containing site information data
                 from a :class:`sites_info.SitesInfo` object.
    :param verbose: Level of print output:

                    * 0 = No extra printing
                    * 1 = Standard level (recommended)
                    * 2 = Extra level for debugging

    :return: :class:`timeseries_plot.TimeSeriesPlot` object.

    >>> import adaq_functions
    >>> ini_dict, sites_data, od, md_list = adaq_functions.get_exampledata()
    ... # doctest: +ELLIPSIS
    Reading inifile ...example_data_1days.ini
    Number of sites:  5

    For this example only, manually change plot directory to save to gallery:

    >>> import config
    >>> ini_dict['plot_dir'] = (config.CODE_DIR +
    ...     "/adaqdocs/figures/adaq_plotting")

    Now produce timeseries:

    >>> tsp = plot_timeseries(ini_dict, od, md_list, 'O3',
    ... sites_data[0]) # doctest: +ELLIPSIS
    Plotting timeseries
    Saved figure  .../adaq_plotting/timeseries/Aberdeen_O3.png

    .. image:: ../adaqdocs/figures/adaq_plotting/timeseries/Aberdeen_O3.png
       :scale: 50%

    >>> print(tsp.title)
    Aberdeen (Urban background)
    mass concentration of ozone in air
    >>> print(tsp.lines[0]['cube']) # doctest: +ELLIPSIS
    mass_concentration_of_ozone_in_air / (ug/m3) (time: 25)
         Dimension coordinates:
              time                                    x
         Scalar coordinates:
              abbrev: ABD
              latitude: 57.15750122 degrees
              longitude: -2.09388876 degrees
              site_altitude: 20 m
              site_id: 35790611.14...
              site_name: Aberdeen
              site_type: URBAN_BACKGROUND
         Attributes:
              Conventions: CF-1.5
              label: Obs
              short_name: O3
              source: AURN
         Cell methods:
              mean: time (1 hour)
    """

    print('Plotting timeseries')

    colours = ini_dict.get('line_colours_list', plotting_functions.COLOURS)

    # Set up time-series-plot
    tsp = TimeSeriesPlot()

    # Create a site cube from the obs ...
    obs_cube = od.extract(short_name=short_name,
                          singlecube=True,
                          abbrev=site['abbrev'])
    # Check not None:
    if obs_cube is not None:
        # Check not all nans - if ok, then can add cube
        if np.isfinite(obs_cube.data).sum() != 0:
            tsp.add_line(obs_cube, colour=colours[0],
                         linestyle='--')

    # Create a site cube from the models ...
    for imd, md in enumerate(md_list):
        mod_cube = md.extract(short_name=short_name,
                              singlecube=True,
                              abbrev=site['abbrev'])
        if mod_cube is not None:
            # Check not all nans - if ok, then can add cube
            if np.isfinite(mod_cube.data).sum() != 0:
                tsp.add_line(mod_cube,
                             colour=colours[imd + 1])

    if not tsp.lines:
        warnings.warn("No timeseries plot for " + short_name + \
                      "(No obs/models)")
        return tsp

    # Create a plot and save
    tsp.plot()
    tsp.save_fig(plotdir=ini_dict['plot_dir'] + '/timeseries/', verbose=verbose)

    return tsp


def plot_timeseries_aggregate_stats(ini_dict, adaqdata_list, short_name,
                                    stat='NANMEAN', gridded=False,
                                    linestyles=None, colours=None,
                                    filesuffix='', verbose=1):
    """
    Produce timeseries plots of aggregated-statistics for a single short-name
    for all ADAQData objects in input list.
    Note this routine aggregates data by collapsing over all coordinates
    other than time. Data in the input list could be observations or model
    data.

    :param ini_dict: Dictionary of an :class:`inifile` object.
                     Must contain 'plot_dir'.
                     May contain 'line_colours_list'
                     (default: :any:`plotting_functions.COLOURS`),
                     although it will be overridden by the keyword argument
                     `colours` if present.
    :param adaqdata_list: list or model and/or observation data in the form of
                          a :class:`adaq_data.ADAQData` objects.
    :param short_name: string to match to short_name attribute in cubes
    :param stat: string, statistic to use for aggregating data. Available
                 statistics are given in `cube_statistics.CUBE_AGGREGATORS`
                 dictionary.
    :param gridded: logical. If False, then the sites_cube_list is used from
                    the :class:`adaq_data.ADAQData` objects. If True, then
                    the gridded_cube_list is instead used.
    :param linestyles: list of matplotlib line styles to use for each adaqdata
                       object - to match ordering in adaqdata_list.
                       If None, then uses solid lines.
    :param colours: list of colours to use for each adaqdata object - to match
                    ordering in adaqdata_list. Uses the ini_dict if not set.
    :param filesuffix: String to add to end of filename for specific naming
                       purposes. By default adds nothing.
    :param verbose: Level of print output:

                    * 0 = No extra printing
                    * 1 = Standard level (recommended)
                    * 2 = Extra level for debugging

    Example of producing timeseries of mean data over all sites:

    >>> import adaq_functions
    >>> ini_dict, sites_data, od, md_list = adaq_functions.get_exampledata(
    ... sites_cube_list=True, gridded_cube_list=True) # doctest: +ELLIPSIS
    Reading inifile ...example_data_1days.ini
    Number of sites:  5

    For this example only, manually change plot directory to save to gallery:

    >>> import config
    >>> ini_dict['plot_dir'] = (config.CODE_DIR +
    ... "/adaqdocs/figures/adaq_plotting")

    >>> tsp = plot_timeseries_aggregate_stats(ini_dict, [od]+md_list, 'O3',
    ... stat='NANMEAN', gridded=False,
    ... linestyles=['--']+['-']*len(md_list))  # doctest: +ELLIPSIS
    Plotting timeseries of aggregated cube statistics
    Saved figure  .../Timeseries_of_sites_nanmean_O3.png

    Note we have used nanmean here as some of the sites in the data do not
    have O3 observations (so have nan in their data arrays), hence using
    'MEAN' would give nan data in the output plots.
    We also modified linestyles to plot observations with dashed lines and
    the model data with solid lines.

    Example output plot:

    .. image:: ../adaqdocs/figures/adaq_plotting/Timeseries_of_sites_nanmean_O3.png
       :scale: 50%

    Example of instead plotting maximum value across the gridded fields
    (note we also here change the colours used to avoid using the initial
    black colour from :any:`plotting_functions.COLOURS`):

    >>> tsp = plot_timeseries_aggregate_stats(ini_dict, md_list, 'O3',
    ... stat='MAX', gridded=True, colours=plotting_functions.COLOURS[1:])
    ... # doctest: +ELLIPSIS
    Plotting timeseries of aggregated cube statistics
    Saved figure  .../Timeseries_of_field_max_O3.png

    Example output plot:

    .. image:: ../adaqdocs/figures/adaq_plotting/Timeseries_of_field_max_O3.png
       :scale: 50%


    """

    print('Plotting timeseries of aggregated cube statistics')

    if colours is None:
        colours = ini_dict.get('line_colours_list', plotting_functions.COLOURS)

    # Set up time-series-plot
    tsp = TimeSeriesPlot()

    for i, adaqdata in enumerate(adaqdata_list):
        cube = adaqdata.extract(short_name=short_name,
                                singlecube=True,
                                gridded=gridded)
        if cube is not None:
            # Collapse over coordinates other than time, using required statistic
            # First determine which these coordinates are:
            coords = [coord for coord in cube.coords(dim_coords=True)
                      if coord.name() != 'time']
            collapsed_cube = cube.collapsed(
                coords, cube_statistics.CUBE_AGGREGATORS[stat])

            # Check not all nans - if ok, then can add cube
            if np.isfinite(collapsed_cube.data).sum() != 0:
                if linestyles is not None:
                    linestyle = linestyles[i]
                else:
                    linestyle = '-'
                tsp.add_line(collapsed_cube, colour=colours[i],
                             linestyle=linestyle)

    if not tsp.lines:
        warnings.warn("No Timeseries of aggregrate stats for " \
                      + stat + ' + ' + short_name)
        return tsp

    cubetype = 'sites'
    if gridded:
        cubetype = 'field'
    tsp.title = (stat.lower().replace('nan', '') + ' across ' + cubetype +
                 '\n' + cube.name().replace('_', ' '))

    tsp.plot()

    # Save figure
    filename = ('Timeseries_of_' + cubetype + '_' + stat.lower() +
                '_' + short_name + filesuffix + '.' + tsp.figformat)
    tsp.save_fig(plotdir=ini_dict['plot_dir'], filename=filename,
                 verbose=verbose)

    return tsp


def plot_timeseries_of_stats(ini_dict, od, md_list, short_name, threshold=None,
                             filesuffix=''):
    """
    Plot timeseries of statistics for a given short_name.

    :param ini_dict: Dictionary of an :class:`inifile` object.
                     Must contain:

                     * 'plot_dir'
                     * 'timeseries_of_stats_list', a list of the required
                       statistics to plot on a timeseries.
                       Valid entries for this list are given in
                       :data:`timeseries_stats.STATS_INFO` (see output from
                       example below).

                     May also contain 'line_colours_list'
                     (default: :any:`plotting_functions.COLOURS`).
    :param od: observation data in the form of an :class:`adaq_data.ADAQData`
               object
    :param md_list: list of model data objects in the form of an
                    :class:`adaq_data.ADAQData` objects.
    :param short_name: string to match to short_name attribute in cubes.
    :param threshold: Value to use as a threshold when calculating
                      threshold-based statistics.
    :param filesuffix: String to add to end of filename for specific naming
                       purposes. By default adds nothing.

    :return: :class:`timeseries_plot.TimeSeriesPlot` object.

    >>> import adaq_functions
    >>> ini_dict, sites_data, od, md_list = adaq_functions.get_exampledata()
    ... # doctest: +ELLIPSIS
    Reading inifile ...example_data_1days.ini
    Number of sites:  5

    For this example only, manually change plot directory to save to gallery:

    >>> import config
    >>> ini_dict['plot_dir'] = (config.CODE_DIR +
    ...     "/adaqdocs/figures/adaq_plotting")

    Add the required statistics to be output to the ini_dict:

    >>> ini_dict['timeseries_of_stats_list'] = ['bias', 'rmse']

    Note the valid entries for this list are given by the keys in
    :data:`timeseries_stats.STATS_INFO`:

    >>> print(sorted(timeseries_stats.STATS_INFO.keys()))
    ... # doctest: +NORMALIZE_WHITESPACE
    ['bias', 'correlation', 'fac2', 'falsealarmrate', 'falsealarmratio',
    'fge', 'hitrate', 'ioa', 'maxmod', 'maxobs', 'mdi', 'meanmod', 'meanobs',
    'mge', 'mnmb', 'nmb', 'nmge', 'npts', 'nsites', 'o<t_m<t', 'o<t_m>=t',
    'o>=t_m<t', 'o>=t_m>=t', 'odds_ratio', 'orss', 'pc95mod', 'pc95obs',
    'perc_correct', 'perc_over', 'perc_under', 'percmod', 'percobs', 'rmse',
    'sdmod', 'sdobs', 'threshold', 'units']

    Then produce the plots:

    >>> plot_timeseries_of_stats(ini_dict, od, md_list, 'O3')
    ... # doctest: +ELLIPSIS
    Plotting timeseries of statistics
    Saved figure .../Timeseries_of_bias_O3.png
    Saved figure .../Timeseries_of_rmse_O3.png

    Example output files:

    .. image:: ../adaqdocs/figures/adaq_plotting/Timeseries_of_bias_O3.png
       :scale: 50%
    .. image:: ../adaqdocs/figures/adaq_plotting/Timeseries_of_rmse_O3.png
       :scale: 50%

    """

    print('Plotting timeseries of statistics')

    colours = ini_dict.get('line_colours_list', plotting_functions.COLOURS)

    obs_cube = od.extract(short_name=short_name, singlecube=True)
    if obs_cube is None:
        # Can't calculate statistics
        warnings.warn("No statistics timeseries plot for "
                      + short_name + "(No obs)")
        return
    if np.sum(np.isfinite(obs_cube.data)) == 0:
        # Can't calculate statistics
        warnings.warn("No statistics timeseries plot for "
                      + short_name + "(No finite obs data)")
        return

    # Set up the coordinates to iterate through combinations of, which
    # should be the equivalent of stepping though time.
    coords = [coord.name() for coord in obs_cube.coords(dim_coords=True)
              if coord.name() != 'time']

    # First calculate statistics for all models
    statscube_list = []  # List of cubes, one per model.
    colours_used = []
    for imd, md in enumerate(md_list):

        mod_cube = md.extract(short_name=short_name, singlecube=True)
        if mod_cube is not None:
            if np.sum(np.isfinite(obs_cube.data)) == 0:
                # All nan data
                continue
            # Set up empty cubelist for this model
            statscube_list_mod = iris.cube.CubeList()

            # Step through combinations of the coordinates setup above.
            # This is the equivalent of stepping through time.
            # (Like for slices_over('time') but need 2 cubes at the same time.)
            for oc, mc in iris.iterate.izip(obs_cube, mod_cube, coords=coords):

                # Check times are the same (best done by converting to datetime
                # first to avoid checking equality between floats:
                oc_tpt = oc.coord('time').units.num2date(
                    mc.coord('time').points)
                mc_tpt = mc.coord('time').units.num2date(
                    mc.coord('time').points)

                if oc_tpt != mc_tpt:
                    print('obs points:', oc.coord('time').units.num2date(
                        oc.coord('time').points))
                    print('mod points:', mc.coord('time').units.num2date(
                        mc.coord('time').points))
                    raise ValueError("Obs and Model cube don't have same times "
                                     "in - run add_missing_times first")

                # Calculate all statistics for this time
                stats = timeseries_stats.TimeSeriesStats(oc, mc, mdi=np.nan,
                                                         threshold=threshold)
                percentile_value = ini_dict.get('percentile', 95)
                stats.get_basic_stats(percentile_value=percentile_value)
                statscube = stats.convert_to_cube()
                # Add cube to cubelist for this model
                statscube_list_mod.append(statscube)

            # Convert model cubelist to a single cube with a time dim and append
            # to the statscube_list
            statscube_list.append(statscube_list_mod.concatenate_cube())
            colours_used.append(colours[imd + 1])

    # Now loop over the required statistics and plot them
    for statistic in ini_dict.get('timeseries_of_stats_list', []):
        tsp_statistic(statscube_list, statistic,
                      plotdir=ini_dict['plot_dir'],
                      filesuffix=filesuffix,
                      colours_list=colours_used)


def plot_tseries_multiple_snames(ini_dict, od, md_list, site, short_name_list,
                                 short_name_group=None, verbose=1):
    """
    Produce timeseries at a given site, displaying multiple short_names
    on the same plot.

    :param ini_dict: Dictionary of an :class:`inifile` object.
                     May contain 'line_colours_list'
                     (default: :any:`plotting_functions.COLOURS`).
    :param od: observation data in the form of an
               :class:`adaq_data.ADAQData` object.
    :param md_list: list of model data objects in the form of an
                    :class:`adaq_data.ADAQData` objects.
    :param site: element from numpy ndarray containing site information data
                 from a :class:`sites_info.SitesInfo` object.
    :param short_name_list: List of strings to match to short_name attribute in
                            cubes - each item will be plotted as a separate
                            line.
    :param short_name_group: String describing the group of short_names passed
                             in - this will be used for titles and for
                             filenames.
    :param verbose: Level of print output:

                    * 0 = No extra printing
                    * 1 = Standard level (recommended)
                    * 2 = Extra level for debugging

    :return: :class:`timeseries_plot.TimeSeriesPlot` object.

    Example of plotting PM10 and PM2p5 on the same plot:

    >>> import adaq_functions
    >>> ini_dict, sites_data, od, md_list = adaq_functions.get_exampledata(
    ...    exampletype="full") # doctest: +ELLIPSIS
    Reading inifile .../example_data_5days.ini
    Number of sites:  5

    For this example only, manually change plot directory to save to gallery:

    >>> import config
    >>> ini_dict['plot_dir'] = (config.CODE_DIR + "/adaqdocs/" +
    ...     "figures/adaq_plotting")

    Now produce plot:

    >>> tsp = plot_tseries_multiple_snames(ini_dict, od, md_list,
    ... sites_data[0], short_name_list = ['PM10','PM2p5'],
    ... short_name_group='PM Components') # doctest: +ELLIPSIS
    Plotting timeseries for multiple short_names
    Saved figure  .../timeseries/Aberdeen_PM_Components_aqum_oper.png
    Saved figure  .../timeseries/Aberdeen_PM_Components_aqum_casestudy.png
    Saved figure  .../timeseries/Aberdeen_PM_Components_Obs.png

    Note that this routine is generally expected to be called by looping over
    keys in ini_dict['timeseries_multiple_short_names_dict'] which correspond
    to the short_name_group and the values in this dictionary then correspond
    to the short_name_list for this group.

    .. image:: ../adaqdocs/figures/adaq_plotting/timeseries/\
Aberdeen_PM_Components_aqum_casestudy.png
       :scale: 50%

    """

    print('Plotting timeseries for multiple short_names')

    if short_name_group is None:
        short_name_group = "_".join(short_name_list)

    colours = ini_dict.get('line_colours_list', plotting_functions.COLOURS)

    # Observations and model data are treated the same
    # - all put into different plots
    for adaq in md_list + [od]:

        # If no cubes in sites_cube_list then can't plot!
        if adaq.sites_cube_list is None:
            continue
        if not adaq.sites_cube_list:
            continue

        # Set up timeseries plot
        tsp = TimeSeriesPlot()
        tsp.phenomena_short_name = short_name_group
        # Set title
        title = ''
        if 'site_name' in site.dtype.names:
            title += site['site_name'].replace('_', ' ') + ' '
        if 'site_type' in site.dtype.names:
            title += '(' + site['site_type'].replace('_', ' ').capitalize()
            title += ')' + '\n'
        label = adaq.sites_cube_list[0].attributes['label']
        title += label + ': ' + short_name_group.replace('_', ' ')
        tsp.title = title

        # Add each short_name on as an individual line
        for isn, short_name in enumerate(short_name_list):

            cube = adaq.extract(short_name=short_name, abbrev=site['abbrev'],
                                singlecube=True)

            if cube is not None:
                if not np.isnan(cube.data).all():  # Don't plot if all nan
                    tsp.add_line(cube,
                                 colour=colours[isn + 1],
                                 label=short_name)

        if not tsp.lines:
            warnings.warn("No Timeseries of multiple short names for " \
                          + label + ' + ' + short_name_group)
            continue

        # Produce plot
        tsp.plot()

        # Save plot
        filename = site['site_name'] + '_' + \
            short_name_group.replace(' ', '_') + '_' + \
            label + '.' + tsp.figformat
        tsp.save_fig(plotdir=ini_dict['plot_dir'] + '/timeseries/',
                     filename=filename, verbose=verbose)


def plot_trajectories(ini_dict, md, verbose=1):
    """
    Code to plot trajectories

    Example trajectory plot

    >>> import config
    >>> import trajectory_data

    Create an ini_dict and add a models_dir_list and a plot_dir

    >>> ini_dict = {}
    >>> ini_dict['plot_dir'] = (config.CODE_DIR + "/adaqdocs/" +
    ...     "figures/adaq_plotting")
    >>> ini_dict['models_dir_list'] = [config.SAMPLE_DATADIR +
    ... "/name_trajectory/Data_Traj_C1_*.txt"]
    >>> ini_dict['marker_interval'] = 12

    Read in data:

    >>> md = trajectory_data.TrajData()
    >>> md.readdata(ini_dict['models_dir_list'])
    [<iris 'Cube' of U Turb / (m/s) (time: 193)>,
    <iris 'Cube' of U Turb / (m/s) (time: 193)>]

    Now produce plot:

    >>> tjhp = plot_trajectories(ini_dict, md) # doctest: +ELLIPSIS
    Saved figure  .../TrajectoryPlot.png

    .. image:: ../adaqdocs/figures/adaq_plotting/TrajectoryPlot.png
       :scale: 50%
    """

    # Use first cube to determine the spacing of the markers on the plot
    # --------------------------------------------------------------------
    time = md.trajectory_cube_list[0].coord('time')
    mk_int = ini_dict.get('marker_interval', 6)
    mk_start, mk_interval = plotting_functions.marker_spacing(time, mk_int)

    # Check for dateline/meridian crossing and compute extent
    # --------------------------------------------------------------------
    dateline = 0
    meridian = 0
    if 'extent_list' in ini_dict:
        # assumes user will use 0:360 across dateline and -180:180 across
        # meridian
        extent = [float(v) for v in ini_dict['extent_list']]
        if extent[0] < 180 and extent[1] > 180:
            dateline = 1
        if extent[0] < 0 and extent[1] > 0:
            meridian = 1

    else:
        for cube in md.trajectory_cube_list:
            lon_points = cube.coord('longitude').points
            if np.any(lon_points > 175) and np.any(lon_points < -175):
                dateline = 1
            if np.any(abs(lon_points) < 5):
                meridian = 1

    clon = 0
    if dateline == 1 and meridian == 0:
        for cube in md.trajectory_cube_list:
            lon_points = cube.coord('longitude').points
            lon_points = [lp + 360 if lp < 0 else lp for lp in lon_points]
            cube.coord('longitude').points = lon_points

            if verbose != 0:
                print('Centering plot on dateline')
            clon = 180

    # Compute extent or use extent from ini_dict
    new_extent = compute_trajectory_extent(md.trajectory_cube_list,
                                           dateline, meridian)
    extent = ini_dict.get('extent_list', new_extent)
    extent = [float(v) for v in extent]

    # Use a TrajectoryPlot for the x-y plot
    tjp = TrajectoryPlot()
    tjp.clon = clon
    tjp.extent = extent
    tjp.mapping = ini_dict.get('mapping', None)
    tjp.mobrand = ini_dict.get('mo_branding', False)
    tjp.release_info = ini_dict.get('release_info_list', None)
    tjp.title = compute_traj_direction(md.trajectory_cube_list)
    tjp.rsmc = ini_dict.get('rsmc', None)
    tjp.annote = ini_dict.get('annote', None)

    for cube in md.trajectory_cube_list:
        longitude = cube.coord('longitude')
        latitude = cube.coord('latitude')

        tjp.add_line(cube, longitude, latitude,
                     wrap=True,
                     marker='*',
                     markeredgecolor='none',
                     markevery=(mk_start, abs(int(mk_interval))),
                     markersize=8)

    tjp.plot()

    # Use a TimeSeriesPlot for the time-z plot
    thp = TimeSeriesPlot()
    rsmc = ini_dict.get('rsmc', False)
    if rsmc:
        thp.fig = plt.subplot2grid((3, 3), (2, 0), colspan=2)
    elif ini_dict.get('annote', False):
        thp.fig = plt.subplot2grid((3, 3), (2, 0), colspan=2)
    else:
        thp.fig = plt.subplot2grid((3, 1), (2, 0))
    thp.legend = False
    thp.title = ''

    for cube in md.trajectory_cube_list:

        z_coord_names = ['Z (Pa)', 'flight_level', 'height', 'altitude']

        # First determine name of z coordinate
        for coord in cube.coords():
            if iris.util.guess_coord_axis(coord) == 'Z':
                z_coord = coord
            elif coord.name() in z_coord_names:
                z_coord = coord

        # If Z coordinate is Pressure - convert to hPa
        if z_coord.name() == 'Z (Pa)':
            cube.coord('Z (Pa)').convert_units('hPa')

        time = cube.coord('time')

        thp.add_line(cube, time, z_coord,
                     wrap=True,
                     marker='*',
                     markeredgecolor='none',
                     markevery=(mk_start, abs(int(mk_interval))),
                     markersize=8)

    thp.ylabel = Z_COORD_YLABEL[z_coord.name()]
    thp.plot()

    # Add annotation if requested
    annote = ini_dict.get('annote', None)
    annote_location = ini_dict.get('annote_location', 'right')
    run_time = ini_dict.get('run_time', None)
    if annote is not None:
        plotting_functions.annotate_plot(cube, annote, annote_location,
                                         run_time=run_time)

    # Add a supertitle if a title is requested
    if 'title' in ini_dict:
        plt.suptitle(ini_dict['title'], fontsize=14)

    # Save plot
    if 'plotname' in ini_dict:
        filename = ini_dict['plotname']
    else:
        filename = 'TrajectoryPlot.' + tjp.figformat
    tjp.save_fig(plotdir=ini_dict['plot_dir'],
                 filename=filename, verbose=verbose)


if __name__ == "__main__":
    import doctest

    doctest.testmod()
