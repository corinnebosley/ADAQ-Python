"""
Plotting of ADAQData objects which have a vertical/height component to them.
These routines all generally take the following arguments:

 * **ini_dict -** Dictionary of a :class:`inifile` object.
   Must contain 'plot_dir'.
 * **md_list -** list of model data objects in the form of an
   :class:`adaq_data.ADAQData` objects.
 * **short_name -** string to match to short_name attribute in cubes

Plus some other routine-specific arguments.
"""

from __future__ import print_function

import ast
import os
import warnings
import cartopy.crs as ccrs
import iris.plot as iplt
import iris.util
import matplotlib.pyplot as plt
import numpy as np

import cube_functions
import field_layer
import field_plot
import inifile
import trajectory_plot
from plotting_functions import COLOURS, units_str


def plot_section_sample_locations(section, plotdir='./', extent=None,
                                  mapping='coastlines', title='',
                                  verbose=1):
    """
    Plot location of sample locations from a section cube. Mainly for use
    from within :func:`plot_md_cross_section()`.

    :param section: Cross section cube with 'i_sample_point', 'latitude' and
                    'longitude' coordinates.
    :param plotdir: Top level output directory for plots.
    :param extent: List of [xmin,xmax,ymin,ymax] in lat/lon
                   coords. Defaults to None if not in ini_dict. Only used
                   for plotting location of waypoints.
    :param mapping: String, for background mapping, defaults to
                    'coastlines', but 'countries','states', 'wms' are also
                    available. Only used for plotting location of waypoints.
    :param title: String to prefix to 'Sample point locations' and use as
                  plot title
    :param verbose: Level of print output:

                    * 0 = No extra printing
                    * 1 = Standard level (recommended)
                    * 2 = Extra level for debugging

    >>> import config
    >>> sample_datadir = config.SAMPLE_DATADIR+'gridded_cube_list/'
    >>> cube = iris.load_cube(sample_datadir+'aqum_oper_1days.nc',
    ... 'mass_concentration_of_ozone_in_air')
    >>> waypoints = [
    ... {'latitude': 50.8, 'longitude': -1.8},
    ... {'latitude': 51.2, 'longitude':	-1.2},
    ... {'latitude': 51.4, 'longitude':	-0.9}]
    >>> section = cube_functions.extract_section(cube, waypoints)
    >>> section.attributes['label'] = 'aqum_section'

    >>> plot_section_sample_locations(section,
    ... plotdir=config.CODE_DIR + "/adaqdocs/figures/adaq_plotting")
    ... # doctest: +ELLIPSIS
    Saved figure  .../CrossSectionSamplePoints_aqum_section.png

    .. image:: ../adaqdocs/figures/adaq_plotting/
               CrossSectionSamplePoints_aqum_section.png

    """

    tp = trajectory_plot.TrajectoryPlot()
    #Extract just the first slice containing all the sample points
    #(eg remove time/Z dimension)
    section_slice = next(section.slices('i_sample_point'), None)

    if section_slice is None:
        return

    if extent is None:
        # Check for dateline/meridian crossing and compute extent
        dateline = 0
        meridian = 0

        lon_points = section_slice.coord('longitude').points
        if np.any(lon_points > 175) and np.any(lon_points < -175):
            dateline = 1
        if np.any(abs(lon_points) < 5):
            meridian = 1

        clon = 0
        if dateline == 1 and meridian == 0:
            lon_points = section_slice.coord('longitude').points
            lon_points = [lp + 360 if lp < 0 else lp for lp in lon_points]
            section_slice.coord('longitude').points = lon_points

            if verbose != 0:
                print('Centering plot on dateline')
            clon = 180

        extent = trajectory_plot.compute_trajectory_extent(
            [section], dateline, meridian, minheight=0.)

    tp.extent = extent
    tp.clon = clon
    tp.mapping = mapping
    if title:
        tp.title = title + '\nSample point locations'
    else:
        tp.title = 'Sample point locations'
    #Set fixed figure size - note compute_trajectory_extent ensures it should
    #always fit this size.
    tp.fig = plt.figure(figsize=[8, 6])
    plt.axes(projection=ccrs.PlateCarree(central_longitude=clon))
    tp_ax = tp.fig.gca()
    tp.add_line(section_slice,
                section.coord('longitude'),
                section.coord('latitude'))
    # Plot the line:
    tp.plot()

    #Plot labelled sample points
    xcoordname, ycoordname = cube_functions.guess_coord_names(section_slice,
                                                              ['X', 'Y'])

    proj = section_slice.coord(xcoordname).coord_system.as_cartopy_projection()
    for sample_slice in section_slice.slices_over('i_sample_point'):
        sample_lat = sample_slice.coord(ycoordname).points[0]
        sample_lon = sample_slice.coord(xcoordname).points[0]
        i_sample_point = sample_slice.coord('i_sample_point').points[0]

        tp_ax.plot(sample_lon, sample_lat,
                   c='black', marker='.', transform=proj)
        tp_ax.text(sample_lon, sample_lat, ' ' + str(i_sample_point),
                   transform=proj,
                   horizontalalignment='left',
                   verticalalignment='center',
                   fontsize='x-small')

    tp.save_fig(plotdir=plotdir,
                filename=('CrossSectionSamplePoints_' +
                          section.attributes['label'] + '.png'),
                verbose=verbose)


def plot_cube_cross_section(cube, waypoints, short_name, ini_dict,
                            bl_depth=None, titleprefix=None,
                            plotdir=None, filesuffix='', tight=False,
                            verbose=1):
    """
    Plot vertical cross section of a cube along a given set of waypoints.

    Defaults (where parameters are set to None or 'default') are set to the
    defaults in :class:`field_layer.FieldLayer` or

    :class:`field_plot.FieldPlot`.

    :param cube: Iris cube which should contain X,Y and Z coordinates, plus
                 'short_name' and 'label' attributes.

    :param waypoints: List of dictionaries, whose keys are 'latitude' and
                      'longitude' and whose values are the lat/lon points that
                      should be included in the cross-section.

    :param short_name: string to match to short_name attribute in cubes.

    :param ini_dict: Dictionary of values used to define appearance of plot.

    :param bl_depth: Iris cube containing boundary layer depth. This should
                     contain X and Y coordinates, plus 'short_name' and 'label'
                     attributes. If given, the boundary layer height will also
                     be plotted on the cross-section plot.

    :param titleprefix: String to prefix to default title (time of plot).

    :param plotdir: Output directory for plot to be saved to.

    :param filesuffix: String to add to end of filename for specific naming
                       purposes. By default adds nothing.

    :param tight: If set to True, will adjust figure size to
                  be tight to bounding box of figure and minimise
                  whitespace. Note this will change size of output file
                  and may not be consistent amongst different plots
                  (particularly if levels on colorbar are different).

    :param verbose: Level of print output:

                    * 0 = No extra printing
                    * 1 = Standard level (recommended)
                    * 2 = Extra level for debugging

    :returns: Iris cube containing the section that was plotted. This will
              have additional coordinate 'i_sample_points' instead of X and Y
              coordinates.

    Load example cubes:

    >>> import config
    >>> import adaq_functions
    >>> sample_datadir = config.SAMPLE_DATADIR + 'aqum_output/oper/3d/'
    >>> tcube = iris.load_cube(sample_datadir +
    ... 'prodm_op_aqum_20170701_18.006.pp', 'air_temperature')
    >>> ini_dict = inifile.get_inidict(defaultfilename=
    ... 'adaqcode/adaq_vertical_plotting.ini') # doctest: +ELLIPSIS
    Reading inifile .../adaqcode/adaq_vertical_plotting.ini

    >>> tcube.attributes['short_name'] = 'T'
    >>> tcube.attributes['label'] = 'aqum'
    >>> tcube.coord('time').points = tcube.coord('time').bounds[:, 1]
    >>> print(tcube.summary(True)) # doctest: +NORMALIZE_WHITESPACE
    air_temperature / (K) (time: 2; model_level_number: 63; grid_latitude:
    182; grid_longitude: 146)
    >>> blcube = iris.load_cube(sample_datadir +
    ... 'prods_op_aqum_20170701_18.000.pp',
    ... 'atmosphere_boundary_layer_thickness')
    >>> blcube.attributes['short_name'] = 'bl_depth'
    >>> blcube.attributes['label'] = 'aqum'

    Set up dictionary of way points and then plot cross section

    >>> waypoints = [{'latitude':51.1, 'longitude':-0.6},
    ... {'latitude':51.7, 'longitude':0.2},
    ... {'latitude':52.0, 'longitude':-1.0}]
    >>> plotdir = config.CODE_DIR + "/adaqdocs/" + "figures/adaq_plotting"
    >>> section = plot_cube_cross_section(tcube, waypoints, 'T', ini_dict,
    ... bl_depth=blcube, plotdir=plotdir) # doctest: +ELLIPSIS
    Saved figure  .../CrossSection_aqum_T_201707020300.png
    Saved figure  .../CrossSection_aqum_T_201707020600.png

    .. image:: ../adaqdocs/figures/adaq_plotting/
               CrossSection_aqum_T_201707020600.png
       :scale: 75%

    >>> print(section.summary(True)) # doctest: +NORMALIZE_WHITESPACE
    air_temperature / (K) (time: 2; level_height: 26; i_sample_point: 30)

    """

    #Get some variables from ini_dict if possible.
    max_height = ini_dict.get('max_height', None)
    if max_height is not None:
        max_height = float(max_height)

    #If levels not already set, get levels from ini_dict if possible
    #(and corresponding number of levels)
    levels_dict = ini_dict.get('levels_dict', {})
    if short_name in levels_dict:
        levels = levels_dict[short_name]
    else:
        levels = ini_dict.get('levels_list', None)
    if levels is not None:
        levels = [float(v) for v in levels]
    if levels is not None:
        nlevels = len(levels)
    else:
        nlevels = ini_dict.get('nlevels', 10)

    cmap = ini_dict.get('cmap', 'YlGnBu')
    cbar_label = ini_dict.get('cbar_label', 'default')
    cbar = ini_dict.get('cbar', True)
    cbar_num_fmt = ini_dict.get('cbar_num_fmt', None)
    line_colours_list = ini_dict.get('line_colours_list', COLOURS)
    line_colour = line_colours_list[0]
    cbar_orientation = ini_dict.get('cbar_orientation', 'vertical')

    # Extract a section along the waypoints
    section = cube_functions.extract_section(cube, waypoints)

    #Plot against a sensible vertical coordinate
    if section.coords('model_level_number') and section.coords('level_height'):
        section.remove_coord('model_level_number')
        iris.util.promote_aux_coord_to_dim_coord(section, 'level_height')

    zcoordname, = cube_functions.guess_coord_names(section, ['Z'])
    if zcoordname is None:
        print('Not plotting cross section for '+cube.attributes['short_name']
              +' from '+cube.attributes['label']+' (no vertical coordinate)')
        return None
    scoordname = 'i_sample_point'

    #Limit section to maximum required height
    if max_height is not None:
        if section.coords('level_height'):
            section = section.extract(iris.Constraint(
                level_height=lambda c: c <= max_height))
        else:
            raise UserWarning('Cannot limit to max_height as level_height'
                              + 'coordinate does not exist')
    if bl_depth:
        bl_depth_section = cube_functions.extract_section(
            bl_depth, waypoints)


    #---
    #Set up field layer
    flr = field_layer.FieldLayer(section)
    flr.set_layerstyle(plottype='pcolormesh',
                       colorscale='linear',
                       levels=levels,
                       nlevels=nlevels,
                       mask=True,
                       cmap=cmap)

    flr.cbar = cbar
    flr.cbar_orientation = cbar_orientation
    flr.cbar_label = cbar_label
    flr.cbar_num_fmt = cbar_num_fmt

    # Set colour for boundary layer line plots
    line_colour = line_colour if line_colour else 'black'

    #---
    #Loop over fields within layer (sample point & Z coords)
    for layer_slice in flr.layer_slice([scoordname, zcoordname]):
        fplt = field_plot.FieldPlot(ini_dict)
        #Don't plot fields which are entirely nan data
        if isinstance(layer_slice.cube.data, np.ma.MaskedArray):
            if np.sum(np.isfinite(layer_slice.cube.data.data)) == 0:
                #All nan data
                warnings.warn('All nan data, no gridded field plot created')
                continue
        else:
            #Not a masked array
            if np.sum(np.isfinite(layer_slice.cube.data)) == 0:
                #All nan data
                warnings.warn('All nan data, no gridded field plot created')
                continue
        fplt.add_layer(layer_slice)

        if titleprefix is None:
            titleprefix = ('Cross section for '
                           + cube.name().replace('_', ' ') + '\n')
        fplt.titleprefix = titleprefix
        # Set a specific figsize for this type of plot:
        fplt.plot(figsize=[15.0, 6.0])

        # Get the coordinate that was actually used on the vertical axis.
        vert_axis_coord = iplt._get_plot_defn(flr.cube,
                                              iris.coords.BOUND_MODE,
                                              ndims=2).coords[0]

        plt.gca().set_ylabel(vert_axis_coord.name().replace('_', ' ')
                             + ' (' + units_str(str(vert_axis_coord.units)) + ')')
        if vert_axis_coord.has_bounds():
            plt.gca().set_ylim([0, vert_axis_coord.bounds[-1, 1]])
        plt.gca().set_xlabel('Section point')
        plt.gca().set_xlim([-0.5, 29.5])  # n_sample_points=30.

        # Also plot BL depth at the same time
        if bl_depth:
            time = layer_slice.cube.coord('time').units.num2date(
                layer_slice.cube.coord('time').points[0])
            bl_depth_slice = bl_depth_section.extract(iris.Constraint
                                                      (time=time))
            if bl_depth_slice:  # found matching time
                # Plot
                iplt.plot(bl_depth_slice, label='Boundary Layer Height',
                          color=line_colour, linestyle='--')
                # Add label
                plt.gca().legend()

        fplt.save_fig(plotdir=plotdir, fileprefix='CrossSection_',
                      filesuffix=filesuffix,
                      tight=tight, verbose=verbose)

    return section


def plot_md_cross_section(ini_dict, md_list, short_name, plot_bl_depth=True,
                          griddedcube=None, filesuffix='', plot_subdir='',
                          titleprefix=None, tight=False, verbose=1):
    """
    Plot vertical cross section along a given set of waypoints. Also plot
    location of waypoints.

    :param ini_dict: Dictionary of an :class:`inifile` object. Could contain:

                      * waypoints_list (compulsory) - List of dictionaries,
                        whose keys are 'latitude' and 'longitude' and whose
                        values are the lat/lon points that should be included
                        in the cross-section.
                      * max_height - Maximum height to display on plot. Only
                        used if cube has a 'level_height' coordinate.
                      * levels_list - list of contouring levels.
                        Defaults to None if not in ini_dict.
                      * levels_dict - dictionary whose keys are short_names
                        and whose values are the list of contouring levels
                        for that short_name. If required short_name is available
                        in dictionary then values from this dictionary are
                        used instead of those from levels_list. Defaults
                        to empty dictionary if in not in ini_dict.
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
                      * plot_dir - Top level output directory for plots.
                      * title - String to specify title
                      * extent_list - list of [xmin,xmax,ymin,ymax] in lat/lon
                        coords. Defaults to None if not in ini_dict. Only used
                        for plotting location of waypoints.
                      * mapping - String, for background mapping, defaults to
                        'coastlines', but 'countries','states', 'wms' are also
                        available. Only used for plotting location of waypoints.

    :param md_list: list of model data objects in the form of
                    :class:`adaq_data.ADAQData` objects.
    :param short_name: string to match to short_name attribute in cubes.
    :param plot_bl_depth: logical. If True (default), plot boundary depth over
                          the top if available as a gridded cube for
                          corresponding model data in md_list.
    :param griddedcube: iris Cube. Overrides md_list and plots gridded cube
                        instead of plotting cubes extracted from
                        md_list.gridded_cube_list. This still allows for
                        values to be used from ini_dict.
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

    Example of plotting cross section.

    First get some example data:

    >>> import adaq_functions
    >>> import config
    >>> ini_dict, sites_data, od, md_list = adaq_functions.get_exampledata(
    ... exampletype='3d', gridded_cube_list=True, sites_cube_list=False)
    ... # doctest: +ELLIPSIS
    Reading inifile .../adaqcode/example_data_3d.ini

    Set up extra variables in ini_dict, including the required waypoints_list:

    >>> ini_dict['waypoints_list'] = [
    ... {'latitude':51.1, 'longitude':-0.6},
    ... {'latitude':51.7, 'longitude':0.2},
    ... {'latitude':52.0, 'longitude':-1.0}]
    >>> ini_dict['max_height'] = 5000.
    >>> ini_dict['plot_dir'] = (config.CODE_DIR + "/adaqdocs/" +
    ... "figures/adaq_plotting")
    >>> ini_dict['mapping'] = 'wms'

    >>> plot_md_cross_section(ini_dict, md_list, 'O3') # doctest: +ELLIPSIS
    Plotting cross sections for O3
    Saved figure  .../CrossSection_aqum_oper_O3_201409012100.png
    Saved figure  .../CrossSection_aqum_oper_O3_201409020000.png
    Saved figure  .../CrossSectionSamplePoints_aqum_oper.png

    .. image:: ../adaqdocs/figures/adaq_plotting/
               CrossSection_aqum_oper_O3_201409012100.png
       :scale: 75%

    .. image:: ../adaqdocs/figures/adaq_plotting/
               CrossSectionSamplePoints_aqum_oper.png
       :scale: 75%



    """

    print('Plotting cross sections for '+short_name)

    #---
    #Set up variables which should be used for all models
    plotdir = ini_dict.get('plot_dir', './')
    if plot_subdir:
        plotdir = os.path.join(plotdir, plot_subdir)
    #Ensure location is a directory by adding '/' on end:
    if plotdir[-1] != os.sep:
        plotdir += os.sep
    title = ini_dict.get('title', None)

    #---
    #Variables which are only required for plotting location of waypoints/
    #sample points

    #Get extent from ini_dict if possible
    extent = ini_dict.get('extent_list', None)
    if extent:
        #Ensure we have floats
        extent = [float(v) for v in extent]

    mapping = ini_dict.get('mapping', 'countries')

    waypoints = ini_dict['waypoints_list']
    for i, waypoint in enumerate(waypoints):
        if not isinstance(waypoint, dict):
            waypoints[i] = ast.literal_eval(waypoint)

    #Loop over all models and extract cube to plot for required short_name
    if griddedcube is not None:
        md_list = [griddedcube]

    for md in md_list:

        if griddedcube is None:
            if not md.gridded_cube_list:
                continue

            cube = md.extract(short_name=short_name, gridded=True,
                              singlecube=True)
            if cube is None:
                print('Not plotting cross section for '+short_name+' from '
                      + md.gridded_cube_list[0].attributes['label']
                      + ' (no cube for this short_name)')
                continue
        else:
            cube = griddedcube

        if plot_bl_depth:
            bl_depth = md.extract(short_name='BL_depth', gridded=True,
                                  singlecube=True)
        else:
            bl_depth = None

        #------------------
        #Plot cross section
        section = plot_cube_cross_section(cube, waypoints, short_name,
                                          ini_dict,
                                          bl_depth=bl_depth,
                                          titleprefix=titleprefix,
                                          plotdir=plotdir,
                                          filesuffix=filesuffix,
                                          tight=tight,
                                          verbose=verbose)

        #-----------------------------------
        #Also plot location of sample points using a trajectory plot

        if section:
            plot_section_sample_locations(section,
                                          plotdir=plotdir,
                                          extent=extent,
                                          mapping=mapping,
                                          title=title,
                                          verbose=verbose)


if __name__ == '__main__':

    import doctest
    doctest.testmod()
