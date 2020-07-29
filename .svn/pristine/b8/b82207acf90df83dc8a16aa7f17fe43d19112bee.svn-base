#!/usr/bin/env python
"""
Script to generate RSMC-style output from NAME fields files. Note that the
RSMC requirements are very specific so most of the ini_dict parameters are
set in the script rather than the inifile.

The script expects to find one set of three files containing integrated air
concentration and deposition, and one set of files (with a different name)
containing air concentration output at a higher temporal resolution
and a set of trajectory files.

Running the code
~~~~~~~~~~~~~~~~
To run the RSMC plotting code you need to create an inifile containing
information about the location of your name output. The plots are then
created by typing

.. code-block:: ksh

   rsmc_plot.py [inifilename]

where inifilename includes the full path to the ini file.

Alternatively, from within other python code, run using:

.. code-block:: python

    from adaqscripts import name_field_plot
    rsmc_plot.rsmc_plotting(inifilename='inifilename')

or if passing an ini_dict rather than and ini file

.. code-block:: python

    from adaqscripts import rsmc_plot
    # Set up ini_dict, e.g.
    ini_dict = {}
    ini_dict['key'] = value
    rsmc_plot.rsmc_plotting(ini_dict=ini_dict)

Contents of the inifile
~~~~~~~~~~~~~~~~~~~~~~~
The inifile needs to contain the following items:

+-----------------+------------------------------------------------------------+
|models_dir_list  |List of the full directory path and start of filename of    |
|                 |the files containing the integrated air concentration and   |
|                 |deposition. Linux wildcards such as * and ? may be used     |
+-----------------+------------------------------------------------------------+
|models_fmt_list  |List of the format of each of the sets of files in the      |
|                 |models_dir_list. e.g. `name` or `traj`                      |
+-----------------+------------------------------------------------------------+
|arrival_time     |Logical (True or False). Default is False - don't calculate |
|                 |arrival time                                                |
+-----------------+------------------------------------------------------------+
|plot_dir         |Full directory path of the location you would like your     |
|                 |plots to be saved to.                                       |
+-----------------+------------------------------------------------------------+

The inifile may also contain any of the following items:

+--------------------+---------------------------------------------------------+
|models_list         |List of names for the model. This can be any text and the|
|                    |text will be added to the filename. So, for example, if  |
|                    |you have done one NAME run with convection on and one    |
|                    |NAME run with convection off you could set the           |
|                    |*models_list* to "no_convection" and "convection"        |
|                    |respectively.                                            |
+--------------------+---------------------------------------------------------+
|source_note         |A note about the origin of the source term to be added   |
|                    |to the annotation next to the plot                       |
|                    |Default is: Results based on initial default values      |
+--------------------+---------------------------------------------------------+
|suptitle            |A string to be used for the suptitle (or super title)    |
|                    |for all plots                                            |
+--------------------+---------------------------------------------------------+

"""
from __future__ import division
from __future__ import print_function

import datetime
import os
import sys

from six.moves.builtins import zip
from six.moves.builtins import range

import iris
import numpy as np
#Find adaq directory -
# should be in directory above the one containing this file
adaq_path = os.path.dirname(os.path.realpath(__file__))+'/../'
sys.path.append(adaq_path)
import adaqcode
from adaqcode import cube_time, name_data, adaq_cmap,\
    field_plot, field_layer
from adaqcode.adaq_plotting import plot_md_gridded_fields, plot_trajectories
import adaqcode.plotting_functions as adaq_plot


def plot_arrival_time(time_cube, ini_dict, plotname):
    """
    Plots the arrival time
    """

    print("Plotting time of arrival")

    # Setup field_layer
    flr = adaqcode.field_layer.FieldLayer(time_cube)
    flr.set_layerstyle(plottype='pcolormesh',
                       colorscale='linear',
                       levels=ini_dict['levels'],
                       mask=True,
                       cmap=adaq_cmap.rsmctoa_cmap())
    flr.cbar_label = ''

    # Plotting
    extent = ini_dict['extent_list']
    clon = extent[0] + \
        (extent[1] - extent[0])/2

    for layer_slice in flr.layer_slice(['longitude', 'latitude']):
        fplt = field_plot.FieldPlot(ini_dict)
        fplt.add_layer(layer_slice)
        fplt.rsmc = True
        fplt.setup_mapping(projection='PlateCarree',
                           extent=extent,
                           central_longitude=clon,
                           mapping='states',
                           gridlines=True)
        fplt.suptitle = ini_dict.get('suptitle', None)

        if layer_slice.cube.coord('Threshold').units == 'Bq s/m3':
            fplt.title = 'Time of Arrival (integrated)'
        else:
            fplt.title = 'Time of Arrival (averaged)'

        fplt.plot(figsize=[11, 8])

        adaq_plot.annotate_plot(flr.cube, 'RSMC', 'right')

        fplt.save_fig(plotdir=ini_dict['plot_dir'], filename=plotname,
                      tight=False, verbose=1)



def rsmc_plotting(inifilename=None, ini_dict=None, verbose=1):
    '''
    Top-level plotting routine to generate field plots

    :param inifilename: String giving filename of ini file.

    :param ini_dict: "inidict" dictionary of preferred values to be read in.

    :param verbose: Level of print output:

                    * 0 = No extra printing
                    * 1 = Standard level (recommended)
                    * 2 = Extra level for debugging

    :returns: (ini_dict, sites_data, od, md_list)

     * **ini_dict** - Dictionary of a :class:`inifile` object
     * **md_list** - Model Data list - each item will be name
       in :class:`adaq_data` format

    **Method:**
     * Reads in inifile
     * Reads in model data
     * Plots data using style given in inifile

    >>> ini_dict, md_list = rsmc_plotting()
    ... # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    Current version of python: ...
    [...]
    Current version of iris: ...
    Beginning at  ...
    Reading inifile .../adaqscripts/rsmc_plot.ini
    Getting model data for  name  at  ...
    Getting model data for  name  at  ...
    Getting model data for  name  at  ...
    <class '....NAMEData'> - Subclass of ADAQData - Contains:
    sites_cube_list:
    < No cubes >
    gridded_cube_list:
    0: CAESIUM-137_DOSAGE / (Bq s/m3)\
    (time: 3; latitude: 480; longitude: 640)
    1: CAESIUM-137_TOTAL_DEPOSITION / (Bq/m2)\
    (time: 3; latitude: 480; longitude: 640)
    trajectory_cube_list:
    < No cubes >
    <class '....NAMEData'> - Subclass of ADAQData - Contains:
    sites_cube_list:
    < No cubes >
    gridded_cube_list:
    0: CAESIUM-137_DOSAGE / (Bq s/m3)\
    (time: 12; latitude: 480; longitude: 640)
    trajectory_cube_list:
    < No cubes >
    <class '....TrajData'> - Subclass of ADAQData - Contains:
    sites_cube_list:
    < No cubes >
    gridded_cube_list:
    < No cubes >
    trajectory_cube_list:
    0: U Turb / (m/s)                      (time: 289)
    1: U Turb / (m/s)                      (time: 289)
    2: U Turb / (m/s)                      (time: 289)
    Plotting gridded fields
    Saved figure  \
    .../Fieldplot_name_CAESIUM-137_DOSAGE_From0_0to500_0magl\
_201702241200_COL.png
    Plotting gridded fields
    Saved figure  \
    .../Fieldplot_name_CAESIUM-137_DOSAGE_From0_0to500_0magl\
_201702241200_BW.png
    Plotting gridded fields
    Saved figure  \
    .../Fieldplot_name_CAESIUM-137_DOSAGE_From0_0to500_0magl\
_201702251200_COL.png
    Plotting gridded fields
    Saved figure  \
    .../Fieldplot_name_CAESIUM-137_DOSAGE_From0_0to500_0magl\
_201702251200_BW.png
    Plotting gridded fields
    Saved figure  \
    .../Fieldplot_name_CAESIUM-137_DOSAGE_From0_0to500_0magl\
_201702261200_COL.png
    Plotting gridded fields
    Saved figure  \
    .../Fieldplot_name_CAESIUM-137_DOSAGE_From0_0to500_0magl\
_201702261200_BW.png
    Plotting gridded fields
    Saved figure  \
    .../Fieldplot_name_CAESIUM-137_TOTAL_DEPOSITION_Boundarylayer\
_201702241200_COL.png
    Plotting gridded fields
    Saved figure  \
    .../Fieldplot_name_CAESIUM-137_TOTAL_DEPOSITION_Boundarylayer\
_201702241200_BW.png
    Plotting gridded fields
    Saved figure  \
    .../Fieldplot_name_CAESIUM-137_TOTAL_DEPOSITION_Boundarylayer\
_201702251200_COL.png
    Plotting gridded fields
    Saved figure  \
    .../Fieldplot_name_CAESIUM-137_TOTAL_DEPOSITION_Boundarylayer\
_201702251200_BW.png
    Plotting gridded fields
    Saved figure  \
    .../Fieldplot_name_CAESIUM-137_TOTAL_DEPOSITION_Boundarylayer\
_201702261200_COL.png
    Plotting gridded fields
    Saved figure  \
    .../Fieldplot_name_CAESIUM-137_TOTAL_DEPOSITION_Boundarylayer\
_201702261200_BW.png
    Plotting time of arrival
    Saved figure  .../Fieldplot_name_Time_of_arrival_of_CAESIUM-137_00.png
    Plotting time of arrival
    Saved figure  .../Fieldplot_name_Time_of_arrival_of_CAESIUM-137_01.png
    Plotting time of arrival
    Saved figure  .../Fieldplot_name_Time_of_arrival_of_CAESIUM-137_02.png
    Saved figure .../TrajectoryPlot.png
    Finished at  ...

    '''

    print('Current version of python:', sys.version)
    print('Current version of iris:', iris.__version__)
    print('Beginning at ', datetime.datetime.now())


    # Read ini file
    if ini_dict is None:
        ini_dict = adaqcode.inifile.get_inidict(
            inifilename=inifilename,
            defaultfilename='adaqscripts/rsmc_plot.ini')
    else:
        if 'field_attribute_dict' not in ini_dict:
            ini_dict['field_attribute_dict'] = {}

    # Add some fixed variables to ini-file
    ini_dict['mapping'] = 'states'
    ini_dict['cbar_label'] = ' '
    ini_dict['title'] = 'name_verbose'
    ini_dict['annote_location'] = 'right'
    ini_dict['annote'] = 'RSMC'
    ini_dict['nlevels'] = 5
    ini_dict['plottype'] = 'contourf_edge'
    ini_dict['cbar_orientation'] = 'horizontal'
    ini_dict['rsmc'] = True
    ini_dict['figsize'] = [11, 8]

    # Read in model data
    md_list = adaqcode.adaq_functions.get_models(ini_dict)
    for md in md_list:
        print(md)

    # Check for presence of short_name_list and populate if not there
    # First cube list only
    if 'short_name_list' not in ini_dict:
        short_name_list = []
        short_name_list.extend([cube.attributes['short_name']
                                for cube in md_list[0].gridded_cube_list])
        ini_dict['short_name_list'] = short_name_list

    # Determine levels and extents of each cube in first cube list
    md = md_list[0]
    extent_list_list = []
    levels_list_list = []
    for short_name in ini_dict['short_name_list']:

        cube = md.extract(short_name=short_name, gridded=True,
                          singlecube=True)

        for itime in range(3):

            cube1 = cube[itime, :, :]
            if np.max(cube1.data) < 1e-23:
                continue

            flr = field_layer.FieldLayer(cube1)
            flr.set_layerstyle(plottype='contourf',
                               colorscale='log',
                               nlevels=5,
                               step=1.0,
                               autozoom=True)
            extent_list_list.append(flr.cube_extent)
            levels_list_list.append(flr.levels[0])

    ini_dict['extent_list'] = adaq_plot.combine_extents(extent_list_list)

    # Plot air concentration and deposition
    for short_name in ini_dict['short_name_list']:
        cube = md_list[0].extract(short_name=short_name, gridded=True,
                                  singlecube=True)
        cube.attributes['source_note'] = ini_dict.get(
            'source_note', "Results based on initial default values")

        for itime in range(3):
            md_s = name_data.NAMEData()
            md_s.gridded_cube_list = iris.cube.CubeList([cube[itime, :, :]])
            if "DEPOSITION" in short_name:
                ini_dict['cmap'] = adaq_cmap.rsmcdep_cmap()
            else:
                ini_dict['cmap'] = adaq_cmap.rsmcac_cmap()

            plot_md_gridded_fields(ini_dict, [md_s], short_name,
                                   filesuffix='_COL',
                                   defaults='NAME', verbose=verbose)

            ini_dict['cmap'] = adaq_cmap.rsmcgrey_cmap()
            plot_md_gridded_fields(ini_dict, [md_s], short_name,
                                   filesuffix='_BW',
                                   defaults='NAME', verbose=verbose)

    # Now read in higher temporal resolution fields and compute time of arrival
    arrival_time = ini_dict.get('arrival_time', False)
    if arrival_time:

        # Plot time of arrival
        plotnums = ['00', '01', '02']
        threshold = np.min(levels_list_list[0:3])
        gcl = md_list[1].gridded_cube_list[0]

        # Modify threshold if cube contains average air concentrations
        # Note that to keep the threshold as a round factor of 10 we divide
        # by 1e6 instead of 86400 (or 24*3600)
        if "average" in gcl.attributes['Time Av or Int']:
            threshold = threshold/1.0e6
        levelses = [[-0.01, 0, 6, 12, 18, 24],
                    [0, 24, 30, 36, 42, 48],
                    [0, 48, 54, 60, 66, 72]]
        # To plot times of arrival in three parts need to only provide the
        # timesteps to be plotted. Otherwise the last contour continues
        # up to the maximum time of arrival on all plots
        maxts = [5, 9, 13]

        for plotnum, levels, maxt  in zip(plotnums, levelses, maxts):
            time_cube = cube_time.plume_arrival_time(gcl[0:maxt, :, :],
                                                     threshold)

            # Add source note to time_cube
            time_cube.attributes['source_note'] = ini_dict.get(
                'source_note', "Results based on initial default values")

            ini_dict['levels'] = levels

            plotname = 'Fieldplot_{}_{}_{}.png'.format(
                ini_dict.get('models_list', ['name'])[1],
                time_cube.name().replace(' ', '_'),
                plotnum)

            plot_arrival_time(time_cube, ini_dict, plotname)

    # Now plot trajectories
    ini_dict['title'] = ini_dict['suptitle']
    ini_dict['annote'] = 'RSMC_trajectories'
    ini_dict['marker_interval'] = 6
    ini_dict['run_time'] = md_s.gridded_cube_list[0].attributes['Run time']
    del ini_dict['extent_list']
    if arrival_time:
        plot_trajectories(ini_dict, md_list[2], verbose=verbose)

    else:
        plot_trajectories(ini_dict, md_list[1], verbose=verbose)

    print('Finished at ', datetime.datetime.now())

    return ini_dict, md_list


if __name__ == '__main__':

    rsmc_plotting()

    #import doctest
    #doctest.testmod()
    #doctest.run_docstring_examples(get_ini, globals())
