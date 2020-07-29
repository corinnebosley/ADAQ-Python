#!/usr/bin/env python
"""
Script to generate postage stamp plots from ensemble NAME gridded files.

Running the code
~~~~~~~~~~~~~~~~
To run the postage stamp plotting code you need to create an inifile containing
information about the location of your name output and the columns you would
like to plot (see next section). The plots are then created by typing

.. code-block:: ksh

   ./name_postage_stamp_plot.py [inifilename]

where inifilename includes the full path to the ini file.

Alternatively, from within other python code, run using:

.. code-block:: python

    from adaqscripts import name_postage_stamp_plot
    name_postage_stamp_plot.postage_stamp_plotting(inifilename='inifilename')

or, if passing an ini_dict rather than an ini file,

.. code-block:: python

    from adaqscripts import name_postage_stamp_plot
    # Set up ini_dict, e.g.
    ini_dict = {}
    ini_dict['key'] = value
    name_postage_stamp_plot.postage_stamp_plotting(ini_dict=ini_dict)

Contents of the inifile
~~~~~~~~~~~~~~~~~~~~~~~

Essential Items
===============

The inifile needs to contain the following items:

+---------------+--------------------------------------------------------------+
|models_dir_list|List of the full directory paths of the data to be read in.   |
|               |NAME output files should be grouped using a separate subfolder|
|               |of a common top-level folder for each ensemble member. These  |
|               |subfolders should be indexed by consecutive integers starting |
|               |from zero with the index labelled as a * in the filepath.     |
|               |Linux wildcards such as * and ? may be used in the path names.|
+---------------+--------------------------------------------------------------+
|models_fmt_list|List of the data formats to be read in. It should be "name"   |
|               |for any output from NAME.                                     |
+---------------+--------------------------------------------------------------+
|models_list    |List of names for the model. This can be any text and the     |
|               |text will be added to the filename. So, for example, if       |
|               |you have done one NAME run with convection on and one         |
|               |NAME run with convection off you could set the                |
|               |*models_list* to "no_convection" and "convection"             |
|               |respectively.                                                 |
+---------------+--------------------------------------------------------------+


Optional Items
==============

The inifile may also contain any of the following items:


+--------------------+---------------------------------------------------------+
|field_attribute_dict|Dictionary of the field attributes of the NAME data to   |
|                    |be plotted. The field attributes are basically extracted |
|                    |from the column headers of the NAME file. See            |
|                    |:ref:`field_attributes_ref` for a more comprehensive list|
+--------------------+---------------------------------------------------------+
|z_level_list        |List of the heights which you would like to plot.        |
|                    |For NAME III format files this is the height given in    |
|                    |the column header  without units. For NAME II format     |
|                    |files you can provide any number between the two numbers |
|                    |given in the column header. For example if the column    |
|                    |header of the column you wish to plot is                 |
|                    |"From 100 - 500 magl" then any number from 100 to 500    |
|                    |inclusive will cause this column to be plotted. Note     |
|                    |however that if two levels are adjacent or overlap then  |
|                    |selecting a number in both ranges will cause both        |
|                    |columns to be plotted.                                   |
+--------------------+---------------------------------------------------------+
|z_leveltype         |If a z_level_list is provided it is also necessary to    |
|                    |provide a z_level type. There are currently four         |
|                    |different types of z level in NAME.                      |
|                    |                                                         |
|                    |   * altitude - height above sea level (m asl)           |
|                    |   * height - height above ground level (m agl)          |
|                    |   * z - text format height such as "vertical integral"  |
|                    |     or "boundary layer"                                 |
|                    |   * flight_level - flight levels                        |
+--------------------+---------------------------------------------------------+
|num_members         |Number of ensemble members to include on postage stamp   |
|                    |plots (current code only supports default value of 18).  |
+--------------------+---------------------------------------------------------+
|levels_list         |List of the contour levels to be used in the plot        |
+--------------------+---------------------------------------------------------+
|nlevels             |Number of contour levels if levels_list is not set.      |
|                    |Defaults to 7.                                           |
+--------------------+---------------------------------------------------------+
|back                |True or False (default is False). If True will allow ADAQ|
|                    |Python to correct the time bound information in the cube |
|                    |for a NAME back run                                      |
+--------------------+---------------------------------------------------------+
|extent_list         |Extent of the plot in longitude and latitude. Should be  |
|                    |a 4-digit list in the order: xmin, xmax, ymin, ymax.     |
+--------------------+---------------------------------------------------------+
|cmap                |Name of a colormap to use in plotting. Default is        |
|                    |'YlGnBu'                                                 |
+--------------------+---------------------------------------------------------+
|mapping             |Mapping background to add to plume plot. Options are     |
|                    |'states', 'countries' (default), 'coastline', 'grey_os'  |
|                    |(UK only), 'wms', 'jam' or 'None'.                       |
|                    |Note that when 'countries' is selected coastlines and the|
|                    |great lakes are plotted and when 'states' are selected   |
|                    |counties are included for zoomed in maps.                |
+--------------------+---------------------------------------------------------+
|projection          |Choice of map projection for output. Options are:        |
|                    |Mercator, PlateCarree (Default), SouthPolarStereo,       |
|                    |NorthPolarStereo, LambertIceland                         |
+--------------------+---------------------------------------------------------+
|plottype            |Select a style of plotting. Options are pcolormesh       |
|                    |(Default, same as pixel in IDL), contourf (filled        |
|                    |contours), contourf_edge (filled contours with a black   |
|                    |edge) or contour                                         |
+--------------------+---------------------------------------------------------+
|cbar_label          |As a default, the QUANTITY is used in the label on the   |
|                    |colourbar. The units are automatically appended in [].   |
|                    |An alternative label can be provided via this keyword.   |
|                    |Set the value to '' to remove the label completely.      |
+--------------------+---------------------------------------------------------+
|cbar_num_fmt        |The colourbar number format. Default for an AQ plot is   |
|                    |'None' and NAME is '%.1e' e - exponential, f - float     |
|                    |For example: '%4.2f' would produce a float with a total  |
|                    |of 4 digits (including the decimal place) with 2 of them |
|                    |after the decimal place, therefore giving  3.14          |
+--------------------+---------------------------------------------------------+
|cbar                |True or False. Show colourbar on postage stamp plots.    |
|                    |If requested, the colourbar is always plotted in an      |
|                    |upper-left position with a horizontal orientation.       |
|                    |Default is True.                                         |
+--------------------+---------------------------------------------------------+
|title               |A string to be used as a plot subtitle on the colourbar. |
|                    |If set to None, a default title will be generated based  |
|                    |on the quantity and species name.                        |
+--------------------+---------------------------------------------------------+
|suptitle            |A string to be used for the main title (super title).    |
|                    |If set to None, a default super title will be generated. |
+--------------------+---------------------------------------------------------+
|percent             |Computes percentage contributions of a horizontal field. |
|                    |Accepted values are True or False.                       |
|                    |Currently cannot be used on average quantities           |
+--------------------+---------------------------------------------------------+
|mobrand             |Add Met Office branding to plots (True or False (default)|
+--------------------+---------------------------------------------------------+
|units               |Option to convert units (e.g. units = 'ug/m3').          |
|                    |Currently applied to all fields requested                |
+--------------------+---------------------------------------------------------+
|figsize             |Sets the figure size. If not set, figsize will be set to |
|                    |the default for postage stamp plots [16, 11].            |
+--------------------+---------------------------------------------------------+
|plot_dir            |Full directory path of the location you would like your  |
|                    |plots to be saved to.                                    |
+--------------------+---------------------------------------------------------+
|short_name_list     |Name of the cube to be used used for calculating the     |
|                    |exceedance. For NAME output files the cube names are     |
|                    |constructed by concatenating the species and the quantity|
|                    |replacing all spaces with underscores and capitalising.  |
|                    |e.g. VOLCANIC_ASH_AIR_CONCENTRATION.  Note that it is    |
|                    |advisable to use only one due to memory constraints.     |
+--------------------+---------------------------------------------------------+

"""
from __future__ import print_function

import datetime
import os
import sys
import warnings
import iris
# Find adaq directory -
# should be in directory above the one containing this file
adaq_path = os.path.dirname(os.path.realpath(__file__))+'/../'
sys.path.append(adaq_path)
import adaqcode.adaq_functions
from adaqcode.adaq_plotting import plot_md_multi_fields
import adaqcode.inifile


def check_ini_dict(ini_dict):
    '''
    Function to check the ini_dict has the correct variables
    to run this script.

    :param ini_dict: Dictionary of the ini file.
    '''

    if not 'models_dir_list' in ini_dict:
        raise ValueError('"models_dir_list" must be set in the inifile.')

    if not 'models_fmt_list' in ini_dict:
        raise ValueError('"models_fmt_list" must be set in the inifile.')

    if not 'models_list' in ini_dict:
        raise ValueError('"models_list" must be set in the inifile.')

    if 'z_level_list' in ini_dict:
        if len(ini_dict['z_level_list']) > 3:
            warnings.warn('The system may not have enough memory for z_level_list length of '
                          +str(len(ini_dict['z_level_list'])))


def postage_stamp_plotting(inifilename=None, ini_dict=None, verbose=1):
    '''
    Top-level plotting routine to generate postage stamps.

    :param inifilename: String giving filename of ini file.
    :param verbose: Level of print output:

                    * 0 = No extra printing
                    * 1 = Standard level (recommended)
                    * 2 = Extra level for debugging

    :returns: (ini_dict, md_list)

     * **ini_dict** - Dictionary of a :class:`inifile` object
     * **md_list** - Model Data list - each item will be name
       in :class:`adaq_data` format

    **Method:**
     * Reads in inifile
     * Reads in model data (multiple ensemble members)
     * Plots postage stamps using style given in inifile

    >>> ini_dict, md_list = postage_stamp_plotting()
    ... # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    Current version of python: ...
    [...]
    Current version of iris: ...
    Beginning at  ...
    Reading inifile .../adaqscripts/name_postage_stamp_plot.ini
    Getting model data for  ...  at  ...
    <class '...NAMEData'> - Subclass of ADAQData - Contains:
    sites_cube_list:
    < No cubes >
    gridded_cube_list:
    0: CAESIUM-137_AIR_CONCENTRATION / (Bq / m^3) \
    (realization: 18; latitude: 300; longitude: 300)
    trajectory_cube_list:
    < No cubes >
    Plotting postage stamps
    Saved figure \
    .../Fieldplot_mogreps-g_name_CAESIUM-137_AIR_CONCENTRATION_Boundarylayeraverage\
_201908290200.png
    Finished at  ...
    '''

    print('Current version of python:', sys.version)
    print('Current version of iris:', iris.__version__)
    print('Beginning at ', datetime.datetime.now())

    # Read ini file
    if ini_dict is None:
        ini_dict = adaqcode.inifile.get_inidict(
            inifilename=inifilename,
            defaultfilename='adaqscripts/name_postage_stamp_plot.ini')

    check_ini_dict(ini_dict)

    md_list = adaqcode.adaq_functions.get_models(ini_dict)

    if 'percent' in ini_dict and ini_dict['percent']:
        md_list = adaqcode.adaq_functions.percent_contribution(md_list)

    # Check for presence of short_name_list and populate if not there
    if 'short_name_list' not in ini_dict:
        short_name_list = []
        for md in md_list:
            short_name_list.extend([cube.attributes['short_name']
                                    for cube in md.gridded_cube_list])
        ini_dict['short_name_list'] = short_name_list

    if len(ini_dict['short_name_list']) > 1:
        warnings.warn('Plotting '+str(len(ini_dict['short_name_list'])) +
                      ' fields may take some time')

    for md in md_list:
        print(md)
        if 'units' in ini_dict:
            for cube in md.gridded_cube_list:
                try:
                    cube.convert_units(ini_dict['units'])
                except ValueError:
                    print("Can't convert units from ", cube.units, \
                          " to ", ini_dict['units'])

    figsize = ini_dict.get('figsize', None)
    if figsize is not None:
        figsize = figsize
    else:
        figsize = [16, 11]

    # Plot data
    for short_name in ini_dict['short_name_list']:
        plot_md_multi_fields(ini_dict, md_list, short_name, defaults='NAME',
                             figsize=figsize, verbose=verbose)

    print('Finished at ', datetime.datetime.now())

    return ini_dict, md_list


if __name__ == '__main__':

    postage_stamp_plotting()

    # import doctest
    # doctest.testmod()
    # doctest.run_docstring_examples(get_ini, globals())
