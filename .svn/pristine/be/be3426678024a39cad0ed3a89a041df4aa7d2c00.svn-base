#!/usr/bin/env python
"""
Generalised script to generate field plots from NAME gridded files

Running the code
~~~~~~~~~~~~~~~~
To run the name field plotting code you need to create an inifile containing
information about the location of your name output and the columns you would
like to plot (see next section). The plots are then created by typing

.. code-block:: ksh

   ./name_field_plot.py [inifilename]

where inifilename includes the full path to the ini file.

Alternatively, from within other python code, run using:

.. code-block:: python

    from adaqscripts import name_field_plot
    name_field_plot.field_plotting(inifilename='inifilename')

or if passing an ini_dict rather than an ini file

.. code-block:: python

    from adaqscripts import name_field_plot
    # Set up ini_dict, e.g.
    ini_dict = {}
    ini_dict['key'] = value
    name_field_plot.field_plotting(ini_dict=ini_dict)

Contents of the inifile
~~~~~~~~~~~~~~~~~~~~~~~

Essential Items
===============

The inifile needs to contain the following items:

+---------------+--------------------------------------------------------------+
|models_dir_list|List of the full directory paths of the data to be read in.   |
|               |Linux wildcards such as * and ? may be used in the path names.|
+---------------+--------------------------------------------------------------+
|models_fmt_list|List of the data formats to be read in. It should be "name"   |
|               |for any output from NAME.                                     |
+---------------+--------------------------------------------------------------+
|plot_dir       |Full directory path of the location you would like your plots |
|               |to be saved to.                                               |
+---------------+--------------------------------------------------------------+
|short_name_list|List of the names of the cubes to be plotted. For NAME output |
|               |files the cube names are constructed by concatenating the     |
|               |species and the quantity, replacing all spaces with           |
|               |underscores and capitalising.                                 |
|               |e.g. VOLCANIC_ASH_AIR_CONCENTRATION.                          |
|               |Remove short_name_list to plot all columns in the file        |
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
|models_list         |List of names for the model. This can be any text and the|
|                    |text will be added to the filename. So, for example, if  |
|                    |you have done one NAME run with convection on and one    |
|                    |NAME run with convection off you could set the           |
|                    |*models_list* to "no_convection" and "convection"        |
|                    |respectively.                                            |
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
|levels_list         |List of the contour levels to be used in the plot        |
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
|cbar_orientaion     |The default for an AQ plot is 'horizontal'. The default  |
|                    |for a NAME plot is 'vertical'.                           |
+--------------------+---------------------------------------------------------+
|cbar_label          |As a default the short name (SPECIES_QUANTITY) is used   |
|                    |as a label on the colourbar. This keyword allows an      |
|                    |alternative label to be provided. To remove the colourbar|
|                    |label completely set the value of this keyword to ''.    |
+--------------------+---------------------------------------------------------+
|cbar_num_fmt        |The colourbar number format. Default for an AQ plot is   |
|                    |'None' and NAME is '%.1e' e - exponential, f - float     |
|                    |For example: '%4.2f' would produce a float with a total  |
|                    |of 4 digits (including the decimal place) with 2 of them |
|                    |after the decimal place, therefore giving  3.14          |
+--------------------+---------------------------------------------------------+
|cbar                |True or False. False will result in plots without        |
|                    |colorbars. Default is True                               |
+--------------------+---------------------------------------------------------+
|title               |A string which will be used as the plot title or the     |
|                    |word 'name_verbose' to request a title constructed from  |
|                    |the column headers in NAME. If set to None will use the  |
|                    |default title 'Valid at <date>'                          |
+--------------------+---------------------------------------------------------+
|suptitle            |A string to be used for the suptitle (or super title)    |
|                    |for all plots                                            |
+--------------------+---------------------------------------------------------+
|annote              |Text to be added to the side or below the plot. Can be   |
|                    |free text or 'NAME_defaults' which will add details of   |
|                    |the release such as the source location and release rate |
|                    |to the page.                                             |
+--------------------+---------------------------------------------------------+
|annote_location     |Select 'right' to put the annotation to the right of     |
|                    |the plot and 'below' to put the annotation below the     |
|                    |plot                                                     |
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
|figsize             |Sets the figure size (Specific size e.g. [4,6],          |
|                    |'Default' (sets to [8,6]) or 'Auto' (default).           |
+--------------------+---------------------------------------------------------+
"""
from __future__ import print_function

import datetime
import os
import sys
import iris
# Find adaq directory -
# should be in directory above the one containing this file
adaq_path = os.path.dirname(os.path.realpath(__file__))+'/../'
sys.path.append(adaq_path)
import adaqcode.adaq_functions
from adaqcode.adaq_plotting import plot_md_gridded_fields
import adaqcode.inifile


def field_plotting(inifilename=None, ini_dict=None, verbose=1):
    '''
    Top-level plotting routine to generate field plots

    :param inifilename: String giving filename of ini file.
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

    >>> ini_dict, md_list = field_plotting()
    ... # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    Current version of python: ...
    [...]
    Current version of iris: ...
    Beginning at  ...
    Reading inifile .../adaqscripts/name_field_plot.ini
    Getting model data for  name  at ...
    <class '...NAMEData'> - Subclass of ADAQData - Contains:
    sites_cube_list:
    < No cubes >
    gridded_cube_list:
    0: VOLCANIC_ASH_AIR_CONCENTRATION / (g/m3) (latitude: 480; longitude: 640)
    trajectory_cube_list:
    < No cubes >
    Plotting gridded fields
    Saved figure \
    .../Fieldplot_name_VOLCANIC_ASH_AIR_CONCENTRATION_FromFL025to\
FL050_201105230000.png
    Finished at  ...

    '''

    print('Current version of python:', sys.version)
    print('Current version of iris:', iris.__version__)
    print('Beginning at ', datetime.datetime.now())

    # Read ini file
    if ini_dict is None:
        ini_dict = adaqcode.inifile.get_inidict(
            inifilename=inifilename,
            defaultfilename='adaqscripts/name_field_plot.ini')

    # Read in model data
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

    for md in md_list:
        print(md)
    if 'units' in ini_dict:
        for ngc in md.gridded_cube_list:
            try:
                ngc.convert_units(ini_dict['units'])
            except ValueError:
                print("Can't convert units from ", ngc.units, \
                      " to ", ini_dict['units'])

    # Plot data
    for short_name in ini_dict['short_name_list']:
        plot_md_gridded_fields(ini_dict, md_list, short_name,
                               defaults='NAME', verbose=verbose)

    print('Finished at ', datetime.datetime.now())

    return ini_dict, md_list


if __name__ == '__main__':

    field_plotting()

    # import doctest
    # doctest.testmod()
    # doctest.run_docstring_examples(get_ini, globals())
