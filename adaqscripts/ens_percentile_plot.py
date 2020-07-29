#!/usr/bin/env python
"""
Generalised script to plot of the percentile of the ensemble members.

Running the code
~~~~~~~~~~~~~~~~
To run the name field plotting code you need to create an inifile containing
information about the location of your name output and the columns you would
like to plot (see next section). The plots are then created by typing

.. code-block:: ksh

   ./ens_percentile_plot.py [inifilename]

where inifilename includes the full path to the ini file.

Alternatively, from within other python code, run using:

.. code-block:: python

    from adaqscripts import ens_percentile_plot
    ens_percentile_plot.percentile_plotting(inifilename='inifilename')

or if passing an ini_dict rather than an ini file

.. code-block:: python

    from adaqscripts import ens_percentile_plot
    # Set up ini_dict, e.g.
    ini_dict = {}
    ini_dict['key'] = value
    ens_percentile_plot.percentile_plotting(ini_dict=ini_dict)

Contents of the inifile
~~~~~~~~~~~~~~~~~~~~~~~

Essential Items
===============

The inifile needs to contain the following items:

+----------------+--------------------------------------------------------------+
|models_dir_list |List of the full directory paths of the data to be read in.   |
|                |Linux wildcards such as * and ? may be used in the path names.|
+----------------+--------------------------------------------------------------+
|models_fmt_list |List of the data formats to be read in. It should be          |
|                |"ens_name" for any esenmble output from NAME.                 |
+----------------+--------------------------------------------------------------+
|models_list     |List of names for the model. This can be any text and the     |
|                |text will be added to the filename. So, for example, if       |
|                |you have done one NAME run with convection on and one         |
|                |NAME run with convection off you could set the                |
|                |*models_list* to "no_convection" and "convection"             |
|                |respectively.                                                 |
+----------------+--------------------------------------------------------------+


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
|                    |'PuBuGn'                                                 |
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
|figsize             |Sets the figure size (Specific size e.g. [4,6],          |
|                    |'Default' (sets to [8,6]) or 'Auto' (default).           |
+--------------------+---------------------------------------------------------+
|plot_dir            |Full directory path of the location you would like your  |
|                    |plots to be saved to.                                    |
+--------------------+---------------------------------------------------------+
|units               |units (e.g. units = 'ug/m3')                             |
+--------------------+---------------------------------------------------------+
|percentile_list     |List of percentiles to calculate from the ensemble.      |
|                    |Default is: [10, 50, 90]                                 |
+--------------------+---------------------------------------------------------+
|short_name_list     |Name of the cube to be used used for calculating the     |
|                    |percentile. For NAME output files the cube names are     |
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
import numpy

# Find adaq directory -
# should be in directory above the one containing this file
adaq_path = os.path.dirname(os.path.realpath(__file__))+'/../'
sys.path.append(adaq_path)
import adaqcode.adaq_functions
from adaqcode.adaq_plotting import plot_md_gridded_fields
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

    # If colormap not in inifile, default to PuBuGn
    if 'cmap' not in ini_dict:
        ini_dict['cmap'] = 'PuBuGn'

def percentile_plotting(inifilename=None, ini_dict=None, verbose=1):
    '''
    Top-level plotting routine to generate percentile plots.


    :param inifilename: String giving filename of ini file.

    :param ini_dict: Dictionary of the ini file.

    :param verbose: Level of print output:

                    * 0 = No extra printing
                    * 1 = Standard level (recommended)
                    * 2 = Extra level for debugging

    :returns: (ini_dict, percs_md_list)

     * **ini_dict** - Dictionary of a :class:`inifile` object
     * **percs_md_list** - Model Data list - each item will be name
       in :class:`adaq_data` format

    **Method:**
     * Reads in inifile
     * Reads in model data
     * Calculates the percentile of the ensemble members
     * Plots data using style given in inifile

    >>> ini_dict, percs_md_list = percentile_plotting()
    ... # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    Current version of python: ...
    [...]
    Current version of iris: ...
    Beginning at  ...
    Reading inifile .../adaqscripts/ens_percentile_plot.ini
    Getting model data for  mogreps-g_name  at ...
    Plotting gridded fields
    Saved figure \
    .../Fieldplot_mogreps-g_name_95TH_PERCENTILE_OF_VOLCANIC_ASH_AIR\
_CONCENTRATION_FromFL100toFL125_201905141600.png
    ...
    Finished at  ...
    ...
    '''

    print(' ')
    print('Current version of python:', sys.version)
    print('Current version of iris:', iris.__version__)
    print('Beginning at ', datetime.datetime.now())

    if ini_dict is None:
        ini_dict = adaqcode.inifile.get_inidict(
            inifilename=inifilename,
            defaultfilename='adaqscripts/ens_percentile_plot.ini')

    check_ini_dict(ini_dict)

    percentiles = ini_dict.get('percentile_list', [10, 50, 90])
    percentiles = [float(v) for v in percentiles]

    md_list = adaqcode.adaq_functions.get_models(ini_dict)

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

    if 'units' in ini_dict:
        for md in md_list:
            for ngc in md.gridded_cube_list:
                try:
                    ngc.convert_units(ini_dict['units'])
                except ValueError:
                    print("Can't convert units from ", ngc.units, \
                          " to ", ini_dict['units'])

    percs_md_list = []
    for md in md_list:                               # List of models
        percs_md = adaqcode.adaq_data.ADAQData()
        for ngc in md.gridded_cube_list:             # List of variables
            ngc.coord('realization').guess_bounds()
            percs_ngc = ngc.collapsed(['realization'],
                                      iris.analysis.PERCENTILE,
                                      percent=percentiles,
                                      fast_percentile_method=True)
            for percentile in percentiles:
                perc_ngc = percs_ngc.extract(iris.Constraint
                                             (percentile_over_realization=
                                              percentile))
                perc_ngc.attributes['short_name'] = str(int(percentile)) + \
                    'TH_PERCENTILE_OF_' + ngc.attributes['short_name']
                perc_ngc.attributes['Quantity'] = str(int(percentile)) + \
                    'th percentile of ' + ngc.attributes['Quantity']
                perc_ngc.attributes['Percentile'] = percentile
                perc_ngc.data = numpy.ma.masked_equal(perc_ngc.data, 0)
                percs_md.gridded_cube_list.append(perc_ngc)

        percs_md_list.append(percs_md)

    del md_list  # Just to free up memory on the system.

    short_name_to_plot_list = []
    for percs_md in percs_md_list:
        short_name_to_plot_list.extend([cube.attributes['short_name']
                                        for cube in percs_md.gridded_cube_list])

    for short_name in short_name_to_plot_list:
        plot_md_gridded_fields(ini_dict,
                               percs_md_list,
                               short_name,
                               defaults='NAME',
                               verbose=verbose)

    print(' ')
    print('Finished at ', datetime.datetime.now())

    return ini_dict, percs_md_list


if __name__ == '__main__':

    percentile_plotting()

    # import doctest
    # doctest.testmod()
    # doctest.run_docstring_examples(get_ini, globals())
