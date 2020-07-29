#!/usr/bin/env python
"""
Generalised script to compare two gridded fields of NAME output.

Running the code
~~~~~~~~~~~~~~~~
To run the name field plotting code you need to create an inifile containing
information about the locations of your name output and the columns you would
like to plot (see next section). The plots are then created by typing

.. code-block:: ksh

   ./name_field_plot_diff.py [inifilename]

where inifilename includes the full path to the ini file.

Alternatively, from within other python code, run using:

.. code-block:: python

    from adaqscripts import name_field_plot_diff
    name_field_plot_diff.field_diff_plotting(inifilename='inifilename')

or if passing an ini_dict rather than an ini file

.. code-block:: python

    from adaqscripts import name_field_plot_diff
    # Set up ini_dict, e.g.
    ini_dict = {}
    ini_dict['key'] = value
    name_field_plot_diff.field_diff_plotting(ini_dict=ini_dict)

Contents of the inifile
~~~~~~~~~~~~~~~~~~~~~~~

Essential Items
===============

The inifile needs to contain the following items:

+-----------------------+-------------------------------------------------------------------------+
|models_dir_case1_list  |List of the full directory paths of the case1 data to be read in.        |
|                       |Linux wildcards such as * and ? may be used in the path names.           |
+-----------------------+-------------------------------------------------------------------------+
|models_dir_case2_list  |List of the full directory paths of the Know Good Output data to be read |
|                       |in. Linux wildcards such as * and ? may be used in the path names.       |
+-----------------------+-------------------------------------------------------------------------+
|models_fmt_list        |List of the data formats to be read in. It should be "name"              |
|                       |for any output from NAME.                                                |
+-----------------------+-------------------------------------------------------------------------+
|plot_dir_montage       |Full directory path of the location you would like your Montage of case1 |
|                       |case2, difference, and relative difference plots.                        |
+-----------------------+-------------------------------------------------------------------------+
|models_list            |List of names for the model. This can be any text and the                |
|                       |text will be added to the filename. So, for example, if                  |
|                       |you have done one NAME run with convection on and one                    |
|                       |NAME run with convection off you could set the                           |
|                       |*models_list* to "no_convection" and "convection"                        |
|                       |respectively.                                                            |
+-----------------------+-------------------------------------------------------------------------+
|suptitle_case1         |A string to be used for the suptitle (or super title) for case 1.        |
+-----------------------+-------------------------------------------------------------------------+
|suptitle_case2         |A string to be used for the suptitle (or super title) for case 2.        |
+-----------------------+-------------------------------------------------------------------------+


Optional Items
==============

The inifile may also contain any of the following items:


+---------------------+---------------------------------------------------------+
|plot_dir_case1       |Full directory path of the location you would like your  |
|                     |case1 plots to be saved to.                              |
+---------------------+---------------------------------------------------------+
|plot_dir_case2       |Full directory path of the location you would like your  |
|                     |case2 plots to be saved to.                              |
+---------------------+---------------------------------------------------------+
|plot_dir_diff        |Full directory path of the location you would like your  |
|                     |Differnce plots to be saved to.                          |
+---------------------+---------------------------------------------------------+
|plot_dir_reldiff     |Full directory path of the location you would like your  |
|                     |Relative Difference plots to be saved to.                |
+---------------------+---------------------------------------------------------+
|field_attribute_dict |Dictionary of the field attributes of the NAME data to   |
|                     |be plotted. The field attributes are basically extracted |
|                     |from the column headers of the NAME file. See            |
|                     |:ref:`field_attributes_ref` for a more comprehensive list|
+---------------------+---------------------------------------------------------+
|z_level_list         |List of the heights which you would like to plot.        |
|                     |For NAME III format files this is the height given in    |
|                     |the column header  without units. For NAME II format     |
|                     |files you can provide any number between the two numbers |
|                     |given in the column header. For example if the column    |
|                     |header of the column you wish to plot is                 |
|                     |"From 100 - 500 magl" then any number from 100 to 500    |
|                     |inclusive will cause this column to be plotted. Note     |
|                     |however that if two levels are adjacent or overlap then  |
|                     |selecting a number in both ranges will cause both        |
|                     |columns to be plotted.                                   |
+---------------------+---------------------------------------------------------+
|z_leveltype          |If a z_level_list is provided it is also necessary to    |
|                     |provide a z_level type. There are currently four         |
|                     |different types of z level in NAME.                      |
|                     |                                                         |
|                     |   * altitude - height above sea level (m asl)           |
|                     |   * height - height above ground level (m agl)          |
|                     |   * z - text format height such as "vertical integral"  |
|                     |     or "boundary layer"                                 |
|                     |   * flight_level - flight levels                        |
+---------------------+---------------------------------------------------------+
|levels_list          |List of the contour levels to be used in the case1 and   |
|                     |case2 plots.                                             |
+---------------------+---------------------------------------------------------+
|abs_diff_levels_list |List of the contour levels to be used in the absolute    |
|                     |difference plots.                                        |
+---------------------+---------------------------------------------------------+
|rel_diff_levels_list |List of the contour levels to be used in the relative    |
|                     |difference plots.                                        |
+---------------------+---------------------------------------------------------+
|back                 |True or False (default is False). If True will allow ADAQ|
|                     |Python to correct the time bound information in the cube |
|                     |for a NAME back run                                      |
+---------------------+---------------------------------------------------------+
|extent_list          |Extent of the plot in longitude and latitude. Should be  |
|                     |a 4-digit list in the order: xmin, xmax, ymin, ymax.     |
+---------------------+---------------------------------------------------------+
|cmap                 |Name of a colormap to use in plotting the case1 and case2|
|                     |plots. The Default is 'YlGnBu'                           |
+---------------------+---------------------------------------------------------+
|mapping              |Mapping background to add to plume plot. Options are     |
|                     |'states', 'countries' (default), 'coastline', 'grey_os'  |
|                     |(UK only), 'wms', 'jam' or 'None'.                       |
|                     |Note that when 'countries' is selected coastlines and the|
|                     |great lakes are plotted and when 'states' are selected   |
|                     |counties are included for zoomed in maps.                |
+---------------------+---------------------------------------------------------+
|projection           |Choice of map projection for output. Options are:        |
|                     |Mercator, PlateCarree (Default), SouthPolarStereo,       |
|                     |NorthPolarStereo, LambertIceland                         |
+---------------------+---------------------------------------------------------+
|plottype             |Select a style of plotting. Options are pcolormesh       |
|                     |(Default, same as pixel in IDL), contourf (filled        |
|                     |contours), contourf_edge (filled contours with a black   |
|                     |edge) or contour                                         |
+---------------------+---------------------------------------------------------+
|cbar_orientaion      |The default for an AQ plot is 'horizontal'. The default  |
|                     |for a NAME plot is 'vertical'.                           |
+---------------------+---------------------------------------------------------+
|cbar_label           |As a default the short name (SPECIES_QUANTITY) is used   |
|                     |as a label on the colourbar. This keyword allows an      |
|                     |alternative label to be provided. To remove the colourbar|
|                     |label completely set the value of this keyword to ''.    |
+---------------------+---------------------------------------------------------+
|cbar_num_fmt         |The colourbar number format. Default for an AQ plot is   |
|                     |'None' and NAME is '%.1e' e - exponential, f - float     |
|                     |For example: '%4.2f' would produce a float with a total  |
|                     |of 4 digits (including the decimal place) with 2 of them |
|                     |after the decimal place, therefore giving  3.14          |
+---------------------+---------------------------------------------------------+
|cbar                 |True or False. False will result in plots without        |
|                     |colorbars. Default is True                               |
+---------------------+---------------------------------------------------------+
|title                |A string which will be used as the plot title or the     |
|                     |word 'name_verbose' to request a title constructed from  |
|                     |the column headers in NAME. If set to None will use the  |
|                     |default title 'Valid at <date>'                          |
+---------------------+---------------------------------------------------------+
|suptitle             |Cannot be used with name_field_plot_diff.py              |
+---------------------+---------------------------------------------------------+
|annote               |Text to be added to the side or below the plot. Can be   |
|                     |free text or 'NAME_defaults' which will add details of   |
|                     |the release such as the source location and release rate |
|                     |to the page.                                             |
+---------------------+---------------------------------------------------------+
|annote_location      |Select 'right' to put the annotation to the right of     |
|                     |the plot and 'below' to put the annotation below the     |
|                     |plot                                                     |
+---------------------+---------------------------------------------------------+
|percent              |Computes percentage contributions of a horizontal field. |
|                     |Accepted values are True or False.                       |
|                     |Currently cannot be used on average quantities           |
+---------------------+---------------------------------------------------------+
|mobrand              |Add Met Office branding to plots (True or False (default)|
+---------------------+---------------------------------------------------------+
|units                |Option to convert units (e.g. units = 'ug/m3').          |
|                     |Currently applied to all fields requested                |
+---------------------+---------------------------------------------------------+
|figsize              |Sets the figure size (Specific size e.g. [4,6],          |
|                     |'Default' (sets to [8,6]) or 'Auto' (default).           |
+---------------------+---------------------------------------------------------+
|short_name_list      |List of the names of the cubes to be plotted. For NAME   |
|                     |output files the cube names are constructed by           |
|                     |concatenating the species and the quantity, replacing    |
|                     |all spaces with underscores and capitalising.            |
|                     |e.g. VOLCANIC_ASH_AIR_CONCENTRATION.                     |
|                     |Remove short_name_list to plot all columns in the file   |
+---------------------+---------------------------------------------------------+
|fail_if_diff         |If True the script will return a positive return code.   |
|                     |The default is false.                                    |
+---------------------+---------------------------------------------------------+
|reldiff_tolerance    |The relative difference tolerance to allow. Any relative |
|                     |difference greater than the value will be flagged as a   |
|                     |difference. If "fail_if_diff" is to True than the script |
|                     |will return a positive return code.                      |
+---------------------+---------------------------------------------------------+


"""
from __future__ import print_function
from __future__ import division

import datetime
import os
import sys
import copy
import math
import iris
import numpy as np

# Find adaq directory -
# should be in directory above the one containing this file
adaq_path = os.path.dirname(os.path.realpath(__file__)) + '/../'
sys.path.append(adaq_path)
from adaqcode import shell_commands, adaq_functions, inifile
from adaqcode.adaq_plotting import plot_md_gridded_fields


def check_ini_dict(ini_dict):
    '''
    Function to check the ini_dict has the correct variables
    to run this script.

    :param ini_dict: Dictionary of the ini file.
    '''

    if not 'models_dir_case1_list' in ini_dict:
        raise ValueError('"models_dir_case1_list" must be set in the inifile.')

    if not 'models_dir_case2_list' in ini_dict:
        raise ValueError('"models_dir_case2_list" must be set in the inifile.')

    if not 'models_fmt_list' in ini_dict:
        raise ValueError('"models_fmt_list" must be set in the inifile.')

    if not 'models_list' in ini_dict:
        raise ValueError('"models_list" must be set in the inifile.')

    if 'plot_dir_montage' in ini_dict:
        if not os.path.exists(ini_dict['plot_dir_montage']):
            os.makedirs(ini_dict['plot_dir_montage'])
    else:
        raise ValueError('"plot_dir_montage" must be set in the inifile.')

    if not 'suptitle_case1' in ini_dict:
        raise ValueError('"suptitle_case1" must be set in the inifile.')

    if not 'suptitle_case2' in ini_dict:
        raise ValueError('"suptitle_case2" must be set in the inifile.')


def calc_symmet_log10_levels(max_abs_val, min_abs_val, n_levels):
    '''
    Function to calculate symmetrical logarithmic levels centered on zero,
    from the maximium absolute value and the number of levels desired.

    :param max_abs_val: float of the maximium absolute value for the levels
    :param min_abs_val: float of the minimium absolute value for the levels
    :param n_levels:    integer number of levels (Must be an even number)

    :returns: (levels_log10_list)

     * **levels_log10_list** - List of the logarithmic levels

    >>> levels_log10_list = calc_symmet_log10_levels(2e36, 1.0, 6)
    >>> levels_log10_list
    [-2.0000000000000004e+36, -1.4142135623730952e+18, -1.0, \
1.0, 1.4142135623730952e+18, 2.0000000000000004e+36]
    >>> levels_log10_list = calc_symmet_log10_levels(1.0e4, 10, 8)
    >>> levels_log10_list
    [-10000.0, -1000.0, -100.0, -10.0, 10.0, 100.0, 1000.0, 10000.0]
    >>> levels_log10_list = calc_symmet_log10_levels(7.0e-37, 7.0e-40, 6)
    >>> levels_log10_list
    [-7.000000000000001e-37, -2.213594362117866e-38, -7.000000000000001e-40, \
7.000000000000001e-40, 2.213594362117866e-38, 7.000000000000001e-37]
    '''

    if (n_levels%2) == 1:
        raise ValueError('The number of levels must be even')

    if min_abs_val >= max_abs_val:
        raise ValueError(' min_abs_val needs to be less than max_abs_val. '+
                         ' min_abs_val: '+str(min_abs_val)+
                         ' max_abs_val: '+str(max_abs_val))

    n_levels_top_half = int(n_levels / 2)
    np_levels_log10_tophalf = np.logspace(math.log10(min_abs_val),
                                          math.log10(max_abs_val),
                                          n_levels_top_half)

    np_levels_log10_bothalf = - np.flipud(np_levels_log10_tophalf)

    np_levels_log10 = np.concatenate([np_levels_log10_bothalf,
                                      np_levels_log10_tophalf])

    levels_log10_list = np_levels_log10.tolist()

    return levels_log10_list


def calc_symmet_linear_levels(max_abs_val, n_levels):
    '''
    Function to calculate symmetrical linear levels centered on zero, from
    the maximium absolute value and the number of levels desired.

    :param max_abs_val: float of the maximium absolute value for the levels
    :param n_levels:    integer number of levels

    :returns: (levels_list)

     * **levels_list** - List of the levels

    >>> levels_list = calc_symmet_linear_levels(2e36, 6)
    >>> levels_list
    [-2e+36, -1.2000000000000001e+36, -4.000000000000001e+35, \
3.9999999999999984e+35, 1.1999999999999998e+36, 2e+36]
    >>> levels_list = calc_symmet_linear_levels(2e-36, 6)
    >>> levels_list
    [-2e-36, -1.2e-36, -4.000000000000001e-37, \
3.999999999999998e-37, 1.1999999999999997e-36, 2e-36]
    '''

    np_levels = np.linspace(-max_abs_val, max_abs_val, n_levels)
    levels_list = np_levels.tolist()

    return levels_list


def differencing_cubes(ini_dict, md_list_case1, md_list_case2):
    '''
    To calulate the difference between two sets of model data lists.

    Diff = case2 - case1
    and
    Relative Diff = Diff * 100 / (case2 + case1)

    :param md_list_case1:  Model Data list - each item will be name in
                           :class:`adaq_data` format

    :param md_list_case2:  Model Data list - each item will be name in
                           :class:`adaq_data` format

    :param ini_dict: Dictionary of properties required for plot type.

    :returns: (diff_count, abs_diff_levels, rel_diff_levels, md_list_diff, md_list_reldiff)

     * **diff_count**      Interger number of cubes which have a difference.
     * **abs_diff_levels** Levels to be used for the difference plotting.
     * **rel_diff_levels**      Levels to be used for the relative difference plotting.
     * **md_list_diff**    Model Data list - of the difference.
     * **md_list_reldiff** Model Data list - of the relative difference.
    '''

    # Options for the script which are not in the ini file.
    n_levels = 6

    # initalise output md_lists
    md_list_diff = copy.deepcopy(md_list_case2)
    md_list_reldiff = copy.deepcopy(md_list_case2)

    diff_count = 0
    max_abs_diff_list = []
    min_abs_diff_exec_zeros_list = []
    max_abs_reldiff_list = []

    for md_case1, md_case2, md_diff, md_reldiff in zip(md_list_case1,
                                                       md_list_case2,
                                                       md_list_diff,
                                                       md_list_reldiff):

        if len(md_case1.gridded_cube_list) != len(md_case2.gridded_cube_list):
            print('md_case1.gridded_cube_list: ', md_case1.gridded_cube_list)
            print('md_case2.gridded_cube_list: ', md_case2.gridded_cube_list)
            raise ValueError('Cannot diff outputs as: \
                             len(md_case1.gridded_cube_list): ',
                             len(md_case1.gridded_cube_list),
                             'len(md_case2.gridded_cube_list): ',
                             len(md_case2.gridded_cube_list))

        for ngc_case1, ngc_case2, ngc_diff, ngc_reldiff in zip(
                md_case1.gridded_cube_list,
                md_case2.gridded_cube_list,
                md_diff.gridded_cube_list,
                md_reldiff.gridded_cube_list):

            if ngc_case1.data.shape != ngc_case2.data.shape:
                print('models_dir_case1_list: ',
                      ini_dict['models_dir_case1_list'])
                print('models_dir_case2_list: ',
                      ini_dict['models_dir_case2_list'])
                raise ValueError('Cannot diff outputs as: ',
                                 'ngc_case1.data.shape: ', ngc_case1.data.shape,
                                 'ngc_case2.data.shape: ', ngc_case2.data.shape)

            ngc_diff.data = ngc_case2.data - ngc_case1.data
            ngc_reldiff.data[:] = 0.0
            ngc_reldiff.units = '%'
            mask1 = ngc_case2.data + ngc_case1.data != 0
            ngc_reldiff.data[mask1] = ngc_diff.data[mask1] * 100.0 / \
                (ngc_case2.data[mask1] + ngc_case1.data[mask1])
            mask2 = ngc_diff.data != 0
            max_abs_diff = np.max(abs(ngc_diff.data))
            min_abs_diff_exec_zeros = np.min(abs(ngc_diff.data[mask2]))
            max_abs_reldiff = np.max(abs(ngc_reldiff.data))
            max_abs_diff_list.append(max_abs_diff)
            min_abs_diff_exec_zeros_list.append(min_abs_diff_exec_zeros)
            max_abs_reldiff_list.append(max_abs_reldiff)
            ngc_diff.data[
                ngc_diff.data == 0] = np.nan  # To blank zeros from plot
            ngc_reldiff.data[
                ngc_reldiff.data == 0] = np.nan  # To blank zeros from plot

            if max_abs_reldiff > float(ini_dict.get('reldiff_tolerance', 0.0)):
                diff_count += 1
                print('*** The two cubes are different ***')
                print('max_abs_reldiff:   ', max_abs_reldiff)
                print('reldiff_tolerance: ',
                      ini_dict.get('reldiff_tolerance', 0.0))
                print('***********************************')
            else:
                print('*** The two cubes are the same ***')
                print('max_abs_reldiff:   ', max_abs_reldiff)
                print('reldiff_tolerance: ',
                      ini_dict.get('reldiff_tolerance', 0.0))
                print('***********************************')

    if diff_count > 0:
        abs_diff_levels = calc_symmet_log10_levels(max(max_abs_diff_list),
                                                   min(min_abs_diff_exec_zeros_list),
                                                   n_levels)
        rel_diff_levels = calc_symmet_linear_levels(max(max_abs_reldiff_list),
                                                    n_levels)
    else:
        abs_diff_levels = [np.nan]
        rel_diff_levels = [np.nan]

    return diff_count, abs_diff_levels, rel_diff_levels, md_list_diff, md_list_reldiff


def check_shape(md_list_case1, md_list_case2):
    '''
    Checking two model data lists have the same length. Raising an exception
    with diagnostic information if there is a difference.

    :param md_list_case1:  Model Data list - each item will be name in
                           :class:`adaq_data` format

    :param md_list_case2:  Model Data list - each item will be name in
                           :class:`adaq_data` format
    '''

    if len(md_list_case1) != len(md_list_case2):
        raise ValueError('Cannot diff outputs as: \
                          len(md_list_case1): ', len(md_list_case1), \
                         'len(md_list_case2): ', len(md_list_case2))


def convert_units(md_list, ini_dict):
    '''
    Function to convert the units of all cubes within a model data list.

    :param md_list:  Model Data list - each item will be name in :class:`adaq_data` format
    :param ini_dict: Dictionary of properties required for plot type.

    :returns: (md_list)

     * **md_list** - Model Data list - each item will be name in :class:`adaq_data` format
    '''

    for md in md_list:
        for ngc in md.gridded_cube_list:
            try:
                ngc.convert_units(ini_dict['units'])
            except ValueError:
                print("Can't convert units from ", ngc.units, \
                      " to ", ini_dict['units'])

    return md_list


def field_diff_plotting(inifilename=None, ini_dict=None, verbose=1):
    '''
    Top-level plotting routine to generate field plots of both
    case1 and case2 together with difference and relative difference plots.

    :param inifilename: String giving filename of ini file.

    :param ini_dict: Dictionary of properties required for plot type.

    :param verbose: Level of print output:

                    * 0 = No extra printing
                    * 1 = Standard level (recommended)
                    * 2 = Extra level for debugging

    :returns: (ini_dict)

     * **ini_dict** - Dictionary of a :class:`inifile` object

    **Method:**
     * Reads in inifile, case1 model data and case2 model data
     * Calculates the differences between the two.
     * Plots both case1 and case 2 with the difference and relative difference in a montage.

    >>> ini_dict = field_diff_plotting()
    ... # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    Current version of python: ...
    [...]
    Current version of iris: ...
    Beginning at  ...
    Reading inifile .../adaqscripts/name_field_plot_diff.ini
    Getting model data for  global-name  at ...
    Getting model data for  global-name  at ...
    *** The two cubes are different ***
    max_abs_reldiff:    99.98737
    reldiff_tolerance:  5
    ***********************************
    Plotting gridded fields
    ...
    Saved figure .../case1/Fieldplot_global-name_NO2_EULERIAN_CONCENTRATION\
_From0_0to100_0magl_201506140000.png
    Plotting gridded fields
    ...
    Saved figure .../case2/Fieldplot_global-name_NO2_EULERIAN_CONCENTRATION\
_From0_0to100_0magl_201506140000.png
    Plotting gridded fields
    ...
    Saved figure .../diff/Fieldplot_global-name_NO2_EULERIAN_CONCENTRATION\
_From0_0to100_0magl_201506140000.png
    Plotting gridded fields
    ...
    Saved figure .../reldiff/Fieldplot_global-name_NO2_EULERIAN_CONCENTRATION\
_From0_0to100_0magl_201506140000.png
    Plotting montage
    ...
    Saved figure .../montage/Fieldplot_global-name_NO2_EULERIAN_CONCENTRATION\
_From0_0to100_0magl_201506140000.png
    Finished at  ...
    '''

    print('Current version of python:', sys.version)
    print('Current version of iris:', iris.__version__)
    print('Beginning at ', datetime.datetime.now())

    # Read ini file and check it
    if ini_dict is None:
        ini_dict = inifile.get_inidict(inifilename=inifilename,
                                       defaultfilename=
                                       'adaqscripts/name_field_plot_diff.ini')

    check_ini_dict(ini_dict)

    # Get Data
    try:
        ini_dict['models_dir_list'] = ini_dict['models_dir_case1_list']
    except:
        raise ValueError("models_dir_case1_list has not be set in inifile")

    md_list_case1 = adaq_functions.get_models(ini_dict)

    try:
        ini_dict['models_dir_list'] = ini_dict['models_dir_case2_list']
    except:
        raise ValueError("models_dir_case2_list has not be set in inifile")

    md_list_case2 = adaq_functions.get_models(ini_dict)

    # Change the units of both md_lists
    if 'units' in ini_dict:
        md_list_case1 = convert_units(md_list_case1, ini_dict)
        md_list_case2 = convert_units(md_list_case2, ini_dict)

    # Check for presence of short_name_list and populate if not there
    if 'short_name_list' not in ini_dict:
        short_name_list = []
        for md_case1 in md_list_case1:
            short_name_list.extend([cube.attributes['short_name']
                                    for cube in md_case1.gridded_cube_list])
        print('short_name_list: ', short_name_list)
        ini_dict['short_name_list'] = short_name_list

    check_shape(md_list_case1, md_list_case2)
    diff_count, abs_diff_levels, rel_diff_levels, md_list_diff, md_list_reldiff = \
        differencing_cubes(ini_dict, md_list_case1, md_list_case2)

    if diff_count > 0:
        for short_name in ini_dict['short_name_list']:

            # Get directories ready
            if 'plot_dir_case1' not in ini_dict:
                ini_dict['plot_dir_case1'] = '/tmp/plot_dir_case1'

            if 'plot_dir_case2' not in ini_dict:
                ini_dict['plot_dir_case2'] = '/tmp/plot_dir_case2'

            if 'plot_dir_diff' not in ini_dict:
                ini_dict['plot_dir_diff'] = '/tmp/plot_dir_diff'

            if 'plot_dir_reldiff' not in ini_dict:
                ini_dict['plot_dir_reldiff'] = '/tmp/plot_dir_reldiff'

            shell_commands.call_shell(
                'rm -f ' + ini_dict['plot_dir_case1'] + '/*')
            shell_commands.call_shell(
                'rm -f ' + ini_dict['plot_dir_case2'] + '/*')
            shell_commands.call_shell(
                'rm -f ' + ini_dict['plot_dir_diff'] + '/*')
            shell_commands.call_shell(
                'rm -f ' + ini_dict['plot_dir_reldiff'] + '/*')

            # Get cbar_labels ready
            if 'cbar_label' in ini_dict:
                cbar_label_diff = 'Absolute Difference in ' + ini_dict[
                    'cbar_label']
                cbar_label_reldiff = 'Relative Difference in ' + ini_dict[
                    'cbar_label']
            else:
                cbar_label_diff = 'Absolute Difference in ' + short_name
                cbar_label_reldiff = 'Relative Difference in ' + short_name

            # Plot the four graphics individually
            # Case 1
            ini_dict['plot_dir'] = ini_dict['plot_dir_case1']
            ini_dict['suptitle'] = 'case1: ' + ini_dict['suptitle_case1']
            plot_md_gridded_fields(ini_dict, md_list_case1, short_name,
                                   defaults='NAME', verbose=verbose)

            # Case 2
            ini_dict['plot_dir'] = ini_dict['plot_dir_case2']
            ini_dict['suptitle'] = 'case2: ' + ini_dict['suptitle_case2']
            plot_md_gridded_fields(ini_dict, md_list_case2, short_name,
                                   defaults='NAME', verbose=verbose)

            # Abs Diff
            if 'abs_diff_levels_list' in ini_dict:
                ini_dict['levels_list'] = ini_dict['abs_diff_levels_list']
            else:
                ini_dict['levels_list'] = abs_diff_levels

            ini_dict['plot_dir'] = ini_dict['plot_dir_diff']
            ini_dict['suptitle'] = 'Absolute Difference = case2 - case1'
            ini_dict['cmap'] = 'PiYG'
            ini_dict['cbar_label'] = cbar_label_diff
            plot_md_gridded_fields(ini_dict, md_list_diff, short_name,
                                   defaults='NAME', verbose=verbose)

            # Rel Diff
            if 'rel_diff_levels_list' in ini_dict:
                ini_dict['levels_list'] = ini_dict['rel_diff_levels_list']
            else:
                ini_dict['levels_list'] = rel_diff_levels

            ini_dict['plot_dir'] = ini_dict['plot_dir_reldiff']
            ini_dict['suptitle'] = 'Relative Difference = ' \
                                   '(case2-case1)*100/(case2+case1)'
            ini_dict['cmap'] = 'bwr'
            ini_dict['cbar_label'] = cbar_label_reldiff
            plot_md_gridded_fields(ini_dict, md_list_reldiff, short_name,
                                   defaults='NAME', verbose=verbose)

            # Create the montage
            print('Plotting montage')
            rc = shell_commands.call_shell("for file in $(ls " +
                                           ini_dict['plot_dir_diff'] +
                                           ") ; do montage " +
                                           ini_dict['plot_dir_case1'] +
                                           "/$file " +
                                           ini_dict['plot_dir_case2'] +
                                           "/$file " +
                                           ini_dict['plot_dir_diff'] +
                                           "/$file " +
                                           ini_dict['plot_dir_reldiff'] +
                                           "/$file -geometry +2+2 " +
                                           ini_dict['plot_dir_montage'] +
                                           "/$file ; "
                                           "echo Saved figure " +
                                           ini_dict['plot_dir_montage'] +
                                           "/$file ; done")

            if rc != 0:
                raise ValueError("Error with script to produce montage")

    print('Finished at ', datetime.datetime.now())

    if diff_count > 0 and ini_dict.get('fail_if_diff', False):
        sys.exit(1)

    return ini_dict


if __name__ == '__main__':
    field_diff_plotting()

    # import doctest
    # doctest.testmod()
    # doctest.run_docstring_examples(calc_symmet_log10_levels, globals())
    # doctest.run_docstring_examples(calc_symmet_linear_levels, globals())
