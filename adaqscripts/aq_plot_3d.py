#!/usr/bin/env python
"""
Generalised script to plot 3D AQ data.

From the 'adaqscripts' directory, run via:

.. code-block:: ksh

    ./aq_plot_3d.py [inifilename]

where inifilename includes the full path to the ini file.

Alternatively, from within other python code, run using:

.. code-block:: python

    from adaqscripts import aq_plot_3d
    aq_plot_3d.aq_plotting_3d(inifilename='inifilename')


"""
from __future__ import print_function

import os
import sys
import datetime
import iris

#Find adaq directory -
# should be in directory above the one containing this file
adaq_path = os.path.dirname(os.path.realpath(__file__))+'/../'
sys.path.append(adaq_path)
import adaqscripts.aq_plot
import adaqcode.adaq_vertical_plotting
import adaqcode.field_layer
import adaqcode.field_plot
import adaqcode.plotting_functions

def aq_plotting_3d(inifilename=None, verbose=1):
    """
    Plotting of 3D air quality data.
    Currently only plots cross sections along waypoints and a map of the
    location of these waypoints.

    :param inifilename: String giving filename of ini file.
    :param verbose: Level of print output:

                * 0 = No extra printing
                * 1 = Standard level (recommended)
                * 2 = Extra level for debugging

    :returns: (ini_dict, md_list)

     * **ini_dict** - Dictionary of an :class:`inifile` object.
     * **md_list** - Model Data list - each item should be in
       :class:`adaq_data.ADAQData` format, or a subclass thereof.

    >>> ini_dict, md_list = aq_plotting_3d(verbose=0)
    ... # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    Current version of python: ...
    [...]
    Current version of iris: ...
    Beginning at  ...
    Reading inifile .../adaqscripts/aq_plot_3d.ini
    Getting model data for  aqum_3d  at  ...
    Adding missing times
    Plotting cross sections for O3
    Plotting cross sections for BL_depth
    Not plotting cross section for BL_depth from aqum_3d \
(no vertical coordinate)
    """
    print('Current version of python:', sys.version)
    print('Current version of iris:', iris.__version__)
    print('Beginning at ', datetime.datetime.now())

    #Read in all required data
    ini_dict, sites_data, od, md_list = adaqscripts.aq_plot.aq_get_data(
        inifilename=inifilename,
        keep_all_data=True,
        surface_only=False,
        defaultinifile='adaqscripts/aq_plot_3d.ini')

    #Prepare all the data as required
    ini_dict, od, md_list, od_stats, md_list_stats = adaqscripts.aq_plot.aq_prepare_data(
        ini_dict, od, md_list)

    #Do plotting
    for short_name in ini_dict['short_name_list']:
        #Plot cross sections
        if ini_dict.get('cross_section', False):
            adaqcode.adaq_vertical_plotting.plot_md_cross_section(
                ini_dict, md_list, short_name,
                plot_bl_depth=True, verbose=verbose)

    return ini_dict, md_list

if __name__ == '__main__':

##    import doctest
##    doctest.testmod()

    aq_plotting_3d()
