#!/usr/bin/env python
'''
Generalised script to generate trajectory plots from NAME trajectory
output

Running the code
~~~~~~~~~~~~~~~~
To run the trajectory plotting code you need to create an inifile containing
information about the location of your name output where you would like
your plots to be saved. The plots are then created by typing

.. code-block:: ksh

   ./plot_trajectory.py [inifilename]

where inifilename includes the full path to the ini file.

Alternatively, from within other python code, run using:

.. code-block:: python

    from adaqscripts import plot_trajectory
    plot_trajectory.trajectory_plotting(inifilename='inifilename')

or if passing an ini_dict rather than an ini file

.. code-block:: python

    from adaqscripts import plot_trajectory
    # Set up ini_dict, e.g.
    ini_dict = {}
    ini_dict['key'] = value
    plot_trajectory.trajectory_plotting(ini_dict=ini_dict)

Contents of the inifile
~~~~~~~~~~~~~~~~~~~~~~~
The inifile needs to contain the following items:

+---------------+--------------------------------------------------------------+
|models_dir_list|List of the full directory paths of the data to be read in.   |
|               |Linux wildcards such as * and ? may be used in the path names.|
+---------------+--------------------------------------------------------------+
|plot_dir       |Full directory path of the location you would like your plots |
|               |to be saved to.                                               |
+---------------+--------------------------------------------------------------+

The inifile may also contain any of the following items:

+--------------------+---------------------------------------------------------+
|extent_list         |Extent of the plot in longitude and latitude. Should be  |
|                    |a 4-digit list in the order: xmin, xmax, ymin, ymax.     |
+--------------------+---------------------------------------------------------+
|title               |A string which will be used as the plot title.           |
+--------------------+---------------------------------------------------------+
|marker_interval     |The distance between the markers along a trajectory in   |
|                    |hours. The default is 6-hours.                           |
+--------------------+---------------------------------------------------------+
|plotname            |Name of the file to save the plot to. Needs to contain a |
|                    |suitable extension (e.g. png, pdf). The default is       |
|                    |TrajectoryPlot.png                                       |
+--------------------+---------------------------------------------------------+
|mapping             |Mapping background to add to plume plot. Options are     |
|                    |'states', 'countries' (default), 'coastline' or 'None'.  |
|                    |Note that when 'countries' is selected coastlines and the|
|                    |great lakes are plotted and when 'states' are selected   |
|                    |counties are included for zoomed in maps.                |
+--------------------+---------------------------------------------------------+
|release_info_list   |List containing the release longitude, the release       |
|                    |latitude and the release time. Any format can be used but|
|                    |there must be three elements in the list.                |
|                    |e.g. ["17.33W", "64.42N", "08:00, 08/02/2017"]           |
+--------------------+---------------------------------------------------------+
|mo_branding         |To add Met Office branding. True or False. Default is    |
|                    |False                                                    |
+--------------------+---------------------------------------------------------+
'''
from __future__ import print_function

import datetime
import os
import sys
import iris

#Find adaq directory -
# should be in directory above the one containing this file
adaq_path = os.path.dirname(os.path.realpath(__file__))+'/../'
sys.path.append(adaq_path)
import adaqcode.adaq_functions
import adaqcode.inifile
import adaqcode.adaq_plotting

def trajectory_plotting(inifilename=None, ini_dict=None, verbose=1):
    '''
    Top-level plotting routine to generate trajectory plots

    :param inifilename: String giving filename of ini file.
    :param inidict: Dictionary of a :class:`inifile` object. If this is given,
                     then it is used instead of reading in from inifilename.
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
     * Reads in model data
     * Plots data using style given in inifile

    >>> ini_dict, md_list = trajectory_plotting()
    ... # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    Current version of python: ...
    [...]
    Current version of iris: ...
    Beginning at  ...
    Reading inifile .../adaqscripts/plot_trajectory.ini
    Getting model data for  Data_Traj_C1_*.txt  at  ...
    <class '...TrajData'> - \
    Subclass of ADAQData - Contains:
    sites_cube_list:
    < No cubes >
    gridded_cube_list:
    < No cubes >
    trajectory_cube_list:
    0: U Turb / (m/s)                      (time: 193)
    1: U Turb / (m/s)                      (time: 193)
    Saved figure  .../TrajectoryPlot.png
    Finished at  ...

    '''

    print('Current version of python:', sys.version)
    print('Current version of iris:', iris.__version__)
    print('Beginning at ', datetime.datetime.now())

   # Read ini file
    if ini_dict is None:
        ini_dict = adaqcode.inifile.get_inidict(
            inifilename=inifilename,
            defaultfilename='adaqscripts/plot_trajectory.ini')

    #Ensure models_list is setup if not already given in inifile
    #Note - requires models_dir_list to be given instead.
    if 'models_list' not in ini_dict:
        ini_dict['models_list'] = [
            os.path.basename(path) if path is not None else None
            for path in ini_dict['models_dir_list']]

    # Read in model data
    ini_dict['models_fmt_list'] = ['traj']
    md_list = adaqcode.adaq_functions.get_models(ini_dict)
    for md in md_list:
        print(md)

        adaqcode.adaq_plotting.plot_trajectories(
            ini_dict, md, verbose=verbose)


    print('Finished at ', datetime.datetime.now())

    return ini_dict, md_list

if __name__ == '__main__':

    trajectory_plotting()

    #import doctest
    #doctest.testmod()
    #doctest.run_docstring_examples(get_ini, globals())
