#!/usr/bin/env python
"""
Generalised script to load netcdf statistics file and plot

Running rolling_stats.py
------------------------

To run plot_rolling_stats() which will load the statistics file
and plot requested items such as csv statistics, timeseries of statistics:

.. code-block:: ksh

  ./rolling_stats.py [inifilename]

where inifilename is an inifile similar to rolling_stats.ini. If not given,
then defaults to rolling_stats.ini in the same directory as the script.

To run aq_create_rolling_stats() which allows netcdf statistics file to be
created by looping over aq_plot code:
Instead run using the -c ("create") argument :

.. code-block:: ksh

  ./rolling_stats.py -c [inifilename]

where inifilename is an inifile similar to aq_plot.ini. If not given,
then defaults to aq_plot.ini in the same directory as the script.

To run aq_create_rolling_stats(), appending to a previously-generated
netcdf file, then also include the --oldnc option:

.. code-block:: ksh

  ./rolling_stats.py -c --oldnc [inifilename]

To display this help on the command line:

.. code-block:: ksh

  ./rolling_stats.py -h

"""
from __future__ import print_function

from six.moves.builtins import range
import datetime
import argparse
import warnings
import os
import sys

import iris
import cf_units

#Find adaq directory -
# should be in directory above the one containing this file
adaq_path = os.path.dirname(os.path.realpath(__file__))+'/../'
sys.path.append(adaq_path)
import adaqcode
from adaqscripts import aq_plot

def aq_create_rolling_stats(inifilename=None, ini_dict=None, sites_data=None,
                            new_ncfile=True):
    """
    Routine to loop over aq_plot.py code and read data/produce netcdf statistic
    file one day at a time. This statistic file will be appended to after being
    created on the first day.
    The created netcdf file contains a 2D array whose values corresponding to
    each statistic for each time (each day)

    :param inifilename: String giving filename of ini file.
    :param ini_dict: Dictionary of a :class:`inifile` object. If this is given,
                     then it is used instead of reading in from inifilename.
    :param sites_data: numpy ndarray containing site information data
                       from a :class:`sites_info.SitesInfo` object.
                       If this is given, then it is used instead of reading in
                       from a file/cube.
    :param new_ncfile: Start a new netcdf file to save to on day1, overwriting
                       any existing file of the same name.

    For example purposes, increase the default range in the aq_plot.ini file,
    and disable aq_plot's own saving:

    >>> ini_dict = adaqcode.inifile.get_inidict(
    ... defaultfilename='adaqscripts/aq_plot.ini') # doctest: +ELLIPSIS
    Reading inifile .../aq_plot.ini
    >>> ini_dict['range_days'] = 3
    >>> ini_dict['save_cubes_nc'] = False

    Now run the code:

    >>> aq_create_rolling_stats(ini_dict=ini_dict) # doctest: +ELLIPSIS
    Date: 2014-03-26
    Number of sites:  5
    Creating observation data at  ...
    Reading obs data files
    Found obs for  .../ABD_20140101_20140818.txt
    ...
    Found obs for  .../YW_20140101_20140818.txt
    Creating obs cubes
    Getting model data for  oper  at  ...
    Getting model data for  name_casestudy  ...
    Removing gridded_cube_list
    Adding missing times
    Converting to DAQI
    Saved to  .../stats_Obs-AURN_Mod-name_casestudy.nc
    Saved to  .../stats_Obs-AURN_Mod-oper.nc
    Statistics saved to  .../stats.csv
    Statistics saved to  .../stats.wiki
    Date: 2014-03-27
    Creating observation data at ...
    Reading obs data files
    Found obs for  .../ABD_20140101_20140818.txt
    ...
    Creating obs cubes
    Getting model data for  oper  at  ...
    Getting model data for  name_casestudy  at ...
    Removing gridded_cube_list
    Adding missing times
    Converting to DAQI
    Saved to  .../stats_Obs-AURN_Mod-name_casestudy.nc
    Saved to  .../stats_Obs-AURN_Mod-oper.nc
    Statistics saved to  .../stats.csv
    Statistics saved to  .../stats.wiki
    Date: 2014-03-28
    Creating observation data at  ...
    Reading obs data files
    Found obs for .../ABD_20140101_20140818.txt
    ...
    Saved to  .../stats_Obs-AURN_Mod-name_casestudy.nc
    Saved to  .../stats_Obs-AURN_Mod-oper.nc
    Statistics saved to  .../stats.csv
    Statistics saved to  .../stats.wiki

    """

    if ini_dict is None:
        ini_dict = adaqcode.inifile.get_inidict(
            inifilename=inifilename,
            defaultfilename='adaqscripts/aq_plot.ini')
    elif inifilename is not None:
        warnings.warn("Inifile not read in - using data from input ini_dict "
                      "instead")
    #Reserve short_name_list for future use (this gets added to later)
    short_name_list = ini_dict['short_name_list'][:]

    #Ensure netcdf output is enabled
    if 'calc_stats_format_list' not in ini_dict:
        ini_dict['calc_stats_format_list'] = []
    if 'nc' not in ini_dict['calc_stats_format_list']:
        ini_dict['calc_stats_format_list'].append('nc')

    #Don't need to produce contours
    #Turn off to ensure that gridded data is not kept in memory
    ini_dict['contours'] = False

    #Calculate statistics for each day and save to combined netcdf file
    #in ini_dict['plot_dir']
    dates = [ini_dict['start_datetime'].date() + datetime.timedelta(days=n)
             for n in range(ini_dict['range_days'])]
    firstdate = True
    for date in dates:
        print('Date:', date)

        #Overwrite value of start_datetime, end_datetime and range_days
        #in ini_dict
        ini_dict['start_datetime'] = datetime.datetime(
            date.year, date.month, date.day, 0)
        ini_dict['end_datetime'] = ini_dict['start_datetime'] + \
                                   datetime.timedelta(days=1)
        ini_dict['range_days'] = 1
        #If need to calculate DAQI, then need to also include previous day
        if ini_dict.get('daqi', False):
            ini_dict['start_datetime'] = ini_dict['start_datetime'] - \
                                         datetime.timedelta(days=1)
            ini_dict['range_days'] = 2

        #Reset values in short_name_list (nb use copy [:], not pointer)
        ini_dict['short_name_list'] = short_name_list[:]

        #Now read in all required data
        ini_dict, sites_data, od, md_list = aq_plot.aq_get_data(
            ini_dict=ini_dict,
            sites_data=sites_data,
            keep_all_data=False)

        #Prepare all the data as required
        ini_dict, od, md_list, od_stats, md_list_stats = \
                  aq_plot.aq_prepare_data(ini_dict, od, md_list)


        #If including DAQI, then extract required date only
        #instead of saving both days
        if ini_dict.get('daqi', False):
            #Set up new cubelist to overwrite with
            new_od_scl = iris.cube.CubeList()
            for cube in od_stats.sites_cube_list:
                new_od_cube = cube.extract(iris.Constraint(
                    time=lambda t: t.point.date() == date))
                if new_od_cube is not None:
                    #Only keep this cube in cube list if date have been
                    #matched, therefore not None
                    #Nb may have only got data for previous day
                    new_od_scl.append(new_od_cube)
            od_stats.sites_cube_list = new_od_scl
            for md in md_list_stats:
            #Set up new cubelist to overwrite with
                new_md_scl = iris.cube.CubeList()
                for cube in md.sites_cube_list:
                    new_md_cube = cube.extract(iris.Constraint(
                        time=lambda t: t.point.date() == date))
                    if new_md_cube is not None:
                        #Only keep this cube in cube list if date have been
                        #matched, therefore not None
                        #Nb may have only got data for previous day
                        new_md_scl.append(new_md_cube)
                md.sites_cube_list = new_md_scl

        #Check that we have some data, otherwise no point in continuing
        if not od.sites_cube_list:
            #Nb print rather than warning, as warning would only be
            #issued for first missing date
            print('No observation data for this date')
            continue

        model_data = False
        for md in md_list_stats:
            if len(md.sites_cube_list) >= 1:
                #Some cubes found
                model_data = True
        if not model_data:
            #Nb print rather than warning, as warning would only be
            #issued for first missing date
            print('No model data for this date')
            continue


        #Fix time-coord to have point in the middle of expected date
        #(otherwise if some times at the beginning/end etc missing, then
        # time chosen will be not as expected, and then later when missing
        # times are added, this adds times in funny unexpected places
        # and can take a very long time/add in a lot of unneeded nan data!)
        tunit = cf_units.Unit('hours since epoch', calendar='gregorian')
        tpt = datetime.datetime(date.year, date.month, date.day, 12)
        tbounds = [ini_dict['start_datetime'], ini_dict['end_datetime']]
        tcoord = iris.coords.DimCoord(
            [tunit.date2num(tpt)],
            standard_name='time',
            bounds=[tunit.date2num(b) for b in tbounds],
            units=tunit)

        #Calculate statistics and save to file
        nc_append = True
        if firstdate and new_ncfile:
            nc_append = False
        adaqcode.adaq_functions.calc_stats(
            ini_dict, od_stats, md_list_stats,
            thresholds=aq_plot.STATS_THRESHOLDS,
            nc_append=nc_append,
            stats_tcoord=tcoord)
        firstdate = False



def calc_stats(ini_dict, stats_cubelist):
    """
    Routine to calculate statistics from stats_cubelist.
    Similar to :func:`adaq_functions.calc_stats`, but takes
    list of statistics cubes and converts them back to
    :class:`TimeSeriesStats` objects, rather than creating
    :class:`TimeSeriesStats` objects from raw obs and model cubes
    derived from od and md_list.

    .. note:: Tested as part of plot_rolling_stats()
    """


    if ini_dict.get('calc_stats', False):
        tsstats_list = []
        startdate, enddate = None, None
        for cube in stats_cubelist:
            if cube.attributes['short_name'] in \
               ini_dict.get('short_name_list', []):
                stats = adaqcode.timeseries_stats.TimeSeriesStats([], [])
                #Calculate a nan-mean over all times
                stats.statscube = cube.collapsed(
                    'time',
                    adaqcode.cube_statistics.CUBE_AGGREGATORS['NANMEAN'])
                #Convert all the information in the cube back into the
                #attributes of the TimeSeriesStats object
                stats.convert_cube_to_attr()
                tsstats_list.append(stats)
                #Extract times ready for header
                startdt, enddt = adaqcode.cube_time.get_startenddt(
                    cube=stats.statscube)
                if startdate is None:
                    startdate = startdt
                if startdt < startdate:
                    startdt = startdate
                if enddate is None:
                    enddate = enddt
                if enddt > enddate:
                    enddate = enddt

        percentile_value = ini_dict.get('percentile', 95)

        header = ('Statistics for '+
                  startdate.strftime("%a %d-%b-%Y %H:%M")+
                  ' - '+enddate.strftime("%a %d-%b-%Y %H:%M"))
        adaqcode.timeseries_stats.save_stats(
            tsstats_list, percentile_value,
            directory=ini_dict['plot_dir'],
            header=header,
            format_list=ini_dict.get('calc_stats_format_list', []))


def plot_rolling_stats(inifilename=None, ini_dict=None):
    """
    Top level routine to calculate and plot rolling statistics from
    pre-saved netcdf files.
    Statistics plotted depend on settings in input inifile.

    :param inifilename: String giving filename of ini file.
                        If set to None: if ini_dict is also not set, then
                        defaults to rolling_stats.ini instead.
                        However if ini_dict is given this filename is ignored
                        (even if set) and ini_dict used instead.

    :param ini_dict: Dictionary of a :class:`inifile` object. If this is given,
                     then it is used instead of reading in from inifilename.
                     Containing:

                     * short_name_list
                     * rolling_stats_file_list (list of netcdf statistics files)
                     * plot_dir
                     * Other optional settings to determine which types of
                       plots to produce, eg 'monthly_stats' to produce monthly
                       versions of the main outputs.

    :returns: (ini_dict, stats_cubelist), where stats_cubelist is
              the iris cubelist containing all the statistics cubes read in from
              netcdf files.

    >>> ini_dict, stats_cubelist = plot_rolling_stats() # doctest: +ELLIPSIS
    Reading inifile .../rolling_stats.ini
    Saved figure  .../Timeseries_of_bias_O3.png
    Saved figure  .../Timeseries_of_bias_O3_monthly.png
    Saved figure  .../Timeseries_of_rmse_O3.png
    Saved figure  .../Timeseries_of_rmse_O3_monthly.png
    Saved figure  .../Timeseries_of_bias_PM2p5.png
    Saved figure  .../Timeseries_of_bias_PM2p5_monthly.png
    Saved figure  .../Timeseries_of_rmse_PM2p5.png
    Saved figure  .../Timeseries_of_rmse_PM2p5_monthly.png
    Statistics saved to  .../stats.csv
    Statistics saved to  .../stats_monthly_aqum_oper.csv
    """

    #Read ini file
    if ini_dict is None:
        ini_dict = adaqcode.inifile.get_inidict(
            inifilename=inifilename,
            defaultfilename='adaqscripts/rolling_stats.ini')
    elif inifilename is not None:
        warnings.warn("Inifile not read in - using data from input ini_dict "
                      "instead")

    #Set up some extra defaults if not set in ini_dict
    ini_dict['plot_dir'] = ini_dict.get('plot_dir', './')
    ini_dict['monthly_stats'] = ini_dict.get('monthly_stats', False)

    #Load in the rolling stats
    stats_cubelist = iris.cube.CubeList()
    for filename in ini_dict['rolling_stats_file_list']:
        stats_cubes = iris.load(filename)
        for cube in stats_cubes:
            stats_cubelist.append(cube)

    #If DAQI needs to be plotted, then append these to short-names
    #Modify short name list to add _DAQI where appropriate and
    #also add DAQI in
    if ini_dict.get('daqi', False):
        short_name_list_with_daqi = []
        for short_name in ini_dict.get('short_name_list', []):
            if short_name in adaqcode.aq_indices.DAQI_LEVELS:
                short_name_list_with_daqi.append(short_name + '_DAQI')
        short_name_list_with_daqi.append('DAQI')
        ini_dict['short_name_list'] += short_name_list_with_daqi

    #Plot timeseries of statistics
    if ini_dict.get('timeseries_of_stats', False):

        if ini_dict['monthly_stats']:
            #Set up dictionary, containing lists of tsstats objects,
            #with keys of model labels
            tsstats_dict = {}

        for short_name in ini_dict.get('short_name_list', []):
            #Extract just the cubelist corresponding to this short_name (sn)
            sn_cubes = stats_cubelist.extract(iris.AttributeConstraint(
                short_name=short_name))

            if ini_dict['monthly_stats']:
                #Convert to monthly, keeping year information
                sn_cubes_monthly = iris.cube.CubeList()

                for cube in sn_cubes:
                    sn_cube_monthly = adaqcode.cube_statistics.aggregate_time(
                        cube, period='yearmonth')

                    sn_cubes_monthly.append(sn_cube_monthly)

            for statistic in ini_dict.get('timeseries_of_stats_list', []):
                #Plot this statistic
                adaqcode.timeseries_plot.tsp_statistic(
                    sn_cubes, statistic, ini_dict.get('plot_dir', './'))

                if ini_dict['monthly_stats']:
                    #Produce monthly plots
                    adaqcode.timeseries_plot.tsp_statistic(
                        sn_cubes_monthly,
                        statistic,
                        ini_dict.get('plot_dir', './'),
                        colours_list=ini_dict.get('line_colours_list', None),
                        filesuffix='_monthly',
                        xcoordname='yearmonthdt')

            if ini_dict['monthly_stats'] and ini_dict.get('calc_stats', False):
                #Need to output to text file
                #For each time, convert to TimeSeriesStats objects
                #treat each time as a different model run
                #Output later using save_stats
                for cube in sn_cubes_monthly:
                    for tcube in cube.slices_over('yearmonth'):
                        #Initalise TimeSeriesStats object with empty lists for
                        #obs and mod as these won't be required
                        stats = adaqcode.timeseries_stats.TimeSeriesStats(
                            [], [],
                            label=' ' + tcube.coord('yearmonthname').points[0])
                        #Set up statistics cube
                        stats.statscube = tcube
                        #Convert to statistics dictionary
                        stats.convert_cube_to_attr()
                        #Add to tsstats_list within tsstats dictionary
                        if tcube.attributes['label'] not in tsstats_dict:
                            tsstats_dict[tcube.attributes['label']] = []
                        tsstats_dict[tcube.attributes['label']].append(stats)


    #Produce statistics file
    if ini_dict.get('calc_stats', False):
        calc_stats(ini_dict, stats_cubelist)

        percentile_value = ini_dict.get('percentile', 95)

        if ini_dict['monthly_stats']:
            #Output monthly statistics to text file
            format_list = ini_dict.get('calc_stats_format_list', ['csv'])
            #Don't allow nc output (doesn't make sense to do this)
            if 'nc' in format_list:
                format_list.remove('nc')

            for label, tsstats_list in tsstats_dict.items():

                adaqcode.timeseries_stats.save_stats(
                    tsstats_list, percentile_value,
                    filename_prefix='stats_monthly_'+label.replace(' ', '_'),
                    directory=ini_dict['plot_dir'],
                    format_list=format_list,
                    header='Monthly Statistics for '+label)



    return ini_dict, stats_cubelist


def parse_commandline_inputs():
    """
    Function to parse the inputs if run on the command line

    As an example add to sys.argv to represent the following command line:

    .. code-block:: ksh

      rolling_stats.py -c myfile.ini

    >>> sys.argv[1:] = ['-c', 'myfile.ini']
    >>> args = parse_commandline_inputs()
    >>> print(args)
    Namespace(create=True, inifilename='myfile.ini', oldnc=False)

    If you wish to use a previously-generated netcdf file to continue
    appending to when creating rolling stats file, then use:

    .. code-block:: ksh

      rolling_stats.py -c --oldnc myfile.ini

    >>> sys.argv[1:] = ['-c', '--oldnc', 'myfile.ini']
    >>> args = parse_commandline_inputs()
    >>> print(args)
    Namespace(create=True, inifilename='myfile.ini', oldnc=True)

    Return sys.argv to defaults so it doesn't confuse other examples:

    >>> sys.argv[1:] = []

    """

    #Set up parser. Use docstring for module as main help string.
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    #Add argument options.
    #Inifilename is optional - by default if not set, then set to None.
    parser.add_argument("inifilename", nargs='?', default=None)
    #Extra argument to convert to calling aq_create_rolling_stats instead.
    parser.add_argument("--create", "-c", help="run aq_create_rolling_stats() "
                        "instead of plot_rolling_stats()",
                        action="store_true")
    parser.add_argument("--oldnc", help="Use old pre-existing netcdf file when "
                        "running aq_create_rolling_stats() rather than starting"
                        "from new file", action="store_true")

    arguments = parser.parse_args()

    return arguments



if __name__ == '__main__':

    #import doctest
    #doctest.testmod()

    #By default run code using command line inputs
    args = parse_commandline_inputs()
    if args.create:
        aq_create_rolling_stats(args.inifilename, new_ncfile=not args.oldnc)
    else:
        plot_rolling_stats(args.inifilename)
