""" Program to read in weather regime data and add to AQUM data cube. """
from __future__ import division
from __future__ import print_function

from six.moves.builtins import str
from six.moves.builtins import range

import os
from datetime import datetime

import iris
import numpy as np
import matplotlib.pyplot as plt

import cube_time
import timeseries_plot
import plotting_functions

def ffdate2date(ffdate):
    """ Converts from fixed format (%Y%m%d) to date format. """
    return datetime.date(datetime.strptime(ffdate, "%Y%m%d"))

def read_regime_txt(filein):
    """
    Read in dates and regimes from text file into numpy 2D array.

    :param filein: .txt file with columns year, month, day, regime.
    :return: numpy 2D array of dates and weather regimes.

    >>> import config
    >>> sample_data_path = config.SAMPLE_DATADIR+'weather_regimes/'
    >>> wr = sample_data_path + 'daily_30regimes_since2010.txt'
    >>> regime_dates = read_regime_txt(wr) # doctest: +ELLIPSIS
    Getting regime data from  .../daily_30regimes_since2010.txt
    >>> print(regime_dates.shape) # doctest: +ELLIPSIS
    (..., 2)
    >>> print(regime_dates[0,0])
    2010-01-01
    >>> print(regime_dates[0,1])
    19
    """

    print('Getting regime data from ', filein)

    #Read in just first 4 columns of .txt file - year, month, day, regime
    data = np.genfromtxt(filein, dtype=None, usecols=(0, 1, 2, 3),
                         names='year, month, day, regime')

    #Initialise empty 2D array for date and weather regime
    regime_dates = np.empty((len(data), 2), dtype=object)

    #Collate year, month, day to form a date string
    for i, datum in enumerate(data):
        date_string = str((datum['year']))
        if datum['month'] < 10:
            date_string = date_string+'0'
        date_string = date_string+str((datum['month']))
        if datum['day'] < 10:
            date_string = date_string+'0'
        date_string = date_string+str((datum['day']))

        #Convert yyyymmdd date string to date format
        #Populate array with formatted date and regime number
        regime_dates[i, 0] = ffdate2date(str(date_string))
        regime_dates[i, 1] = datum['regime']

    return regime_dates


def add_regime_coord(cube, regime_data, data_dims=1):
    """
    Add weather regime as an auxiliary coordinate to the input cube.

    :param regime_data: numpy 2D array of dates and weather regimes.
    :param data_dims: Must point at the time coordinate of the cube.
                      Typically found at 1 for sites_cube_list, 0 for
                      gridded_cube_list.

    >>> import config
    >>> import adaq_functions
    >>> sample_data_path = config.SAMPLE_DATADIR+'weather_regimes/'
    >>> wr = sample_data_path + 'daily_30regimes_since2010.txt'
    >>> regime_dates = read_regime_txt(wr) # doctest: +ELLIPSIS
    Getting regime data from  .../daily_30regimes_since2010.txt
    >>> ini_data, sites_data, od, md_list = adaq_functions.get_exampledata()
    ... # doctest: +ELLIPSIS
    Reading inifile ...example_data_1days.ini
    Number of sites:  5
    >>> for cube in od.sites_cube_list:
    ...     cube = add_regime_coord(cube, regime_dates)
    >>> for coord in od.sites_cube_list[0].coords():
    ...     print(coord.name())
    site_id
    time
    abbrev
    latitude
    longitude
    site_altitude
    site_name
    site_type
    regime
    >>> for md in md_list:
    ...     for cube in md.sites_cube_list:
    ...         cube = add_regime_coord(cube, regime_dates)
    >>> scl0 = md_list[0].sites_cube_list[0]
    >>> names = [str(coord.name()) for coord in scl0.coords()]
    >>> print('regime' in  names)
    True
    >>> print(len(scl0.coord('time').points) == len(scl0.coord('regime').points))
    True
    """

    #Get cube's datetime points
    dtpoints = cube_time.cube_tpoints_dt(cube)

    regime_list = []
    #Add regime number to list where date in cube matches the regime date
    for point in dtpoints:
        for i in range(len(regime_data)):
            if datetime.date(point) == regime_data[i, 0]:
                regime_list.append(regime_data[i, 1])

    #Add list of regimes to cube as a coordinate
    cube.add_aux_coord(iris.coords.AuxCoord(regime_list, long_name='regime'),
                       data_dims)

    return cube


def plot_regime_bar(cube, plotdir='./'):
    """
    Plot bar chart of the frequency of weather regime occurrence.

    :param cube: Iris cube, which has a 'regime' coordinate and a 'time'
                 coordinate.
    :param plotdir: Directory to save resulting plot in.

    Get regime data:

    >>> import config
    >>> wr = config.SAMPLE_DATADIR + \
    'weather_regimes/daily_30regimes_since2010.txt'
    >>> regime_dates = read_regime_txt(wr) # doctest: +ELLIPSIS
    Getting regime data from  .../daily_30regimes_since2010.txt

    Load an example sites cube:

    >>> data_path = config.SAMPLE_DATADIR+'sites_cube_list/'
    >>> cube = iris.load_cube(data_path + 'aurn_5days.nc',
    ... iris.AttributeConstraint(short_name='PM2p5'))

    Give it a regime coordinate:

    >>> cube = add_regime_coord(cube, regime_dates)

    Plot this as a time-series for the first site

    >>> plotdir = config.CODE_DIR + "/adaqdocs/figures"
    >>> tsp = plot_regime_bar(cube[0], plotdir=plotdir) # doctest: +ELLIPSIS
    Plotting regime bar chart
    Saved figure  .../Regime_Bar.png

    .. image:: ../adaqdocs/figures/Regime_Bar.png
       :scale: 50%

    """

    print("Plotting regime bar chart")

    #List of regime numbers for plotting on x axis
    x = list(range(1, 31))
    #Initialise empty list for plotting on y axis
    y = []

    #Get cube's datetime points
    dtpoints = cube_time.cube_tpoints_dt(cube)

    #Record number of days of occurrence for each regime in y axis list
    for iregime in x:
        try:
            filter_cube = cube.extract(iris.Constraint(regime=iregime))
            nhours = len(filter_cube.coord('time').points)
            y.append(nhours//24)
        except:
            #If regime does not occur in this period, mark as 0
            y.append(int(0))

    #Plot
    plt.bar(x, y, align='center', width=0.75, color='teal', alpha=0.9)
    plt.xticks(list(range(1, 31)), x, ha='center', fontsize='x-small')
    plt.xlabel('Weather Regime')
    plt.ylabel('Number of Days')
    plt.suptitle('Frequency of Weather Regime Occurrence', fontsize=16)
    plt.title(str(dtpoints[0].strftime("%d/%m/%Y %H:%M"))+' to '+
              str(dtpoints[-1].strftime("%d/%m/%Y %H:%M")), fontsize=12)

    #Save
    if plotdir[-1] != '/':
        plotdir += '/'
    if not os.path.isdir(plotdir):
        print('Creating output directory:', plotdir)
        os.makedirs(plotdir)

    filename = 'Regime_Bar.png'
    plt.savefig(plotdir+filename)
    print('Saved figure ', plotdir+filename)
    plt.close()


def plot_regime_timeseries(cube, plotdir='./'):
    """
    Produce time series plot illustrating the daily weather regime
    classification.

    :param cube: Iris cube, which has a 'regime' coordinate and a 'time'
                 coordinate.
    :param plotdir: Directory to save resulting plot in.

    :return: :class:`timeseries_plot.TimeSeriesPlot` object

    >>> import config
    >>> wr = config.SAMPLE_DATADIR + \
    'weather_regimes/daily_30regimes_since2010.txt'
    >>> regime_dates = read_regime_txt(wr) # doctest: +ELLIPSIS
    Getting regime data from  .../daily_30regimes_since2010.txt

    Load an example sites cube:

    >>> data_path = config.SAMPLE_DATADIR+'sites_cube_list/'
    >>> cube = iris.load_cube(data_path + 'aurn_5days.nc',
    ... iris.AttributeConstraint(short_name='PM2p5'))

    Give it a regime coordinate:

    >>> cube = add_regime_coord(cube, regime_dates, 1)

    Plot this as a time-series for the first site

    >>> plotdir = config.CODE_DIR + "/adaqdocs/figures"
    >>> tsp = plot_regime_timeseries(cube[0], plotdir=plotdir)
    ... # doctest: +ELLIPSIS
    Plotting regime time-series
    Saved figure  .../Regime_Timeseries.png

    .. image:: ../adaqdocs/figures/Regime_Timeseries.png
       :scale: 50%

    """

    print('Plotting regime time-series')

    #Set up time-series-plot
    tsp = timeseries_plot.TimeSeriesPlot()

    tsp.add_line(cube, y=cube.coord('regime'), colour='teal')

    #Create a plot and save
    tsp.title = 'Daily Weather Regime Classification'
    tsp.ylabel = 'Weather Regime'
    tsp.legend = False
    tsp.plot()
    ax = plt.gca()
    ax.set_yticks(list(range(2, 31, 2)))
    tsp.save_fig(plotdir=plotdir, filename='Regime_Timeseries.png')

    return tsp


def plot_regime_stats(ini_dict,
                      xstat_list=None,
                      ystat_list=None,
                      statsfile=None):
    """
    Produce scatter plots of each xstat vs each ystat, labelling points
    by corresponding weather regime.
    By default, picks out meanobs vs a selected group of useful stats.

    :param ini_dict: Dictionary of a :class:`inifile` object. Should contain:

                      * 'plot_dir' - location of statistics file to read in,
                        plus output location
                      * 'short_name_list'

    :param xstat_list: List of stats to plot on the x axis, in turn.
                       If set to None, defaults to ['meanobs']
    :param ystat_list: List of stats to plot on the y axis, in turn.
                       If set to None, defaults to
                       ['mnmb', 'fge', 'bias', 'meanmod']

    *Available stats for plotting:*
    'nsites': 'Number of sites',
    'n': 'Number of points',
    'correlation': 'Correlation',
    'bias': 'Bias',
    'nmb': 'Normalised Mean Bias',
    'mnmb': 'Modified Normalised Mean Bias',
    'nmge': 'Normalised Mean Gross Error',
    'fge': 'Fractional Gross Error',
    'rmse': 'Root Mean Square Error',
    'fac2': 'Factor of 2',
    'ioa': 'Index of Agreement',
    'threshold': 'Threshold',
    'ORSS': 'Odds Ratio Skill Score',
    'hitrate': 'Hitrate',
    'falsealarmrate': 'False Alarm Rate',
    'falsealarmratio': 'False Alarm Ratio',
    'o>=t_m>=t': 'Number Obs >= Threshold and Model >= Threshold',
    'o<t_m>=t': 'Number Obs < Threshold and Model >= Threshold',
    'o>=t_m<t': 'Number Obs >= Threshold and Model < Threshold',
    'o<t_m<t': 'Number Obs < Threshold and Model < Threshold',
    'maxobs': 'Maximum Observation Value',
    'maxmod': 'Maximum Model Value',
    'meanobs': 'Mean Observation Value',
    'meanmod': 'Mean Model Value',
    'sdobs': 'Standard Deviation of Observations',
    'sdmod': 'Standard Deviation of Model',
    'perc_correct': 'Percentage of Correct values',
    'perc_over': 'Percentage of Over-predicions',
    'perc_under': 'Percentage of Under-predictions',
    'units': 'Units'

    Firstly get some example data, and ensure that plot_dir is set.
    Also shorten the short_names required to be plotted.

    >>> import adaq_functions
    >>> ini_dict, sites_data, od, md_list = adaq_functions.get_exampledata(
    ... exampletype="full") # doctest: +ELLIPSIS
    Reading inifile .../example_data_5days.ini
    Number of sites:  5
    >>> import config
    >>> ini_dict['plot_dir'] = config.CODE_DIR + "/adaqdocs/figures/"
    >>> ini_dict['short_name_list'] = ['O3']
    >>> ini_dict['calc_stats_format_list'] = ['csv']

    Before any plotting can be done, first need to calculate some statistics:

    >>> statsfile = 'regime_stats'
    >>> tsstats = adaq_functions.calc_stats(
    ... ini_dict, od, md_list, statsfile_prefix=statsfile) # doctest: +ELLIPSIS
    Statistics saved to  .../regime_stats.csv

    Can now read in the file which contains statistics and plot these.
    >>> plot_regime_stats(ini_dict, ystat_list=['bias'],
    ... statsfile=statsfile) # doctest: +ELLIPSIS
    Statistics found at .../regime_stats.csv
    Saved figure  .../meanobs_vs_bias_O3.png
    """

    #Set up defaults
    if xstat_list is None:
        xstat_list = ['meanobs']
    if ystat_list is None:
        ystat_list = ['mnmb', 'fge', 'bias', 'meanmod']

    plotdir = ini_dict['plot_dir']
    if plotdir[-1] != '/':
        plotdir += '/'
    if statsfile is None:
        filein = plotdir+'stats.csv'
    else:
        filein = plotdir + statsfile+'.csv'
    try:
        #Read in stats.csv file from plot directory
        stats = np.genfromtxt(filein, dtype=str, skip_header=1,
                              autostrip=True, delimiter=',')
        print('Statistics found at ', filein)
    except:
        #If stats file not found in directory, raise error
        raise IOError("No statistics file found in plotdir for reading")

    #Initialise certain stats to carry units
    needs_units = ['meanobs', 'meanmod', 'maxobs', 'maxmod', 'bias', 'rmse']

    regimes = []
    #Pick out regime numbers from model label
    #Ignore first two columns, 'Phenomenon' and 'Statistic'
    for iregime in stats[0][2:]:
        if iregime[-2] in ['1', '2', '3']:
            regimes.append(iregime[-2:])
        else:
            regimes.append(iregime[-1])

    #Produce species-specific plots
    for short_name in ini_dict['short_name_list']:
        for xstat in xstat_list:
            for ystat in ystat_list:
                for item in stats:
                    #Create list of selected stat
                    if item[0] == short_name and item[1].split()[0] == xstat:
                        x = item[2:]
                    elif item[0] == short_name and item[1].split()[0] == ystat:
                        y = item[2:]
                    elif item[0] == short_name and item[1] == 'units':
                        units = item[2]
                        units = plotting_functions.units_str(units)

                plt.figure(figsize=(8, 7.25))
                plt.scatter(x, y, s=0)

                #Add units to axes labels if appropriate
                if xstat in needs_units:
                    xunits = ' ('+units+')'
                else:
                    xunits = ''
                if ystat in needs_units:
                    yunits = ' ('+units+')'
                else:
                    yunits = ''
                plt.xlabel(xstat+xunits)
                plt.ylabel(ystat+yunits)
                plt.title(short_name)

                #Use regime number instead of dot to indicate data point
                for i, regime in enumerate(regimes):
                    plt.annotate(regime, (x[i], y[i]),
                                 horizontalalignment='center',
                                 verticalalignment='center',
                                 color='teal')

                ymin, ymax = plt.ylim()
                xmin, xmax = plt.xlim()
                #Add axis line at y=0 if there's +ve and -ve data
                if ymin < 0 < ymax:
                    plt.axhline(color='gray')

                #Equalise axes when plotting similar stats
                if xstat == 'meanobs' and ystat == 'meanmod':
                    plt.ylim(min(xmin, ymin), max(xmax, ymax))
                    plt.xlim(min(xmin, ymin), max(xmax, ymax))

                filename = xstat+'_vs_'+ystat+'_'+short_name+'.png'
                plt.savefig(plotdir+filename)

                print('Saved figure ', plotdir+filename)
                plt.close()



if __name__ == '__main__':

    import doctest
    doctest.testmod()
