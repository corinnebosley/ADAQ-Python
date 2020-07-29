"""
Code for producing time-series plots.
Contains class for producing the plot, as well as
generic functions and colours variable which
may be useful for other plotting.
"""
import os
import datetime
import warnings
import pytz

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mpldates
import iris
import iris.plot as iplt
import numpy as np

import line_plot
import plotting_functions
import timeseries_stats

if not os.getenv('DISPLAY'):
    #Enable DISPLAY if running under cron
    mpl.use('Agg')


def set_xaxis_date_fmt(ax):
    """
    Set the dates on an x-axis nicely, spread over 2 lines
    If x-axis does not include dates, then has no effect.

    :param ax: matplotlib axis instance

    """
    # work around for bug in matplotlib 1.5.3
    # see: https://github.com/matplotlib/matplotlib/issues/7630
    mpldates.UTC = pytz.UTC

    xax = ax.get_xaxis() # get the x-axis

    #Get hold of existing locator and formatter objects
    major_locator = xax.get_major_locator()
    major_formatter = xax.get_major_formatter()

    if not isinstance(major_formatter, mpldates.AutoDateFormatter):
        #Not using dates, so just return as can't set up dates formats
        return ax

    #Use this to figure out the date range over which the xaxis has
    xaxis_range = major_locator.viewlim_to_dt()[1] - \
                  major_locator.viewlim_to_dt()[0]

    #Now use this range to set up the date locations and formatters sensibly.
    #Minor locator - location of upper datestrings
    #Minor formatter - format of upper datestrings
    #Major locator - location of lower datestrings (larger intervals than minor)
    #Major formatter - format of lower datestrings
    #Major - used for setting gridlines, so set to smaller intervals than minor
    #AutoDateLocator - automatically finds appropriate intervals, generally
    # requesting a minimum of 3 tick locations and a maximum number to ensure
    # text fits.
    # interval_multiples is set to ensure these are put at nice locations
    #MonthLocator etc finds every month (not just nicely intervaled ones)
    if xaxis_range <= datetime.timedelta(hours=3):
        #<=3 hour range - display minutes, eg 03:30UTC \n 02/03/2014
        major_locator = mpldates.AutoDateLocator(minticks=3, maxticks=6,
                                                 interval_multiples=True)
        major_formatter = mpldates.DateFormatter('%H:%M%Z')
        minor_locator = mpldates.WeekdayLocator()
        minor_formatter = mpldates.DateFormatter('\n%d/%m/%Y')
    elif xaxis_range <= datetime.timedelta(days=1):
        #<=1day range - display hours eg 03UTC \n 02/03/2014
        major_locator = mpldates.AutoDateLocator(minticks=3, maxticks=6,
                                                 interval_multiples=True)
        #Don't allow 4-hour interval
        major_locator.intervald[mpldates.HOURLY] = [1, 2, 3, 6, 12]
        major_formatter = mpldates.DateFormatter('%H%Z')
        minor_locator = mpldates.WeekdayLocator()
        minor_formatter = mpldates.DateFormatter('\n%d/%m/%Y')
    elif xaxis_range <= datetime.timedelta(days=3):
        #<=3 day range - display every 12 hours and dates eg 12UTC \n 02/03/2014
        major_locator = mpldates.HourLocator(byhour=[0, 12])
        major_formatter = mpldates.DateFormatter('%H%Z')
        minor_locator = mpldates.AutoDateLocator(minticks=1, maxticks=4,
                                                 interval_multiples=True)
        minor_formatter = mpldates.DateFormatter('\n%d/%m/%Y')
    elif xaxis_range <= datetime.timedelta(days=14):
        #<=2 week range - display days and months, eg Sun 02 \n Mar 2014
        major_locator = mpldates.AutoDateLocator(minticks=3, maxticks=15,
                                                 interval_multiples=True)
        major_formatter = mpldates.DateFormatter('%a %d')
        minor_locator = mpldates.MonthLocator()
        minor_formatter = mpldates.DateFormatter('\n%b %Y')
    elif xaxis_range <= datetime.timedelta(days=62):
        #<= ~2 month range - display date and month, eg 02 \n Mar 2014
        major_locator = mpldates.AutoDateLocator(minticks=3, maxticks=63,
                                                 interval_multiples=True)
        major_formatter = mpldates.DateFormatter('%d')
        minor_locator = mpldates.MonthLocator()
        minor_formatter = mpldates.DateFormatter('\n%b %Y')
    elif xaxis_range <= datetime.timedelta(days=750):
        #<= ~2years range - display month and year, eg Mar \n 2014
        major_locator = mpldates.MonthLocator()
        major_formatter = mpldates.DateFormatter('%b')
        minor_locator = mpldates.AutoDateLocator(minticks=1, maxticks=3,
                                                 interval_multiples=True)
        minor_formatter = mpldates.DateFormatter('\n%Y')
    else:
        #> ~ 2years range - display year only eg 2014
        major_locator = mpldates.YearLocator()
        major_formatter = mpldates.DateFormatter('%Y')
        minor_locator = mpldates.YearLocator()
        minor_formatter = mpldates.DateFormatter('\n')

    #Can now set these new locator and formatter objects into the x-axis:
    xax.set_minor_locator(minor_locator)
    xax.set_minor_formatter(minor_formatter)
    xax.set_major_locator(major_locator)
    xax.set_major_formatter(major_formatter)

    return ax


#--- Class ---
class TimeSeriesPlot(line_plot.LinePlot):
    """
    Class for plotting time-series of cubes, a subclass of LinePlot.

    As an example, try loading observations and model data cubes:

    >>> import adaq_data
    >>> import adaq_functions
    >>> import config
    >>> import inifile
    >>> import sites_info
    >>> sample_data_path = config.SAMPLE_DATADIR+'sites_cube_list/'
    >>> OD = adaq_data.ADAQData()
    >>> scl = OD.load_ts(sample_data_path+'aurn_1days.nc')
    >>> obs_cube = OD.extract(short_name='O3', abbrev='HAR', singlecube=True)
    >>> print(obs_cube) # doctest: +ELLIPSIS
    mass_concentration_of_ozone_in_air / (ug/m3) (time: 25)
         Dimension coordinates:
              time                                    x
         Scalar coordinates:
              abbrev: HAR
              latitude: 51.57110977 degrees
              longitude: -1.326666594 degrees
              site_altitude: 137 m
              site_id: 35867333.14...
              site_name: Harwell
              site_type: RURAL
         Attributes:
              Conventions: CF-1.5
              label: Obs
              short_name: O3
              source: AURN
         Cell methods:
              mean: time (1 hour)

    >>> MD = adaq_data.ADAQData()
    >>> scl = MD.load_ts(sample_data_path+'aqum_oper_1days.nc')
    >>> mod_cube = MD.extract(short_name='O3', abbrev='HAR', singlecube=True)

    Load postprocessed aqum (sppo) data

    >>> inifilename = 'adaqcode/example_data_1days.ini'
    >>> ini_dict = inifile.get_inidict(defaultfilename=inifilename)
    ...  # doctest: +ELLIPSIS
    Reading inifile .../example_data_1days.ini

    Include required variables to load sppo data

    >>> ini_dict['models_fmt_list'] = ['nimrod']
    >>> ini_dict['models_dir_list'] = [config.SAMPLE_DATADIR +
    ... 'aqum_output/nimrod/']
    >>> ini_dict['models_list'] = ['aqeur_sppo']
    >>> ini_dict['forecast_day'] = 'latest'

    Get sites info

    >>> sites_data = sites_info.get_siteinfo(ini_dict)
    Number of sites:  5

    Load sppo data

    >>> sppo_md_list = adaq_functions.get_models(
    ...  ini_dict, sites_data) # doctest: +ELLIPSIS
    Getting model data for  aqeur_sppo  at  ...

    Extract required cubes

    >>> sppo_cube = sppo_md_list[0].extract(
    ...     site_name='Harwell', short_name='O3', singlecube=True)

    Set up time-series plotting class:

    >>> TSP = TimeSeriesPlot()

    Add cubes as individual lines:

    >>> TSP.add_line(obs_cube, label='obs', colour='k', linestyle='--')

    If linestyle not given, will default to solid line, while if label
    not given, will try to extract from cube attributes:

    >>> TSP.add_line(mod_cube, colour='r')
    >>> TSP.add_line(sppo_cube, colour='orange')

    If colour not given, will chose from given distinct list of colours,
    but one that has not yet been used:

    >>> TSP.add_line(mod_cube*0.5, label='mod*0.5')

    Produce plot

    >>> fig = TSP.plot()

    .. Note:: To display figure to screen, just use plt.show()

    Save figure. If filename not given, will chose sensible name.
    If plotdir not given, will default to current directory.

    >>> plotdir = config.CODE_DIR + "/adaqdocs/figures"
    >>> TSP.save_fig(plotdir=plotdir) # doctest: +ELLIPSIS
    Saved figure  .../Harwell_O3.png

    .. image:: ../adaqdocs/figures/Harwell_O3.png
      :scale: 75%

    To view all attributes now set in class, use tsp.__dict__
    """

    def __init__(self):
        """
        Initiates class as a subset of line_plot.LinePlot
        """

        line_plot.LinePlot.__init__(self)
        self.horiz_y = None #Value on y axis to plot a horizontal line at
        self.xticks = None #Override automatic x-axis ticks
        self.yticks = None #Override automatic y-axis ticks

    def plot(self):
        """
        Produce time-series plot.

        Returns fig object for further plotting if needed.

        """

        if not self.lines:
            raise ValueError("TimeSeriesPlot: no lines have been added")

        if self.fig is None:
            self.fig = plt.figure()
        ax = plt.gca()

        # Setup log y-axis if requested
        if self.semilogy:
            plt.semilogy()

        for line in self.lines:

            # Add extra settings into a dictionary
            add_settings = {}
            if 'add_settings' in line:
                add_settings = line['add_settings']

            #Check plotting all lines in same units
            if self.units == line['cube'].units:
                #Plot line
                iplt.plot(line['x'], line['y'], color=line['colour'],
                          label=line['label'], linestyle=line['linestyle'],
                          linewidth=line['linewidth'], marker=line['marker'],
                          **add_settings
                         )
            else:
                warnings.warn('Line not plotted - units of '+
                              str(line['cube'].units) +' for '+
                              str(line['cube'].name())+
                              ' do not match plot units of '+
                              str(self.units))
                # Units are a mismatch, so doesn't make sense to be on the same plot
                continue

        #Sort out dates on x-axis
        ax = set_xaxis_date_fmt(ax)

        #Set x/ylimits
        if self.xlim is not None:
            ax.set_xlim(self.xlim)
        if self.ylim is not None:
            ax.set_ylim(self.ylim)

        #Set x/yticks
        if self.xticks is not None:
            ax.set_xticks(self.xticks)
        if self.yticks is not None:
            ax.set_yticks(self.yticks)

        #Rotate labels on x-axis
        labels = ax.get_xticklabels(minor=True) #Minor x-axis labels
        for label in labels:
            label.set_fontsize('x-small')
        labels = ax.get_xticklabels()
        for label in labels:
            #label.set_rotation(30) #if want rotated at some point, uncomment
            label.set_fontsize('x-small')

        #Add horizontal y line (don't extend axis if not already shown)
        if self.horiz_y is not None:
            ylim = ax.get_ylim()
            if ylim[0] <= self.horiz_y <= ylim[1]:
                plt.axhline(self.horiz_y, color='k', linestyle='--')

        #Add legend
        if self.legend:
            plotting_functions.add_legend_belowaxes()

        #Add title and axis labels
        if self.title is None:
            self.title = self.gen_title()

        ax.set_title(self.title)
        if self.xlabel is None:
            self.xlabel = 'Time'
        ax.set_xlabel(self.xlabel)
        if self.ylabel is None:
            self.ylabel = self.gen_ylabel()
        ax.set_ylabel(self.ylabel)

        #Add gridlines
        if self.gridlines:
            plt.grid()

        # Apply branding
        if self.mobrand:
            line_plot.add_mobranding()

        return self.fig

    def gen_title(self):
        """
        Generate automatic title
        Uses attributes already set in class

        >>> TSP = TimeSeriesPlot()
        >>> TSP.site_name = 'Harwell'
        >>> TSP.site_type = 'RURAL'
        >>> TSP.phenomena_name = \
        'mass_concentration_of_pm10_ambient_aerosol_in_air'
        >>> title = TSP.gen_title()
        >>> print(title)
        Harwell (Rural)
        mass concentration of pm10 ambient aerosol in air
        """
        title = ''
        if self.site_name is not None:
            title += self.site_name.replace('_', ' ')
            if self.site_type is not None:
                title += ' ('
                title += self.site_type.replace('_', ' ').capitalize()
                title += ')'
            title += '\n'
        if self.phenomena_name is not None:
            title += self.phenomena_name.replace('_', ' ')

        return title

    def gen_ylabel(self):
        r"""
        Generate automatic label for y-axis
        Uses attributes already set in class

        >>> TSP = TimeSeriesPlot()
        >>> TSP.phenomena_short_name = 'O3'
        >>> TSP.units = 'ug/m3'
        >>> label = TSP.gen_ylabel()
        >>> print(label)
        O3 ($\mu g\ m^{-3}$)
        """
        ylabel = ''
        if self.phenomena_short_name is not None:
            ylabel += self.phenomena_short_name
        if self.units is not None:
            units = plotting_functions.units_str(self.units)
            if self.phenomena_short_name is not None:
                ylabel += ' ('
            ylabel += units
            if self.phenomena_short_name is not None:
                ylabel += ')'

        return ylabel

    def gen_filename(self):
        """
        Generate automatic filename (doesn't include path)
        Uses attributes already set in class

        >>> TSP = TimeSeriesPlot()
        >>> TSP.site_name = 'Harwell'
        >>> TSP.phenomena_short_name = 'O3'
        >>> filename = TSP.gen_filename()
        >>> print(filename)
        Harwell_O3.png
        """

        filename = ''
        if self.site_name is not None:
            if filename != '':
                filename += '_'
            filename += self.site_name
        if self.phenomena_short_name is not None:
            if filename != '':
                filename += '_'
            filename += self.phenomena_short_name

        #Add figure format
        filename += '.' + self.figformat

        return filename


def tsp_statistic(statscube_list, statistic, plotdir='./', filesuffix='',
                  colours_list=None, xcoordname=None):
    """
    Function to plot timeseries of a statistic, given a list of statistics
    cubes from :class:`TimeSeriesStats`.

    Note this routine is also called from the higher level
    :func:`adaq_plotting.plot_timeseries_of_stats`

    :param statscube_list: List of statistics cubes created by \
:meth:`timeseries_stats.TimeSeriesStats.convert_to_cube`
                           and then merged such that each cube has more than
                           one time coordinate.
    :param statistic: Name of statistic to plot. Should be one of the keys of
                      the :data:`timeseries_stats.STATS_INFO` dictionary
    :param plotdir: String - directory for plot to be saved.
    :param filesuffix: String to add to end of filename for specific naming
                       purposes. By default adds nothing.
    :param colours_list: List of colours to use to match in order against
                         cubes in statscube_list. If not set, defaults to
                         plotting_functions.COLOURS
    :param xcoordname: Name of coordinate to use on x-axis. If not set, defaults
                       to first dimension coordinate (usually time).

    Example usage:

    >>> import adaq_functions
    >>> ini_dict, stats_cubes = adaq_functions.get_exampledata(
    ... exampletype='stats') # doctest: +ELLIPSIS
    Reading inifile ...example_data_10days.ini
    >>> import config
    >>> plotdir = config.CODE_DIR + "/adaqdocs/figures"

    Extract a list of just the O3 cubes:

    >>> stats_cubes_o3 = []
    >>> for cubelist in stats_cubes:
    ...     stat_cube_o3 = cubelist.extract(iris.AttributeConstraint(
    ... short_name='O3'), strict=True)
    ...     stats_cubes_o3.append(stat_cube_o3)

    >>> print(stats_cubes_o3)
    [<iris 'Cube' of mass_concentration_of_ozone_in_air / (1) \
(istatistic: 32; time: 10)>]
    >>> print(stats_cubes_o3[0])
    mass_concentration_of_ozone_in_air / (1) (istatistic: 32; time: 10)
         Dimension coordinates:
              istatistic                                x         -
              time                                      -         x
         Auxiliary coordinates:
              statistic                                 x         -
              statistic_long_name                       x         -
              statistic_units                           x         -
         Attributes:
              Conventions: CF-1.5
              label: aqum_oper
              obs: AURN
              short_name: O3
    >>> print(stats_cubes_o3[0].coord('statistic').points) \
# doctest: +NORMALIZE_WHITESPACE
    ['mdi' 'nsites' 'npts' 'correlation' 'bias' 'nmb' 'mnmb' 'mge' 'nmge' \
'fge' 'rmse' 'fac2' 'ioa' 'threshold' 'orss' 'odds_ratio' 'hitrate' \
'falsealarmrate' 'falsealarmratio' 'o>=t_m>=t' 'o<t_m>=t' 'o>=t_m<t' \
'o<t_m<t' 'maxobs' 'maxmod' 'meanobs' 'meanmod' 'sdobs' 'sdmod' \
'perc_correct' 'perc_over' 'perc_under']

    Choose 'bias' from the list of available statistics and plot this:

    >>> tsp = tsp_statistic(stats_cubes_o3, 'bias', plotdir=plotdir)
    ... # doctest: +ELLIPSIS
    Saved figure  .../adaqdocs/figures/Timeseries_of_bias_O3.png

    .. image:: ../adaqdocs/figures/Timeseries_of_bias_O3.png
       :scale: 50%

    """

    if not colours_list:
        colours_list = plotting_functions.COLOURS[1:] #Ignore black

    #Information dictionary about statistic
    statsinfo = timeseries_stats.STATS_INFO[statistic]

    #Set up TimeSeriesPlot
    tsp = TimeSeriesPlot()

    short_name = ''

    #Loop through each cube, adding it as a line to the plot
    for icube, statscube in enumerate(statscube_list):
        #Extract statistic
        cube = statscube.extract(iris.Constraint(statistic=statistic))
        if cube is not None:
            short_name = cube.attributes['short_name']
            if sum(np.isfinite(cube.data)) >= 2:
                #Need at least 2 valid (non-nan) points to produce
                # sensible plot (To draw a line between two points)

                #Label according to model only
                label = cube.attributes['label']

                if xcoordname is not None:
                    xcoord = cube.coord(xcoordname)
                elif 'time' in [c.name() for c in cube.coords()]:
                    #Default to time coord if possible
                    xcoord = cube.coord('time')
                else:
                    xcoord = None

                tsp.add_line(cube, x=xcoord, label=label,
                             colour=colours_list[icube])

    if not tsp.lines:
        warnings.warn("No Timeseries of " + statistic + " plot for " + \
                      short_name + " (Not enough obs&model points)")
        #Don't plot this statistic
        return tsp

    #Add statistic and units to y axis label
    #Get units from last cube
    units = cube.coord('statistic_units').points[0]
    tsp.ylabel = statsinfo['long_name']
    if units is not None and units != '1':
        #Add units to ylabel, converting to latex units if in LATEX,
        #otherwise leaving as units.
        tsp.ylabel += ' (' + plotting_functions.units_str(
            units) + ')'

    #Also add a perfect-value line if appropriate
    perfect_value = statsinfo.get('perfect_value', None)
    if perfect_value is not None:
        tsp.horiz_y = perfect_value

    #Set up title
    tsp.title = 'Time series of ' + statsinfo['long_name']
    tsp.title += '\n ' + tsp.phenomena_name.replace('_', ' ')

    #Generate plot
    tsp.plot()

    #Save plot
    filename = 'Timeseries_of_' + statistic + '_' + \
               short_name + filesuffix+'.png'
    tsp.save_fig(plotdir=plotdir,
                 filename=filename)

    return tsp


if __name__ == '__main__':

    import doctest
    doctest.testmod()
