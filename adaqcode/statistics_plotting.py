"""
Various routines for plotting and displaying statistics
"""
from __future__ import division
from __future__ import print_function

from six.moves.builtins import zip

import warnings

import numpy as np
import matplotlib.pyplot as plt
import iris.coords
import iris.plot as iplt

import cube_time
import line_plot
import plotting_functions
import timeseries_stats

class Histogram(line_plot.LinePlot):
    """
    Class for plotting histograms, subclass of LinePlot.

    As an example, try loading observations and model data cubes:

    >>> import adaq_functions
    >>> ini_dict, sites_data, od, md_list = adaq_functions.get_exampledata(
    ... exampletype="full") # doctest: +ELLIPSIS
    Reading inifile .../example_data_5days.ini
    Number of sites:  5
    >>> obs_cube = od.extract(short_name='O3', singlecube=True)
    >>> mod_cube = md_list[0].extract(short_name='O3', singlecube=True)

    Set up Histogram plotting class:

    >>> Hist = Histogram()

    Add cubes as individual lines:

    >>> Hist.add_line(obs_cube,label='obs',colour='k')

    If colour not given, will choose from distinct list of colours

    >>> Hist.add_line(mod_cube)

    To manually set bin size:

    >>> Hist.binsize = 2.

    .. Note:: binsize will be set automatically otherwise

    Produce plot

    >>> fig = Hist.plot()

    .. Note:: To display figure to screen, just use plt.show()

    Save, using automatically generated filename to required directory

    >>> import config
    >>> plotdir = config.CODE_DIR + "/adaqdocs/figures"
    >>> Hist.save_fig(plotdir=plotdir) # doctest: +ELLIPSIS
    Saved figure  .../Histogram_O3.png

    .. image:: ../adaqdocs/figures/Histogram_O3.png
      :scale: 75%

    To check edges of bins:

    >>> np.set_printoptions(formatter={'float':lambda x: '{:3.0f}'.format(x)})
    >>> print(Hist.bin_edges)
    [  0   2   4   6   8  10  12  14  16  18  20  22  24  26  28  30  32  34
      36  38  40  42  44  46  48  50  52  54  56  58  60  62  64  66  68  70
      72  74  76  78  80  82  84  86  88  90  92  94  96  98 100 102 104 106
     108 110]
    >>> np.set_printoptions()

    The plots are automatically normalised, but this can be removed by:

    >>> Hist.normed = False

    Type of histogram plot can also be varied:

    >>> Hist.histtype = 'bar'

     * **step** [Default] generates a lineplot that is by default unfilled.
     * **bar** is a traditional bar-type histogram. If multiple data are given
       the bars are aranged side by side.
     * **barstacked** is a bar-type histogram where multiple data are
       stacked on top of each other.
     * **stepfilled** generates a lineplot that is by default filled.

    The maximum value displayed on the x-axis can be set based on a percentile
    of the data, for example setting to 99.5 will cut off the extreme maximum
    values:

    >>> Hist.maxperc = 99.5

    """

    def __init__(self):
        """
        Initiates class as a subset of line_plot.LinePlot,
        plus other histogram specific details
        """

        line_plot.LinePlot.__init__(self)
        self.binsize = None
        self.bin_edges = None
        self.alldata = None

        #Plotting options
        self.histtype = 'step' #Unfilled line-plot
        self.normed = True #Normalise histogram
        self.maxperc = 100 #Maximum percentile to include

    def __gen_alldata(self):
        """
        Get all data values and reshape to a 1D array
        """
        self.alldata = np.array((0))
        for line in self.lines:
            cube = line['cube']
            self.alldata = np.append(self.alldata, np.reshape(cube.data, -1))

    def gen_binsize(self):
        """
        Calculate sensible binsize if none set.
        This will generally be such that 10 bins are created,
        or for a larger dataset, creates 50 bins.
        These numbers should give a nice looking plot.
        Returns binsize.
        """
        if self.binsize is None:

            if self.alldata is None:
                self.__gen_alldata()

            minval = np.nanmin(self.alldata)
            maxval = np.nanmax(self.alldata)

            #May want to be modified for larger/smaller numbers
            # in future.
            if len(self.alldata) > 500:
                #50 bins
                self.binsize = (maxval - minval) / 50.
            else:
                #10 bins
                self.binsize = (maxval - minval) / 10.

        return self.binsize

    def gen_binedges(self):
        """
        Calculate edges of bins.
        All but the last (righthand-most) bin is half-open.
        In other words, if bins is:
        [1, 2, 3, 4]
        then the first bin is [1, 2) (including 1, but excluding 2)
        and the second [2, 3).
        The last bin, however, is [3, 4], which includes 4.

        >>> Hist = Histogram()
        >>> Hist.binsize = 2.
        >>> Hist.alldata = np.arange(10)
        >>> bin_edges = Hist.gen_binedges()
        >>> print(Hist.bin_edges) #doctest: +NORMALIZE_WHITESPACE
        [ 0.  2.  4.  6.  8. 10.]

        Returns bin_edges.
        """

        if self.binsize is None:
            self.gen_binsize()
        if self.alldata is None:
            self.__gen_alldata()

        minval = np.nanmin(self.alldata)
        maxval = np.nanmax(self.alldata)
        if minval == maxval:
            #Put some sensible numbers in either side of this value
            self.bin_edges = np.array([minval-0.5, minval+0.5])
            #Ensure that the max value (=min value) is plotted.
            self.maxperc = 100.
        else:
            #Max & min bin edges should be one binsize greater/less than
            #max/min values to ensure all data is plotted.
            bin_edges_min = minval - minval % self.binsize
            bin_edges_max = maxval + self.binsize
            #Generate bin edges
            self.bin_edges = np.arange(bin_edges_min, bin_edges_max, self.binsize)

        return self.bin_edges

    def plot(self):
        """
        Plot histogram.
        Returns figure object for further plotting if needed
        """

        if not self.lines:
            raise ValueError("Histogram: no lines have been added")

        if self.fig is None:
            self.fig = plt.figure()
        ax = self.fig.add_subplot(111)

        if self.bin_edges is None:
            self.gen_binedges()

        max_x = None

        for line in self.lines:
            #Get a 1d array of data
            data = np.reshape(line['cube'].data, -1)
            if sum(np.isfinite(data)) == 0:
                #All nan data, so don't include this line
                warnings.warn('All nan data')
                continue
            returned_tuple = ax.hist(data, bins=self.bin_edges,
                                     histtype=self.histtype, normed=self.normed,
                                     color=line['colour'], label=line['label'],
                                     linewidth=line['linewidth'],
                                     range=[np.nanmin(data), np.nanmax(data)])

            if self.maxperc != 100:
                cumulative_sum = np.cumsum(returned_tuple[0])
                #Get maximum value to be 1
                cumulative_sum /= cumulative_sum[-1]
                indices = np.where(cumulative_sum*100. > self.maxperc)
                if indices[0].size:
                    max_x_line = returned_tuple[1][indices[0][0]]
                    if max_x is None:
                        max_x = max_x_line
                    else:
                        max_x = np.max([max_x, max_x_line])

        #Add legend
        if self.legend:
            plotting_functions.add_legend_belowaxes()

        #Add title
        if self.title is None:
            self.gen_title()
        ax.set_title(self.title)

        #Add ylabel
        if self.ylabel is None:
            self.ylabel = 'Frequency'
        ax.set_ylabel(self.ylabel)

        #Add xlabel (units)
        if self.xlabel is None:
            if self.units is not None:
                units = plotting_functions.units_str(self.units)
                self.xlabel = units
            else:
                self.xlabel = ''
        ax.set_xlabel(self.xlabel)

        #Set maximum x limit
        if max_x is not None:
            ax.set_xlim(right=max_x)

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
        """
        self.title = 'Frequency Distribution'
        #Add short_name
        if self.phenomena_short_name is not None:
            self.title += ' of '+self.phenomena_short_name
        #Add date/time if same for all cubes
        startstr, endstr = self.get_startend_str()
        if startstr is not None and endstr is not None:
            self.title += '\n '+startstr+' to '+endstr

    def gen_filename(self):
        """
        Generate automatic filename
        Uses attributes already set in class

        >>> Hist = Histogram()
        >>> Hist.phenomena_short_name = 'O3'
        >>> filename = Hist.gen_filename()
        >>> print(filename)
        Histogram_O3.png
        """

        filename = 'Histogram'
        if self.phenomena_short_name is not None:
            filename += '_' + self.phenomena_short_name
        filename += '.' + self.figformat
        return filename

class SoccerPlot(line_plot.LinePlot):
    """
    Class for plotting Soccer Plots, subclass of LinePlot.

    **Example**

    Firstly, get some ADAQData...

    >>> import adaq_functions
    >>> ini_dict, sites_data, od, md_list = \
    adaq_functions.get_exampledata()  # doctest: +ELLIPSIS
    Reading inifile .../example_data_1days.ini
    Number of sites:  5
    >>> obs_cube = od.extract(short_name='O3', singlecube=True)
    >>> mod_cube = md_list[0].extract(short_name='O3', singlecube=True)
    >>> mod2_cube = md_list[1].extract(short_name='O3', singlecube=True)

    Now set up SoccerPlot class:

    >>> SP = SoccerPlot()

    Need to calculate statistics by comparing observations and model data:

    >>> stats_cube = SP.get_stats(obs_cube,mod_cube)

    This returns a statistics cube, which is the same as a model cube, but
    has two extra coordinates which contain statistics data:

    >>> print(stats_cube) # doctest: +ELLIPSIS
    mass_concentration_of_ozone_in_air / (ug/m3) (site_id: 5; time: 25)
         Dimension coordinates:
              site_id                                    x        -
              time                                       -        x
         Auxiliary coordinates:
              abbrev                                     x        -
              fge                                        x        -
              grid_latitude                              x        -
              grid_longitude                             x        -
              latitude                                   x        -
              longitude                                  x        -
              mnmb                                       x        -
              site_altitude                              x        -
              site_name                                  x        -
              site_type                                  x        -
              surface_altitude                           x        -
              forecast_period                            -        x
         Scalar coordinates:
              forecast_day: 1.0 Days
              level_height: 20.000... m, bound=(0.0, 49.998...) m
              model_level_number: 1
              sigma: 0.997..., bound=(1.0, 0.994...)
         Attributes:
              Conventions: CF-1.5
              STASH: m01s34i001
              label: aqum_oper
              short_name: O3
              source: Data from Met Office Unified Model
         Cell methods:
              mean: time (1 hour)

    To have a look at the values:

    >>> np.set_printoptions(formatter={'float':lambda x: '{:5.2f}'.format(x)})
    >>> print(stats_cube.coord('fge').points)
    [ 0.26  0.41  0.32  0.19  0.80]

    >>> print(stats_cube.coord('mnmb').points)
    [ 0.24 -0.41 -0.20 -0.08 -0.41]
    >>> np.set_printoptions()

    Now add this cube as a line (using the line_plot.LinePlot.add_line method:

    >>> SP.add_line(stats_cube)

    Similarly, can add a second set of model data:

    >>> stats_cube = SP.get_stats(obs_cube,mod2_cube)

    This time, give the line a different marker shape (a square) and a
    different colour:

    >>> SP.add_line(stats_cube,marker='s',colour='g')

    Now plot the data:

    >>> fig = SP.plot()

    Save, using automatically generated filename to required directory

    >>> import config
    >>> plotdir = config.CODE_DIR + "/adaqdocs/figures"
    >>> SP.save_fig(plotdir=plotdir) # doctest: +ELLIPSIS
    Saved figure  .../Soccer_Plot_O3.png

    .. image:: ../adaqdocs/figures/Soccer_Plot_O3.png
      :scale: 75%

    .. Note:: Code is set up to change which statistics are plotted on x and y
              axis by changing SP.stat_xaxis='rmse' etc. However no goals or
              ranges are currently defined automatically in get_goals() and
              get_ranges() - these should be added into these functions as
              required in the future.
    """

    def __init__(self):
        """
        Initiates class as a subset of line_plot.LinePlot,
        plus other histogram specific details
        """

        line_plot.LinePlot.__init__(self)
        self.stat_xaxis = 'mnmb'
        self.stat_yaxis = 'fge'
        self.stat_xgoal = None
        self.stat_ygoal = None

    def get_stats(self, obs_cube, mod_cube):
        """
        Produce a cube similar to the model cube, but with two
        extra coordinates which contain their statistics as points
        """

        #Only calculate statistics for periods where times match
        obs_cube, mod_cube = cube_time.intersect_cubetime(
            [obs_cube, mod_cube])

        stats_cube = mod_cube.copy()

        #Add statistic which will appear on x-axis as a coordinate
        stats_x_coord = stats_cube.coord('site_name').copy()
        stats_x_coord.rename(self.stat_xaxis)
        if len(stats_x_coord.points) > 1:
            #Assumes that site coordinates are in 1st dimension.
            #(is there a way to check this?)
            stats_cube.add_aux_coord(stats_x_coord, 0)
        else:
            #Scalar coordinate
            stats_cube.add_aux_coord(stats_x_coord)


        #Add statistic which will appear on y-axis as a coordinate
        stats_y_coord = stats_cube.coord('site_name').copy()
        stats_y_coord.rename(self.stat_yaxis)
        if len(stats_y_coord.points) > 1:
            stats_cube.add_aux_coord(stats_y_coord, 0)
        else:
            #Scalar coordinate
            stats_cube.add_aux_coord(stats_y_coord)


        stats_cube_list = iris.cube.CubeList()
        #Loop through individual sites and calculate their statistics
        for site_name in stats_cube.coord('site_name').points:
            obs_sitecube = obs_cube.extract(iris.Constraint(
                site_name=site_name))
            mod_sitecube = mod_cube.extract(iris.Constraint(
                site_name=site_name))

            #Check site has been found. If site is not found, it
            # is not included in returned stats_cube.
            if obs_sitecube is None:
                print('Site not found in obs_cube: ', site_name)
                continue
            if mod_sitecube is None:
                print('Site not found in mod_cube: ', site_name)
                continue

            stats_sitecube = stats_cube.extract(iris.Constraint(
                site_name=site_name))

            #Calculate statistics
            stats = timeseries_stats.TimeSeriesStats(obs_sitecube, mod_sitecube,
                                                     mdi=np.nan)

            #Add x-axis statistics as a coordinate point
            # Use getattr to call method from stats
            # with the name self.stat_xaxis
            stats_x = getattr(stats, self.stat_xaxis)()
            stats_sitecube.coord(self.stat_xaxis).points = stats_x

            #Add y-axis statistics as a coordinate point
            stats_y = getattr(stats, self.stat_yaxis)()
            stats_sitecube.coord(self.stat_yaxis).points = stats_y

            stats_cube_list.append(stats_sitecube)
        stats_cube = stats_cube_list.merge_cube()

        return stats_cube


    def plot(self, legendcols=None):
        """
        Produce plot.
        :param legendcols: Number of columns to include in legend.
        """

        if not self.lines:
            raise ValueError("SoccerPlot: no lines have been added")

        if self.fig is None:
            self.fig = plt.figure()
        ax = self.fig.add_subplot(111)


        #Scatter Plot
        for line in self.lines:
            if line['marker'] is None:
                line['marker'] = 'o' #Default value for a scatter plot
            #Check has coordinates
            if line['cube'].coords(self.stat_xaxis) and \
               line['cube'].coords(self.stat_yaxis):
                #Check that data are not all NaNs:
                xaxispts = line['cube'].coord(self.stat_xaxis).points
                yaxispts = line['cube'].coord(self.stat_yaxis).points
                if not np.isnan(np.nanmax(xaxispts)) and \
                   not np.isnan(np.nanmax(yaxispts)):
                    iplt.scatter(line['cube'].coord(self.stat_xaxis),
                                 line['cube'].coord(self.stat_yaxis),
                                 color=np.atleast_1d(line['colour']),
                                 label=line['label'],
                                 marker=line['marker'],
                                 edgecolor='k', s=30
                                )
            else:
                raise ValueError("Cube does not have statistics coordinates \n"+
                                 "May need to run get_stats() first")

        #Set x & y axis limits
        range_x = self.get_range(self.stat_xaxis)
        if range_x is not None:
            ax.set_xlim(range_x)

        range_y = self.get_range(self.stat_yaxis)
        if range_y is not None:
            ax.set_ylim(range_y)

        #Plot goal regions
        if self.stat_xgoal is None:
            #Get goal if not already set
            self.stat_xgoal = self.get_goals(self.stat_xaxis)
        if self.stat_ygoal is None:
            #Get goal if not already set
            self.stat_ygoal = self.get_goals(self.stat_yaxis)
        if self.stat_xgoal is not None and self.stat_ygoal is not None:
            #Can plot goals - plot as a square
            for xgoal, ygoal in zip(self.stat_xgoal, self.stat_ygoal):
                xpoints = [xgoal, -xgoal, -xgoal, xgoal, xgoal]
                ypoints = [ygoal, ygoal, -ygoal, -ygoal, ygoal]
                ax.plot(xpoints, ypoints, 'k--')

        #Add lines through zero
        ax.plot(ax.get_xlim(), [0, 0], 'k')
        ax.plot([0, 0], ax.get_ylim(), 'k')

        #Add legend
        if self.legend:
            if legendcols is None:
                plotting_functions.add_legend_belowaxes(scatterpoints=1)
            else:
                plotting_functions.add_legend_belowaxes(scatterpoints=1,
                                                        ncol=legendcols)

        #Add title
        if self.title is None:
            self.gen_title()
        ax.set_title(self.title)

        #Add x and y labels
        if self.xlabel is None:
            if self.stat_xaxis in timeseries_stats.STATS_INFO:
                self.xlabel = timeseries_stats.STATS_INFO[self.stat_xaxis][
                    'long_name']
            else:
                self.xlabel = self.stat_xaxis
            ax.set_xlabel(self.xlabel)

        if self.ylabel is None:
            if self.stat_yaxis in timeseries_stats.STATS_INFO:
                self.ylabel = timeseries_stats.STATS_INFO[self.stat_yaxis][
                    'long_name']
            else:
                self.ylabel = self.stat_yaxis
            ax.set_ylabel(self.ylabel)

        #Add gridlines
        if self.gridlines:
            plt.grid()

        # Apply branding
        if self.mobrand:
            line_plot.add_mobranding()

        return self.fig


    def get_range(self, statname):
        """
        Get the range to display for each statistic

        >>> SP = SoccerPlot()
        >>> print(SP.get_range('fge'))
        [0.0, 1.75]

        If statistic is not yet included:

        >>> print(SP.get_range('rmse'))
        None
        """

        ranges = {'mnmb': [-1.50, 1.50],
                  'fge' : [0.0, 1.75]}

        if statname in ranges:
            stat_range = ranges[statname]
        else:
            stat_range = None

        return stat_range

    def get_goals(self, statname):
        """
        Get goal critera for each statistic

        >>> SP = SoccerPlot()
        >>> print(SP.get_goals('mnmb'))
        [0.15, 0.3, 0.6]
        """

        goals = {'mnmb':[0.15, 0.30, 0.60],
                 'fge' :[0.35, 0.50, 0.75]}

        if statname in goals:
            stat_goals = goals[statname]
        else:
            stat_goals = None

        return stat_goals

    def gen_title(self):
        """
        Generate automatic title
        Uses attributes already set in class
        """
        self.title = 'Soccer Plot'
        #Add short_name
        if self.phenomena_short_name is not None:
            self.title += ' of '+self.phenomena_short_name
        #Add date/time if same for all cubes
        startstr, endstr = self.get_startend_str()
        if startstr is not None and endstr is not None:
            self.title += '\n '+startstr+' to '+endstr

    def gen_filename(self):
        """
        Generate automatic filename
        Uses attributes already set in class

        >>> SP = SoccerPlot()
        >>> SP.phenomena_short_name = 'O3'
        >>> filename = SP.gen_filename()
        >>> print(filename)
        Soccer_Plot_O3.png

        """

        filename = 'Soccer_Plot'
        if self.phenomena_short_name is not None:
            filename += '_' + self.phenomena_short_name
        filename += '.' + self.figformat
        return filename


class QQPlot(line_plot.LinePlot):
    """
    Class for plotting quantile-quantile plots, subclass of LinePlot.

    A quantile-quantile plot takes two distributions (or sets of data)
    and plots the quantiles from each distribution against each other.
    For example the 10th percentile value from one distribution might
    be 3.0,which is plotted on the x-axis, against the 10th percentile
    value from the second distribution (eg 3.5) plotted on the y-axis.
    If both datasets come from the same distribution, it would be
    expected for the plotted values to fall on the 1-1 line.

    As an example, try loading observations and model data cubes:

    >>> import adaq_functions
    >>> ini_dict, sites_data, od, md_list = \
    adaq_functions.get_exampledata()  # doctest: +ELLIPSIS
    Reading inifile .../example_data_1days.ini
    Number of sites:  5
    >>> obs_cube = od.extract(short_name='O3', singlecube=True)
    >>> mod_cube = md_list[0].extract(short_name='O3', singlecube=True)
    >>> mod2_cube = md_list[1].extract(short_name='O3', singlecube=True)

    Set up Quantile-Quantile plotting class:

    >>> qq = QQPlot()

    The data that is required for the x-axis needs to be added separately.
    This can be done either on class initialisation:

    >>> qq = QQPlot(xcube = obs_cube)

    Or later:

    >>> qq.xcube = obs_cube

    The data required for the y-axis can then be added as individual lines,
    and will all compare against the same data on the x-axis:

    >>> qq.add_line(mod_cube)
    >>> qq.add_line(mod2_cube, marker='+')

    Can then produce the plot:

    >>> fig = qq.plot()

    Save to file

    >>> import config
    >>> qq.save_fig(plotdir=config.CODE_DIR + "/adaqdocs/figures")
    ... # doctest: +ELLIPSIS
    Saved figure  .../adaqdocs/figures/Quantile-Quantile_O3.png

    .. image:: ../adaqdocs/figures/Quantile-Quantile_O3.png
      :scale: 75%

    """
    def __init__(self, xcube=None):
        """
        Initiates class as a subset of line_plot.LinePlot,
        plus other quantile-quantile specific details
        """

        line_plot.LinePlot.__init__(self)
        self.quantiles = list(np.arange(101)) #List of quantiles
        self.xcube = xcube
        self.xpercentiles = None
        self.one2one = True #Plot 1-1 line


    def get_percentiles(self):
        """
        Calculate percentiles for x-axis values, and for each line data.
        Removes any nans before percentiles are calculated.
        """

        if self.xcube is None:
            self.xpercentiles = self.quantiles
        else:
            data = self.xcube.data
            if np.sum(np.isfinite(data)) == 0:
                #All data is nan, cannot calculate percentiles
                raise ValueError('Data for x-axis all nan')
            self.xpercentiles = np.percentile(data[~np.isnan(data)],
                                              q=self.quantiles)

        for line in self.lines:
            data = line['cube'].data
            if np.sum(np.isfinite(data)) == 0:
                #All data is nan, cannot calculate percentiles
                warnings.warn('All nan data')
                continue
            line['ypercentiles'] = np.percentile(data[~np.isnan(data)],
                                                 q=self.quantiles)

        return self.lines

    def plot(self):
        """
        Plot quantile-quantile plot.
        Returns figure object for further plotting if needed.
        """

        if not self.lines:
            raise ValueError("Quantile-Quantile plot: no lines have been added")

        self.get_percentiles()

        self.fig = plt.figure()
        ax = plt.gca()

        axis_max = np.nanmax(self.xpercentiles)
        axis_min = np.nanmin(self.xpercentiles)

        for line in self.lines:

            if line['marker'] is None:
                #Ensure a marker is set
                line['marker'] = 'o'

            if 'ypercentiles' not in line:
                #ypercentiles were not calculated due to all nan data
                #therefore cannot plot this line.
                continue

            plt.scatter(self.xpercentiles, line['ypercentiles'],
                        color=np.atleast_1d(line['colour']), label=line['label'],
                        linewidth=line['linewidth'], marker=line['marker']
                       )
            axis_max = np.nanmax([axis_max, np.nanmax(line['ypercentiles'])])
            axis_min = np.nanmin([axis_min, np.nanmin(line['ypercentiles'])])

        #Set x/y lims
        if axis_min > 0 and axis_min < axis_max / 25.:
            #Ensure zero is plotted if close to zero.
            axis_min = 0
        ax.set_xlim([axis_min, axis_max])
        ax.set_ylim([axis_min, axis_max])

        #Add 1-1 line
        if self.one2one:
            plt.plot([axis_min, axis_max], [axis_min, axis_max],
                     color='k', linestyle='-')

        #Add title
        if self.title is None:
            self.gen_title()
        ax.set_title(self.title)

        #Add axis labels
        if self.xlabel is None:
            self.xlabel = self.xcube.attributes['label']
            if self.units is not None:
                units = plotting_functions.units_str(self.units)
                self.xlabel += ' ('+units+')'
        ax.set_xlabel(self.xlabel)

        if self.ylabel is None:
            if self.units is not None:
                units = plotting_functions.units_str(self.units)
            if len(self.lines) == 1:
                self.ylabel = self.lines[0]['cube'].attributes['label']
                if self.units is not None:
                    self.ylabel += ' ('+units+')'
            else:
                self.ylabel = self.phenomena_short_name +' ('+units+')'
        ax.set_ylabel(self.ylabel)

        #Add legend if needed
        if len(self.lines) > 1:
            plotting_functions.add_legend_belowaxes(scatterpoints=1)

        #Add gridlines
        if self.gridlines:
            plt.grid()

        # Apply branding
        if self.mobrand:
            line_plot.add_mobranding()

        return self.fig


    def gen_title(self):
        """
        Generate automatic title.
        Uses attributes already set in class.
        """
        self.title = 'Quantile-Quantile Plot'
        #Add short_name
        if self.phenomena_short_name is not None:
            self.title += ' of '+self.phenomena_short_name
        #Add date/time if same for all cubes
        startstr, endstr = self.get_startend_str()
        if startstr is not None and endstr is not None:
            self.title += '\n '+startstr+' to '+endstr


    def gen_filename(self):
        """
        Generate automatic filename.
        Uses attributes already set in class.

        >>> qq = QQPlot()
        >>> qq.phenomena_short_name = 'O3'
        >>> filename = qq.gen_filename()
        >>> print(filename)
        Quantile-Quantile_O3.png
        """

        filename = 'Quantile-Quantile'
        if self.phenomena_short_name is not None:
            filename += '_' + self.phenomena_short_name
        filename += '.' + self.figformat
        return filename


if __name__ == '__main__':

    import doctest
    doctest.testmod()
