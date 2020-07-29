'''
Includes class to set up and calculate common time-series statistics.
Also includes routines to deal with multiple sets of model data
including saving statistics to file (save_stats)

.. to do::
    Rescale model and observations to same precision before calculating
    statistics
'''
from __future__ import division
from __future__ import print_function

from six.moves.builtins import zip
from six.moves.builtins import str
from six.moves.builtins import range
from six.moves.builtins import object

import os
import warnings
from operator import itemgetter
import copy

import scipy.stats
import numpy as np
import numpy.ma as ma
import iris
import cf_units

import plotting_functions
import cube_time

#: Information about each statistic in the form of a dictionary for each
#: statistic.
#: Each dictionary must contain 'long_name' and 'istat'.
#: NOTE: the istat number must NOT be changed for each statistic.
#: - any new statistics should have a new unique istat number which has not been
#: used previously for any other statistic.
#: If units = 'units', then this means use the same units as the input data
#: Also available for some statistics (where appropriate):
#: 'perfect_value', 'min_value', 'max_value' which gives the 'best' score
#: and the theoretic range of values for this statistic.
#: 'ranking_method' should be given to indicate how to produce the 'best' score
#: options are 'max', 'min', 'absmax', 'absmin'
STATS_INFO = {
    'mdi': {'long_name': 'Missing data indicator',
            'istat': 1},
    'nsites': {'long_name': 'Number of sites',
               'min_value': 0,
               'ranking_method': 'max',
               'istat': 2},
    'npts': {'long_name': 'Number of points',
             'min_value': 0,
             'ranking_method': 'max',
             'istat': 3},
    'correlation': {'long_name': 'Correlation',
                    'perfect_value': 1.0,
                    'min_value': -1.0,
                    'max_value': 1.0,
                    'ranking_method': 'max',
                    'istat': 4},
    'bias': {'long_name':'Bias',
             'units': 'units',
             'perfect_value': 0.0,
             'ranking_method': 'absmin',
             'istat': 5},
    'nmb': {'long_name':'Normalised Mean Bias',
            'perfect_value': 0.0,
            'min_value': -1.0,
            'max_value': 1.0,
            'ranking_method': 'absmin',
            'istat': 6},
    'mnmb': {'long_name':'Modified Normalised Mean Bias',
             'perfect_value': 0.0,
             'min_value': -2.0,
             'max_value': 2.0,
             'ranking_method': 'absmin',
             'istat': 7},
    'mge': {'long_name': 'Mean Gross Error',
            'units': 'units',
            'perfect_value': 0.0,
            'min_value': 0.0,
            'ranking_method': 'min',
            'istat': 8},
    'nmge': {'long_name':'Normalised Mean Gross Error',
             'perfect_value': 0.0,
             'min_value': 0.0,
             'max_value': 1.0,
             'ranking_method': 'min',
             'istat': 9},
    'fge': {'long_name':'Fractional Gross Error',
            'perfect_value': 0.0,
            'min_value': 0.0,
            'max_value': 2.0,
            'ranking_method': 'min',
            'istat': 10},
    'rmse': {'long_name':'Root Mean Square Error',
             'units': 'units',
             'perfect_value': 0.0,
             'min_value': 0.0,
             'ranking_method': 'min',
             'istat': 11},
    'fac2': {'long_name':'Factor of 2',
             'perfect_value': 1.0,
             'min_value': 0.0,
             'max_value': 1.0,
             'ranking_method': 'max',
             'istat': 12},
    'ioa': {'long_name':'Index of Agreement',
            'perfect_value': 1.0,
            'min_value': -1.0,
            'max_value': 1.0,
            'ranking_method': 'max',
            'istat': 13},
    'threshold': {'long_name': 'Threshold',
                  'units': 'units',
                  'istat': 14},
    'orss': {'long_name': 'Odds Ratio Skill Score',
             'perfect_value': 1.0,
             'min_value': -1.0,
             'max_value': 1.0,
             'ranking_method': 'max',
             'istat': 15},
    'odds_ratio': {'long_name': 'Odds Ratio',
                   'min_value': 0.0,
                   'ranking_method': 'max',
                   'istat': 16},
    'hitrate': {'long_name': 'Hitrate',
                'perfect_value': 1.0,
                'min_value': 0.0,
                'max_value': 1.0,
                'ranking_method': 'max',
                'istat': 17},
    'falsealarmrate': {'long_name': 'False Alarm Rate',
                       'perfect_value': 0.0,
                       'min_value': 0.0,
                       'max_value': 1.0,
                       'ranking_method': 'min',
                       'istat': 18},
    'falsealarmratio': {'long_name': 'False Alarm Ratio',
                        'perfect_value': 0.0,
                        'min_value': 0.0,
                        'max_value': 1.0,
                        'ranking_method': 'min',
                        'istat': 19},
    'o>=t_m>=t': {'long_name': 'Number Obs >= Threshold and Model >= Threshold',
                  'istat': 20},
    'o<t_m>=t': {'long_name': 'Number Obs < Threshold and Model >= Threshold',
                 'istat': 21},
    'o>=t_m<t': {'long_name': 'Number Obs >= Threshold and Model < Threshold',
                 'istat': 22},
    'o<t_m<t': {'long_name': 'Number Obs < Threshold and Model < Threshold',
                'istat': 23},
    'maxobs': {'long_name': 'Maximum Observation Value',
               'units': 'units',
               'istat': 24},
    'maxmod': {'long_name': 'Maximum Model Value',
               'units': 'units',
               'istat': 25},
    'meanobs': {'long_name': 'Mean Observation Value',
                'units': 'units',
                'istat': 26},
    'meanmod': {'long_name': 'Mean Model Value',
                'units': 'units',
                'istat': 27},
    'sdobs': {'long_name': 'Standard Deviation of Observations',
              'units': 'units',
              'istat': 28},
    'sdmod': {'long_name': 'Standard Deviation of Model',
              'units': 'units',
              'istat': 29},
    'perc_correct': {'long_name': 'Percentage of Correct values',
                     'units': '%',
                     'perfect_value': 100.0,
                     'istat': 30},
    'perc_over': {'long_name': 'Percentage of Over-predictions',
                  'units': '%',
                  'perfect_value': 0.0,
                  'istat': 31},
    'perc_under': {'long_name':'Percentage of Under-predictions',
                   'units': '%',
                   'perfect_value': 0.0,
                   'istat': 32},
    'percobs': {'long_name': '{n}th Percentile Observation Value',
                'units': 'units',
                'istat': 33},
    'percmod': {'long_name': '{n}th Percentile Model Value',
                'units': 'units',
                'istat': 34},
    'units': {'long_name': 'Units',
              'istat': 35}
    }




class TimeSeriesStats(object):
    """
    Class to set up and for calculating basic time-series statistics

    >>> modellist = [131.12137, 131.65209, 94.16319, 121.12923, 139.73375,
    ...     104.43431, 145.28893, 95.45047]
    >>> obslist   = [None,           None, 99.80224, 133.56458,  95.24243,
    ...     107.61688, 211.10431, 115.73505]

    Optional input:
      * threshold - value to use for ORSS, hitrate and FAR calculations
      * mdi - value given to any missing data in modellist and obslist input
        data eg -999.
        Note None and np.nan is always assumed an mdi as well.
        Same value given to any undefined stats
      * label - string identifier to identify obs/model combination
      * phenomenon - string identifier to identify the phenomenon whose
        statistics are being considered

    >>> stats = TimeSeriesStats(obslist,modellist,mdi=None,threshold=100.)

    Get a single statistic:

    >>> correlation = stats.correlation()
    >>> print('{:.3f}'.format(correlation))
    0.580

    Get a full dictionary of common basic statistics:

    >>> statsdict = stats.get_basic_stats(percentile_value=95)
    >>> statsdict == {'sdmod': 22.259362602603872, 'pc95obs': 191.7193775,
    ... 'mdi': None, 'bias': -10.47760166666667, 'nmge': 0.19899766401439542,
    ... 'threshold': 100.0, 'meanobs': 127.17758166666665, 'orss': 0.5,
    ... 'perc_over': 16.666666666666664, 'mnmb': -0.06142778873417744,
    ... 'maxobs': 211.10431, 'pc95mod': 143.90013499999998, 'o>=t_m<t': 1,
    ... 'odds_ratio': 3.0, 'fac2': 1.0, 'hitrate': 0.75, 'o<t_m<t': 1,
    ... 'maxmod': 145.28893, 'o<t_m>=t': 1, 'mge': 25.308041666666668,
    ... 'sdobs': 43.283496393941846, 'meanmod': 116.69997999999998,
    ... 'nmb': -0.08238560231573311, 'perc_under': 83.33333333333334,
    ... 'falsealarmrate': 0.5, 'fge': 0.18765709418973464,
    ... 'ioa': 0.5796645327225636, 'rmse': 33.95872875514507,
    ... 'perc_correct': 0.0, 'falsealarmratio': 0.25, 'o>=t_m>=t': 3,
    ... 'correlation': 0.5795978100474246, 'npts': 6}
    True

    Input can also be in the form of iris cubes:

    >>> import adaq_functions
    >>> ini_dict, sites_data, od, md_list = adaq_functions.get_exampledata()
    ... # doctest: +ELLIPSIS
    Reading inifile .../example_data_1days.ini
    Number of sites:  5
    >>> obs_cube = od.extract(short_name='O3', singlecube=True)

    >>> mod_cube = md_list[0].extract(short_name='O3', singlecube=True)
    >>> stats = TimeSeriesStats(obs_cube, mod_cube)

    Calculate just correlation and bias

    >>> corr = stats.correlation()
    >>> bias = stats.bias()

    Check that values are as expected (albeit possibly in a different order)
    >>> expected_stats = {'mdi': None, 'nsites': 5, 'units': 'ug/m3',
    ... 'bias': -6.204036978363991, 'correlation': 0.64052159576610845,
    ... 'threshold': None, 'npts': 125}
    >>> for kv_pair in stats.statsdict:
    ...    assert(kv_pair in expected_stats)

    # >>> stats.statsdict == {'mdi': None, 'nsites': 5, 'units': 'ug/m3',
    # ... 'bias': -6.204036978363991, 'correlation': 0.64052159576610845,
    # ... 'threshold': None, 'npts': 125}
    # True

    As the input to TimeSeriesStats was two cubes, the output statistics
    can be converted into a cube:

    >>> cube = stats.convert_to_cube()
    >>> print(cube)
    mass_concentration_of_ozone_in_air / (1) (istatistic: 6; time: 1)
         Dimension coordinates:
              istatistic                                x        -
              time                                      -        x
         Auxiliary coordinates:
              statistic                                 x        -
              statistic_long_name                       x        -
              statistic_units                           x        -
         Attributes:
              label: aqum_oper
              obs: AURN
              short_name: O3

    This has a time coordinate - taken from meaning over the observation
    cube times:

    >>> print(cube.coord('time'))  # doctest: +ELLIPSIS
    DimCoord([2014-04-02 11:30:00], \
bounds=[[2014-04-01 23:00:00, 2014-04-03 00:00:00]], \
standard_name=...'time', calendar=...'gregorian', var_name='time_O3')

    Given on the istatistic dimension are the statistic details:

    >>> print(cube.coord('statistic').points)
    ['mdi' 'nsites' 'npts' 'correlation' 'bias' 'threshold']
    >>> print(cube.coord('statistic_long_name').points)
    ['Missing data indicator' 'Number of sites' 'Number of points'
     'Correlation' 'Bias' 'Threshold']
    >>> print(cube.coord('statistic_units').points)
    ['1' '1' '1' '1' 'ug/m3' 'ug/m3']

    To access the data for each 'bias':

    >>> bias_cube = cube.extract(iris.Constraint(statistic='bias'))
    >>> print(bias_cube)
    mass_concentration_of_ozone_in_air / (1) (time: 1)
         Dimension coordinates:
              time                                x
         Scalar coordinates:
              istatistic: 5
              statistic: bias
              statistic_long_name: Bias
              statistic_units: ug/m3
         Attributes:
              label: aqum_oper
              obs: AURN
              short_name: O3
    >>> np.set_printoptions(precision=3)
    >>> print(bias_cube.data)
    [-6.204]
    >>> np.set_printoptions()

    """

    def __init__(self, obs, mod, mdi=None, threshold=None,
                 label=None, phenomenon=None, stats_tcoord=None):
        """
        Initialise class with observations and model data, converting this
        data to np.arrays.

        :param obs: Observation data, in the form of an iris cube, numpy array
                    or a list.
        :param mod: Model data, in the form of an iris cube, numpy array
                    or a list.
        :param mdi: Missing Data Indicator - all data equal to this is removed
                    from model and observation arrays. Also any undefined
                    statistics are given this value.
        :param threshold: Threshold value to be used for threshold-based
                          statistics
        :param label: String identifier to identify obs/model combination
        :param phenomenon:  String identifier to identify the phenomenon whose
                            statistics are being considered
        :param stats_tcoord: iris DimCoord representing time, which can be used
                             to overwrite default setting of the time coordinate
                             in the statistics cubes and which will be therefore
                             be used for creating netcdf statistics files from.
        """

        self.obs = None
        self.model = None
        self.statsdict = {'mdi':mdi,
                          'threshold':threshold}
        self.label = label
        self.phenomenon = phenomenon
        self.n = 0
        self.diff = None
        self.nzobs = None
        self.count_nzobs = None
        self.nz = None
        self.count_nz = None
        self.statscube = None
        self.tcoord = None
        self.attributes = None

        #If both obs and mod are iris cubes, then use these to prepare for
        #creating a new stats cube.
        if isinstance(obs, iris.cube.Cube) and isinstance(mod, iris.cube.Cube):
            #Prepare for creating a new cube
            if stats_tcoord is None:
                if len(obs.coord('time').points) > 1:
                    #If multiple time points, then it is sensible to choose the
                    #mean of these points to store as the time in the new cube.
                    #It also preserves the max and min times through the
                    #resulting time bounds.
                    self.tcoord = obs.collapsed(
                        'time', iris.analysis.MEAN).coord('time')
                else:
                    self.tcoord = obs.coord('time')
            else:
                self.tcoord = copy.deepcopy(stats_tcoord)
            #Prepare for saving to netcdf:
            #Give time variable name a unique name defined by the short_name
            #Otherwise if multiple times in netcdf file for different
            #short_names, all covering different ranges, then they get renamed
            #as time, time_0, time_1 etc and can get confused between the
            #different cubes.
            self.tcoord.var_name = 'time'+'_'+obs.attributes['short_name']
            #Save attributes that will be required to build cube from
            self.attributes = {'short_name': obs.attributes['short_name'],
                               'obs': obs.attributes['source'],
                               'label': mod.attributes['label'],
                               'long_name': obs.name(),
                               'units': obs.units}

        #Check input types and convert to 1d list
        #Input is in cubes:
        if isinstance(obs, iris.cube.Cube):
            #Can use cube information to calculate number of sites
            # - Check each site and check has some non-nan data
            self.statsdict['nsites'] = 0
            for subcube in obs.slices_over(['site_id']):
                #print subcube
                maxval = np.nanmax(subcube.data)
                if not np.isnan(maxval) and maxval != self.statsdict['mdi']:
                    self.statsdict['nsites'] += 1

            obs = obs.data #Converts to np.array, which is converted below

        if isinstance(mod, iris.cube.Cube):
            #Use cube information to extract units and label
            self.statsdict['units'] = mod.units.__str__()
            if self.label is None:
                self.label = mod.attributes['label']
            mod = mod.data #Converts to np.array, which is converted below

        #Input is in np arrays:
        if isinstance(obs, np.ndarray):
            obs = list(obs.reshape(obs.size))
        if isinstance(mod, np.ndarray):
            mod = list(mod.reshape(mod.size))

        #Check now a list
        if not isinstance(obs, list):
            raise ValueError("obs is not a list")
        if not isinstance(mod, list):
            raise ValueError("mod is not a list")
        obslist = obs
        modellist = mod
        if len(obslist) != len(modellist):
            raise ValueError("len(obs) != len(mod)")

        #Remove None, np.nan values and mdi values
        obs = []
        model = []
        for i in range(len(obslist)):
            if obslist[i] is not None and modellist[i] is not None:
                if not np.math.isnan(obslist[i]):
                    if not np.math.isnan(modellist[i]):
                        if obslist[i] != mdi and modellist[i] != mdi:
                            obs.append(obslist[i])
                            model.append(modellist[i])

        self.obs = np.asarray(obs)
        self.model = np.asarray(model)
        self.n = float(len(self.obs))
        self.statsdict['npts'] = int(self.n)

        if self.n == 0:
            warnings.warn('No valid data')

    def get_basic_stats(self, percentile_value=95):
        """
        Calculate all basic time-series statistics and return dictionary of them
        """
        self.correlation()
        self.bias()
        self.nmb()
        self.mnmb()
        self.mge()
        self.nmge()
        self.fge()
        self.rmse()
        self.fac2()
        self.ioa()
        self.contingency_table()
        self.maxobs()
        self.maxmod()
        self.meanobs()
        self.meanmod()
        self.sdobs()
        self.sdmod()
        self.percentage()
        self.percentileobs(percentile_value)
        self.percentilemod(percentile_value)
        return self.statsdict

    def convert_to_cube(self):
        """
        Convert self.statsdict into an iris Cube, stored in self.statscube.
        The resulting cube will have two dimensions.

        The first dimension coordinate is istatistic. This contains a unique
        integer number for each statistic which is in self.statsdict. This
        ensures we can create a DimCoord which needs strictly monotonic points.
        Also stored on this same dimension are other coordinates whose points
        relate to statistics.

        The second dimension is time. This has a single point in it (which is
        generated from the mean of the input cubes as part of __init__).
        Having time as a dimension however will allow it to be concenated
        with other statistics cubes for other times during other code.

        The values in the cube correspond to the statistic values for each
        statistic in the istatistic dimension and for the time in the time
        dimension.
        """

        if self.tcoord is None:
            warnings.warn('Cannot convert to cube: input data are not cubes')
            return self.statscube #None

        #Get list of statistics keys, such that they will be in the order
        #to produce a monotonic list of istat values
        istat_stat = {}
        for stat in self.statsdict:
            istat_stat[stat] = STATS_INFO[stat]['istat']
        stats = []
        for k, _v in sorted(list(istat_stat.items()), key=itemgetter(1)):
            stats.append(k)

        #Get lists of the values and coordinate points
        #for each statistic, by looping over the list of statistic names
        #which are now in istat order.
        istat_points = []
        stat_points = []
        stat_longname_points = []
        stat_units_points = []
        values = []
        #Loop over statistic names
        for stat in stats:
            if stat == 'units':
                #Don't include units in the list
                #This is the units for the overall input data
                #The units for the individual statistic may be different
                #so are added below depending on the information about
                #units in STATS_INFO.
                continue
            #Add unique integer number
            istat_points.append(STATS_INFO[stat]['istat'])
            #Add statistic name
            stat_points.append(stat)
            #Add statistic long name
            stat_longname_points.append(STATS_INFO[stat]['long_name'])
            #Add statistic units (this depends on the statistic itself,
            # for example might be normalised, %, or raw data units)
            stat_units = STATS_INFO[stat].get('units', '1')
            if stat_units == 'units':
                stat_units = self.attributes['units']
            stat_units_points.append(str(stat_units))
            #Add values for data array
            if self.statsdict[stat] is None:
                values.append(np.nan) #As None not good for cubes
            else:
                values.append(self.statsdict[stat])
        #Now set up netcdf var_name and the coordinate for each required
        #coordinate
        istats_var_name = 'istatistic'+'_'+self.attributes['short_name']
        istats_coord = iris.coords.DimCoord(
            istat_points,
            long_name='istatistic',
            var_name=istats_var_name,
            units=1)
        stats_var_name = 'statistic'+'_'+self.attributes['short_name']
        stats_coord = iris.coords.AuxCoord(
            np.array(stat_points, dtype=(str, 20)),
            long_name='statistic',
            var_name=stats_var_name,
            units=cf_units.Unit('no_unit'))
        stats_longname_var_name = 'statistic_long_name' + '_' + \
                                  self.attributes['short_name']
        stats_longname_coord = iris.coords.AuxCoord(
            np.array(stat_longname_points, dtype=(str, 50)),
            long_name='statistic_long_name',
            var_name=stats_longname_var_name,
            units=cf_units.Unit('no_unit'))
        stats_units_var_name = 'statistic_units' + '_' + \
                               self.attributes['short_name']
        stats_units_coord = iris.coords.AuxCoord(
            np.array(stat_units_points, dtype=(str, 20)),
            long_name='statistic_units',
            var_name=stats_units_var_name,
            units=cf_units.Unit('no_unit'))

        #Convert values to np.array and add time dimension
        #(of shape (1,) to values array = > shape (nstats=2,ntimes=1)
        values = np.array(values)[:, np.newaxis]

        #Set up attributes dictionary - use the dictionary set up in __init__
        #but with some keys removed as these are used differently
        attributes = {}
        for attr, value in self.attributes.items():
            if attr != 'units' and attr != 'long_name':
                attributes[attr] = value

        #Now can create the iris Cube.
        #The dim coords are time and istats
        #istat is added on the dim=0 dimension (istats_coord, 0)
        #time is added on the dim=1 dimension (self.tcoord, 1)
        #The other coords contain the statistic name, long_name and units.
        #- these are added on the dim=0 dimension, to match istat.
        self.statscube = iris.cube.Cube(
            np.array(values),
            long_name=self.attributes['long_name'],
            units=1,
            attributes=attributes,
            var_name=self.attributes['long_name'],
            dim_coords_and_dims=[(istats_coord, 0),
                                 (self.tcoord, 1)],
            aux_coords_and_dims=[(stats_coord, 0),
                                 (stats_longname_coord, 0),
                                 (stats_units_coord, 0)])
        return self.statscube

    def convert_cube_to_attr(self):
        """
        Convert self.statscube into the attributes required for self, including
        the dictionary for self.statsdict and setting self.label and
        self.phenomenon if not already set
        """

        if self.statscube is None:
            return self.statsdict

        for stat_slice in self.statscube.slices_over('statistic'):
            stat_name = stat_slice.coord('statistic').points[0]
            self.statsdict[stat_name] = float(stat_slice.data)
            stat_unit = stat_slice.coord('statistic_units').points[0]
            if 'units' not in self.statsdict:
                #Set units using first available statistic units that
                #are a possible unit, ie not 1 or %
                if stat_unit != '1' and stat_unit != '%':
                    self.statsdict['units'] = stat_unit

        #if units are still not set, then default to '1'
        if 'units' not in self.statsdict:
            self.statsdict['units'] = '1'

        if self.label is None:
            self.label = self.statscube.attributes['label']

        if self.phenomenon is None:
            self.phenomenon = self.statscube.attributes['short_name']

        return self.statsdict



    #--------------------------------------------------------------------------#
    # All statistics calculations

    def calcdiff(self):
        """ Difference between model and obs - useful intermediate value"""
        self.diff = self.model - self.obs

    def nonzero(self):
        """ Find non-zero elements to avoid dividing by zero"""
        self.nzobs = self.obs.nonzero()
        self.count_nzobs = np.count_nonzero(self.obs)
        self.nz = (self.obs+self.model).nonzero()
        self.count_nz = np.count_nonzero(self.obs+self.model)

    def correlation(self):
        """
        Calculate pearson r correlation coefficient
          * -1 <= r <= +1
          * 1 => +ve correlation (good!)
        """
        if self.n > 0:
            self.statsdict['correlation'] = scipy.stats.pearsonr(self.model,
                                                                 self.obs)[0]
        else:
            self.statsdict['correlation'] = self.statsdict['mdi']

        return self.statsdict['correlation']

    def bias(self):
        """
        Calculate Bias - mean bias
          * Same units as obs
          * 0 => No bias (good!)
          * Positive bias for wind dir => model wind is clockwise of obs
        """
        if self.diff is None:
            self.calcdiff()
        if self.n > 0:
            self.statsdict['bias'] = np.sum(self.diff)/self.n
        else:
            self.statsdict['bias'] = self.statsdict['mdi']

        return self.statsdict['bias']

    def nmb(self):
        """
        NMB - Normalised mean bias
          * -1 <= nmb <= 1.0
          * 0 => No bias (good!)
          * 'Acceptable' if -0.2 < nmb < 0.2
        """
        if self.nzobs is None:
            self.nonzero()
        if self.count_nzobs > 0:
            self.statsdict['nmb'] = np.sum(self.diff) \
                                     / np.sum(self.obs[self.nzobs])
        else:
            self.statsdict['nmb'] = self.statsdict['mdi']

        return self.statsdict['nmb']

    def mnmb(self):
        """
        MNMB - Modified Normalised Mean Bias**
          * -2 <= mnmb <- +2
          * 0 => No bias (good!)
        """
        if self.nzobs is None:
            self.nonzero()
        if self.count_nz > 0:
            self.statsdict['mnmb'] = (2. / self.n) * \
                             np.sum((self.model[self.nz] - self.obs[self.nz]) \
                                    /(self.model[self.nz] + self.obs[self.nz]))
        else:
            self.statsdict['mnmb'] = self.statsdict['mdi']

        return self.statsdict['mnmb']

    def mge(self):
        """
        MGE - Mean gross error
          * Same units as obs
          * 0 <= mge
          * 0 => no error (good!)
        """
        if self.diff is None:
            self.calcdiff()
        if self.n > 0:
            self.statsdict['mge'] = np.sum(np.fabs(self.diff)) / self.n
        else:
            self.statsdict['mge'] = self.statsdict['mdi']

        return self.statsdict['mge']

    def nmge(self):
        """
        NMGE - Normalised mean gross error
          * 0 <= nmge <= 1
          * 0 => no error (good!)
        """
        if self.diff is None:
            self.calcdiff()
        if self.nzobs is None:
            self.nonzero()
        if self.count_nzobs > 0:
            self.statsdict['nmge'] = np.sum(np.fabs(self.diff)) \
                                    / np.sum(self.obs[self.nzobs])
        else:
            self.statsdict['nmge'] = self.statsdict['mdi']

        return self.statsdict['nmge']

    def fge(self):
        """
        FGE - Fractional gross error
          * 0 <= fge <= 2
          * 0 => no error (good!)
        """
        if self.count_nz is None:
            self.nonzero()
        if self.count_nz > 0:
            self.statsdict['fge'] = (2. / self.n) * \
                                    np.sum(np.fabs((self.model[self.nz] \
                                                       -self.obs[self.nz]) \
                                                     / (self.model[self.nz] \
                                                       + self.obs[self.nz])))
        else:
            self.statsdict['fge'] = self.statsdict['mdi']

        return self.statsdict['fge']

    def rmse(self):
        """
        RMSE - Root mean square error
          * 0 <= rmse
          * 0 => no error (good!)
        """
        if self.diff is None:
            self.calcdiff()
        if self.n > 0:
            self.statsdict['rmse'] = np.sqrt(np.sum(np.power(self.diff, 2)) \
                                             / self.n)
        else:
            self.statsdict['rmse'] = self.statsdict['mdi']

        return self.statsdict['rmse']

    def fac2(self):
        """
        FAC2 - Fraction of model within factor 2 of observations
          * 0 <= fac2 <= 1.0
          * 'Acceptable' for air quality if fac2 > 0.5
        """
        if self.count_nzobs is None:
            self.nonzero()
        if self.count_nzobs > 0:
            factor = self.model[self.nzobs] / self.obs[self.nzobs]
            fac2_indices = np.where((factor >= 0.5) & (factor <= 2.0))
            self.statsdict['fac2'] = float(np.size(fac2_indices)) / self.n
        else:
            self.statsdict['fac2'] = self.statsdict['mdi']

        return self.statsdict['fac2']

    def ioa(self):
        """
        IOA - Index of Agreement**
          * -1 <= IOA <= +1
          * +1 => good model performance
          * As defined in Willmott, C. J., Robeson, S. M., Matsuura, K., 2011.
              A refined index of model performance.
              International Journal of Climatology. 7, 64 :
              Interpretation of IOA is relatively straightforward. It indicates
              the sum of the magnitudes of the differences between the
              model-predicted and observed deviations about the observed mean
              relative to the sum of the magnitudes of the perfect-model
              (Mi = Oi, for all i) and observed deviations about the observed
              mean.
              A value of IOA of 0.5, for example, indicates that the sum of the
              error-magnitudes is one half of the sum of the
              perfect-model-deviation and observed-deviation magnitudes.
              When IOA = 0.0, it signifies that the sum of the magnitudes of the
              errors and the sum of the perfect-model-deviation and
              observed-deviation magnitudes are equivalent.
              When IOA = -0.5, it indicates that the sum of the error-magnitudes
              is twice the sum of the perfect-model-deviation and
              observed-deviation magnitudes.
              Values of IOA near -1.0 can mean that the model-estimated
              deviations about mean(obs) are poor estimates of the observed
              deviations; but, they also can mean that there simply is little
              observed variability. As the lower limit of IOA is approached,
              interpretations should be made cautiously.
         """
        if self.n > 0:
            ioa_component1 = np.sum(np.fabs(self.model - self.obs))
            ioa_component2 = 2. * np.sum(np.fabs(self.obs \
                                                     - np.mean(self.obs)))
            if ioa_component1 == 0:
                self.statsdict['ioa'] = 1.
            elif ioa_component1 <= ioa_component2:
                self.statsdict['ioa'] = 1. - (ioa_component1 / ioa_component2)
            else:
                self.statsdict['ioa'] = (ioa_component2 / ioa_component1) - 1.
        else:
            self.statsdict['ioa'] = self.statsdict['mdi']

        return self.statsdict['ioa']

    def contingency_table(self, threshold=None):
        """
        Contingency table:

        * +--------------+-------------+------------+
          |         Obs  | >=threshold | <threshold |
          +--------------+-------------+------------+
          | Model        |             |            |
          +--------------+-------------+------------+
          | >=threshold  |      a      |     b      |
          +--------------+-------------+------------+
          | < threshold  |      c      |     d      |
          +--------------+-------------+------------+
        * ct_a - Contigency table value a
               - 'o>=t_m>=t'
               - Number of times when obs ge threshold and model ge threshold
        * ct_b - Contigency table value b
               - 'o<t_m>=t'
               - Number of times when obs lt threshold and model ge threshold
        * ct_c - Contigency table value c
               - 'o>=t_m<t'
               - Number of times when obs ge threshold and model lt threshold
        * ct_d - Contigency table value d
               - 'o<t_m<t'
               - Number of times when obs lt threshold and model lt threshold
        * odds_ratio
            - 0 <= odds_ratio <= infinity
            - ratio of odds of a hit, to odds of a false alarm
            -  0 => no skill
            - >0 => Increasing skill as odds_ratio increases
        * orss - Odds ratio skills score
            - -1 <= ORSS <= +1
            - -1 => Strong -ve association with obs (bad!)
            -  0 => Random forecast
            - +1 => Strong +ve association with obs (good!)
        * Hitrate
            - probability that an event is forecast given that is observed.
            - 0 <= hitrate <= 1
            - 1 => event always forecast when observed (good!)
            - 0 => event never forecast when observed (bad!)
        * Falsealarmrate
            - probability that an event is forecast given that is was
              not observed.
            - 0 <= falsealarmrate <= 1
            - 1 => event is always forecast when it is not observed (bad!)
            - 0 => event is never forecast when it is not observed (good!)
        * Falsealarmratio
            - Fraction of the forecasts of the event associated with
              non-occurrences
        * Return value is (ct_a, ct_b, ct_c, ct_d)
        """
        if threshold is not None:
            self.statsdict['threshold'] = threshold
        #Set defaults
        self.statsdict['o>=t_m>=t'] = self.statsdict['mdi']
        self.statsdict['o<t_m>=t'] = self.statsdict['mdi']
        self.statsdict['o>=t_m<t'] = self.statsdict['mdi']
        self.statsdict['o<t_m<t'] = self.statsdict['mdi']
        self.statsdict['odds_ratio'] = self.statsdict['mdi']
        self.statsdict['orss'] = self.statsdict['mdi']
        self.statsdict['hitrate'] = self.statsdict['mdi']
        self.statsdict['falsealarmrate'] = self.statsdict['mdi']
        self.statsdict['falsealarmratio'] = self.statsdict['mdi']

        if self.statsdict['threshold'] is not None:
            threshold = self.statsdict['threshold']
            thresindex = np.where((self.obs >= threshold) &
                                  (self.model >= threshold))
            ct_a = np.size(thresindex)
            thresindex = np.where((self.obs < threshold) &
                                  (self.model >= threshold))
            ct_b = np.size(thresindex)
            thresindex = np.where((self.obs >= threshold) &
                                  (self.model < threshold))
            ct_c = np.size(thresindex)
            thresindex = np.where((self.obs < threshold) &
                                  (self.model < threshold))
            ct_d = np.size(thresindex)

            if ct_b*ct_c > 0:
                self.statsdict['odds_ratio'] = (float(ct_a)*float(ct_d)) \
                                           / (float(ct_b)*float(ct_c))
                self.statsdict['orss'] = (self.statsdict['odds_ratio']-1) \
                                           / (self.statsdict['odds_ratio']+1)

            if ct_a+ct_c > 0:
                self.statsdict['hitrate'] = float(ct_a) \
                                             / (float(ct_a)+float(ct_c))

            if ct_b+ct_d > 0:
                self.statsdict['falsealarmrate'] = float(ct_b) \
                                                   / (float(ct_b)+float(ct_d))

            if ct_a+ct_b > 0:
                self.statsdict['falsealarmratio'] = float(ct_b) \
                                                    / (float(ct_a)+float(ct_b))

            self.statsdict['o>=t_m>=t'] = ct_a
            self.statsdict['o<t_m>=t'] = ct_b
            self.statsdict['o>=t_m<t'] = ct_c
            self.statsdict['o<t_m<t'] = ct_d

        return self.statsdict['o>=t_m>=t'], self.statsdict['o<t_m>=t'], \
                self.statsdict['o>=t_m<t'], self.statsdict['o<t_m<t']


    def maxobs(self):
        """
        maxobs - max observed value
        """
        if self.n > 0:
            self.statsdict['maxobs'] = np.max(self.obs)
        else:
            self.statsdict['maxobs'] = self.statsdict['mdi']

        return self.statsdict['maxobs']

    def maxmod(self):
        """
        maxmod - max model value when obs avaliable
        """
        if self.n > 0:
            self.statsdict['maxmod'] = np.max(self.model)
        else:
            self.statsdict['maxmod'] = self.statsdict['mdi']

        return self.statsdict['maxmod']

    def meanobs(self):
        """
        meanobs - mean observed value
        """
        if self.n > 0:
            self.statsdict['meanobs'] = np.mean(self.obs)
        else:
            self.statsdict['meanobs'] = self.statsdict['mdi']

        return self.statsdict['meanobs']

    def meanmod(self):
        """
        meanmod - mean model value when obs avaliable
        """
        if self.n > 0:
            self.statsdict['meanmod'] = np.mean(self.model)
        else:
            self.statsdict['meanmod'] = self.statsdict['mdi']

        return self.statsdict['meanmod']

    def sdobs(self):
        """
        sdobs - standard deviation of observations
        """
        if self.n > 0:
            self.statsdict['sdobs'] = np.nanstd(self.obs, ddof=1)
        else:
            self.statsdict['sdobs'] = self.statsdict['mdi']

        return self.statsdict['sdobs']

    def sdmod(self):
        """
        sdmod - standard deviation of model data when obs avaliable
        """
        if self.n > 0:
            self.statsdict['sdmod'] = np.nanstd(self.model, ddof=1)
        else:
            self.statsdict['sdmod'] = self.statsdict['mdi']

        return self.statsdict['sdmod']

    def percentage(self):
        """
        * perc_correct - percentage of correct forecasts
            - best used for integer data (eg DAQI)
        * perc_under - percentage of under forecasts
            - best used for integer data (eg DAQI)
        * perc_over - percentage of over forecasts
            - best used for integer data (eg DAQI)

        returns (perc_correct, perc_under, perc_over)
        """
        if self.n > 0:
            ncorrect = 0
            nover = 0
            nunder = 0
            for mod, obs in zip(self.model, self.obs):
                if mod == obs:
                    ncorrect += 1
                elif mod > obs:
                    nover += 1
                elif mod < obs:
                    nunder += 1
            totaln = ncorrect+nover+nunder
            self.statsdict['perc_correct'] = 100. * (float(ncorrect) \
                                                     / float(totaln))
            self.statsdict['perc_over'] = 100. * (float(nover) \
                                                  / float(totaln))
            self.statsdict['perc_under'] = 100. * (float(nunder) \
                                                   / float(totaln))
        else:
            self.statsdict['perc_correct'] = self.statsdict['mdi']
            self.statsdict['perc_over'] = self.statsdict['mdi']
            self.statsdict['perc_under'] = self.statsdict['mdi']

        return self.statsdict['perc_correct'], self.statsdict['perc_over'], \
               self.statsdict['perc_under']

    def percentilemod(self, percentile_value):
        """
        User defined percentile of the model, defaults to 95th percentile
        """
        percentile_mod_name = 'pc' + str(percentile_value) + 'mod'
        percentile_value = int(percentile_value)

        if percentile_mod_name not in STATS_INFO:
            STATS_INFO[percentile_mod_name] = dict(STATS_INFO["percmod"])
            info = STATS_INFO["percmod"]["long_name"]
            name = info.format(n=percentile_value)
            STATS_INFO[percentile_mod_name]["long_name"] = name


        if self.n > 0:
            self.statsdict[percentile_mod_name] = np.percentile(self.model, percentile_value)
        else:
            self.statsdict[percentile_mod_name] = self.statsdict['mdi']

        return self.statsdict[percentile_mod_name]

    def percentileobs(self, percentile_value):
        """
        User defined percentile of the obs, defaults to 95th percentile
        """
        percentile_obs_name = 'pc' + str(percentile_value) + 'obs'
        percentile_value = int(percentile_value)

        if percentile_obs_name not in STATS_INFO:
            STATS_INFO[percentile_obs_name] = dict(STATS_INFO["percobs"])
            info = STATS_INFO["percobs"]["long_name"]
            name = info.format(n=percentile_value)
            STATS_INFO[percentile_obs_name]["long_name"] = name


        if self.n > 0:
            self.statsdict[percentile_obs_name] = np.percentile(self.obs, percentile_value)
        else:
            self.statsdict[percentile_obs_name] = self.statsdict['mdi']

        return self.statsdict[percentile_obs_name]


    #--------------------------------------------------------------------------#

def best_worst_indices(values, stat_name):
    """
    Calculate the indices of the 'best' and 'worst' values within the input
    list 'values' for given statistic.
    Note if two or more values are equivalently 'best' or 'worst' this indice
    is not set.

    :param values: input list of values
    :param stat_name: statistics name (that should match a key in STATS_INFO)

    :returns: (best_indice, worst_indice). If these can not be set - as statistic
              cannot sensibly be ranked, or because two or more values match the
              'best' or 'worst' value, then this these are set to None.

    >>> values = [0.2, -0.9, 0.01, 0.8, -0.5]
    >>> best_indice, worst_indice = best_worst_indices(values, 'correlation')
    >>> print(best_indice, worst_indice)
    3 1
    >>> print(values[best_indice], values[worst_indice])
    0.8 -0.9

    >>> values = [0.2, -0.9, 0.01, 0.8, -0.5, -0.9]
    >>> best_indice, worst_indice = best_worst_indices(values, 'bias')
    >>> print(best_indice, worst_indice)
    2 None
    >>> print(values[best_indice])
    0.01

    """

    best_indice = None
    worst_indice = None

    if len(values) > 1:

        ranking_method = STATS_INFO[stat_name].get('ranking_method', None)
        if ranking_method is not None:
            #Ensure all values are floats and in a numpy array
            npvalues = np.array([np.nan if str(v).strip() == 'None'
                                 else float(v) for v in values])

            #Use argsort - note this returns the index of the lowest value first
            if ranking_method == 'max':
                indices = np.argsort(npvalues)[::-1]
            elif ranking_method == 'min':
                indices = np.argsort(npvalues)
            elif ranking_method == 'absmax':
                indices = np.argsort(np.abs(npvalues))[::-1]
            elif ranking_method == 'absmin':
                indices = np.argsort(np.abs(npvalues))

            #Ignore nan (aka None) values
            finite_indices = indices[np.isfinite(npvalues[indices])]

            if len(finite_indices) > 1:

                #Set the indice of the best and worst value
                if values[finite_indices[0]] != values[finite_indices[1]]:
                    best_indice = finite_indices[0]
                else:
                    best_indice = None
                if values[finite_indices[-1]] != values[finite_indices[-2]]:
                    worst_indice = finite_indices[-1]
                else:
                    worst_indice = None
                if best_indice == worst_indice:
                    best_indice = None
                    worst_indice = None


    return (best_indice, worst_indice)


def _get_stat_name(species, stat_name, statsdict_all, stat_abbrev=False):
    """
    Extract stat name with units for stats_string* functions

    :param species: short species name
    :param stat_name: statistic abbreviation name
    :param statsdict_all: dictionary that stores statistics for each species
                          and each model.
    :param stat_abbrev: logical value to choose between long or abbreviation
                        stat names.

    :returns: stat_name_unit: string with the complete stat name with the unit.
    """

    # Get stat name with units
    if stat_abbrev:
        stat_name_unit = stat_name
    else:
        stat_name_unit = STATS_INFO[stat_name]['long_name']


    # Check if the current stat_name has units checking
    # its entry in the STATS_info dictionary
    units = STATS_INFO[stat_name].get('units', None)

    if units == 'units':
        models = list(statsdict_all[species].keys())
        units = statsdict_all[species][models[0]]['units']
        if units != '1':
            stat_name_unit += " (%s)" % units

    # if units contains some predefined units (e.g %), add it
    elif units is not None:
        stat_name_unit += " (%s)" % str(units)

    return stat_name_unit

def _get_stat_value(statsdict_species, stat_name, model):
    """
    Function to check the stats dictionary, set up the NAN value and extract the
    stat value in all stats_string* functions.

    :param statsdict_species: a dictionary with the stat for all model for one
                              species.
    :param stat_name: string name of the statistic.
    :param model: string containing the current model for the current species.

    :returns: value: a string containing the stat value.
    """

    value = None

    # Firstly check if the model is included in the
    #  stats dictionary for the current species, if not,
    #  stat is equal to None or preloaded NAN
    if model in statsdict_species:
        statsdict = statsdict_species[model]

        # Check if the current stat_name is inside the
        #  stats dictionary of the current species;
        #  if not, stat is equal to None
        if stat_name in statsdict:
            # check the type of stats:
            # add as float, or as string (e.g. when stat_name=units)
            value = statsdict[stat_name]
            if isinstance(value, float):

                # Check if stat value is under the limit, if not
                # write in scientific notation formatting
                if (-999.0 <= value <= -0.01) or (0.01 <= value <= 999.0) \
                    or (value == 0.0):
                    value = "%8.2f" % value
                else:
                    value = "%8.2e" % value

            else:
                value = "%8s" % str(value)
        else:
            # Check for preloaded Not A Number values
            #  in the stats dictionary
            if 'nanval' in statsdict:
                nanval = str(statsdict['nanval'])
            else:
                nanval = 'None'

            value = nanval

    return value


def _stats_string_csv(statsdict_all, allmodels, allspecies,
                      statsnames, header):
    """
    Given a sorted dictionary of species/model name dictionaries create an
    output string to write stats in a csv format.

    :param statsdict_all: dictionary that stores statistics for each species
                          and each model.
    :param allmodels: list containing all models evaluated for all species.
    :param allspecies: list containing all species evaluated.
    :param statsnames: list containing all statistics abbreviation.
    :param header: string passed in to be used as a header.

    :returns: output_string: string with all statistics, to be written
                             in csv format.

    Example:

    >>> import config
    >>> statsfile = config.CODE_DIR + 'adaqdocs/figures/stats.csv'
    >>> with open(statsfile, "r") as fin:
    ...     for iline, line in enumerate(fin):
    ...         print(line.strip())
    ...         if iline > 10: break
    Statistics for Tue 25-Mar-2014 23:00 - Thu 27-Mar-2014 00:00,,,
    Phenomemon,           Statistic,   model1,   model2
    <BLANKLINE>
    O3,              nsites,      None,      None
    O3,                npts,         5,         4
    O3,         correlation,     -0.21,     -0.62
    O3,        bias (ug/m3),     -1.31,      1.98
    O3,                 nmb,     -0.01,      0.02
    O3,                mnmb, -8.47e-03,      0.02
    O3,                nmge,      0.23,      0.20
    O3,                 fge,      0.22,      0.19
    O3,        rmse (ug/m3),     30.02,     26.46
    """

    statsdict_species = {}
    #Generate output strings
    output_string = ''
    if header is not None:
        # Remove any commas in the header
        header = header.replace(',', '')
        # Add the commas at the end of the header, to ensure
        # same number of fields in every line in output file.
        commas = (len(allmodels) + 1) * ','
        output_string += header + commas + '\n'

    output_string += "%15s," % "Phenomemon"
    output_string += "%20s" % "Statistic"
    output_string += len(allmodels)*",%9s" % tuple(allmodels)
    output_string += '\n'

    # This is a three nested loops to access to the nested dictionaries in
    #  statsdict_all that contains all the statistics for each species
    #  and each model.
    # To access to one statistic we previously need the species key, and the
    #  model key, e.g.:  statsdict_all[species][model][stat_name]
    # The final output_string contains stats sorted by species, statistics
    #  and model in this order.

    # Loop through species keys
    for species in allspecies:
        output_string += '\n'

        # Loop through stat_name
        for stat_name in statsnames:
            output_string += "%15s," % species

            # Get stat name with units
            stat_name_unit = _get_stat_name(species, stat_name, statsdict_all,
                                            stat_abbrev=True)
            output_string += "%20s" % stat_name_unit

            # loop through all models
            for model in allmodels:
                statsdict_species = statsdict_all[species]

                # Obtain stats value, checking for existence in the
                # stats dictionary or writing NAN value if not
                value = _get_stat_value(statsdict_species, stat_name, model)

                output_string += ", %9s" % value

            output_string += '\n'


    return output_string

def _stats_string_wiki(statsdict_all, allmodels, allspecies,
                       statsnames, header, colour=False):
    """
    Given a sorted dictionary of species/model name dictionaries create an
    output string to write stats in trac's wiki format.

    :param statsdict_all: dictionary that stores statistics for each species
                          and each model.
    :param allmodels: list containing all models evaluated for all species.
    :param allspecies: list containing all species evaluated.
    :param statsnames: list containing all statistics abbreviation.
    :param header: string passed in to be used as a header.
    :param colour: logical, if True, then highlights best model in green and
                   worst in red for each statistic.

    :returns: output_string: string with all statistics, to be written in
                             Trac's wiki format.

    Example:

    >>> import config
    >>> statsfile = config.CODE_DIR + 'adaqdocs/figures/stats.wiki-colour'
    >>> with open(statsfile, "r") as fin:
    ...     for iline, line in enumerate(fin):
    ...         print(line.strip())
    ...         if iline > 10: break # doctest: +NORMALIZE_WHITESPACE
    === Statistics for Tue 25-Mar-2014 23:00 - Thu 27-Mar-2014 00:00 ===
    ||  ||= '''model1''' =||= '''model2''' =||
    ||||||= '''O3''' =||
    ||=Number of sites                =|| \
None||                                   None||
    ||=Number of points               =|| \
[[span(style=color: green,        5)]]|| [[span(style=color: red  ,        4)]]||
    ||=Correlation                    =|| \
[[span(style=color: green,    -0.21)]]|| [[span(style=color: red  ,    -0.62)]]||
    ||=Bias (ug/m3)                   =|| \
[[span(style=color: green,    -1.31)]]|| [[span(style=color: red  ,     1.98)]]||
    ||=Normalised Mean Bias           =|| \
[[span(style=color: green,    -0.01)]]|| [[span(style=color: red  ,     0.02)]]||
    ||=Modified Normalised Mean Bias  =|| \
[[span(style=color: green,-8.47e-03)]]|| [[span(style=color: red  ,     0.02)]]||
    ||=Normalised Mean Gross Error    =|| \
[[span(style=color: red  ,     0.23)]]|| [[span(style=color: green,     0.20)]]||
    ||=Fractional Gross Error         =|| \
[[span(style=color: red  ,     0.22)]]|| [[span(style=color: green,     0.19)]]||
    ||=Root Mean Square Error (ug/m3) =|| \
    [[span(style=color: red  ,    30.02)]]|| [[span(style=color: green,    26.46)]]||
    """

    statsdict_species = {}
    #Generate output strings
    output_string = ''
    if header is not None:
        output_string += "=== {:s} ===\n".format(header)
    output_string += '||  ||'
    output_string += len(allmodels)*"= '''%s''' =||" % tuple(allmodels)
    output_string += '\n'

    # This is a three nested loops to access to the nested dictionaries in
    #  statsdict_all that contains all the statistics for each species
    #  and each model.
    # To access to one statistic we previously need the species key, and the
    #  model key, e.g.:  statsdict_all[species][model][stat_name]
    # The final output_string contains stats sorted by species, statistics
    #  and model in this order.

    # Loop through species keys
    for species in allspecies:
        output_string += '||'*len(allmodels)
        output_string += "||= '''%s''' =||" % species
        output_string += '\n'

        # Loop through stat_name
        for stat_name in statsnames:

            # Get stat name and unit
            stat_name_unit = _get_stat_name(species, stat_name, statsdict_all)

            output_string += "||=%-46s =||" % stat_name_unit
            values = []

            # loop through all models
            for model in allmodels:
                statsdict_species = statsdict_all[species]

                # Obtain stats value, checking for existence in the
                # stats dictionary or writing NAN value if not
                value = _get_stat_value(statsdict_species, stat_name, model)
                values.append(value)

                if not colour:
                    output_string += " %9s||" % value

            if colour:
                #Assign red or green colours depending on the best and worst values
                best_indice, worst_indice = best_worst_indices(values, stat_name)
                for i, value in enumerate(values):
                    if i == best_indice:
                        output_string += " [[span(style=color: green,%9s)]]||" % value
                    elif i == worst_indice:
                        output_string += " [[span(style=color: red  ,%9s)]]||" % value
                    else:
                        output_string += " %38s||" % value

            output_string += '\n'


    return output_string

def _stats_string_html(statsdict_all, allmodels, allspecies,
                       statsnames, header):
    """
    Given a sorted dictionary of species/model name dictionaries create an
    output string to write stats in html format.

    :param statsdict_all: dictionary that stores statistics for each species
                          and each model.
    :param allmodels: list containing all models evaluated for all species.
    :param allspecies: list containing all species evaluated.
    :param statsnames: list containing all statistics abbreviation.
    :param header: string passed in to be used as a header..

    :returns: output_string: string with all statistics, to be written
                             in html format.

    Example:

    >>> import config
    >>> statsfile = config.CODE_DIR + 'adaqdocs/figures/stats.html'
    >>> with open(statsfile, "r") as fin:
    ...     print('reading', statsfile)
    ...     for iline, line in enumerate(fin):
    ...         print(line.strip())
    ...         if iline > 19: break  # doctest: +ELLIPSIS
    reading .../stats.html
    <h3 style="text-align:center;">Statistics for Tue 25-Mar-2014 23:00 - Thu \
27-Mar-2014 00:00</h>
    <br>
    <table style="width: 50%"; border="1"; table align="center";>
    <tr>
    <th>  </th>
    <th align="center";>    model1 </th>
    <th align="center";>    model2 </th>
    </tr>
    <tr>
    <th colspan="3"; align="center";>         O3 </th>
    </tr>
    <tr>
    <td align="left"> Number of sites </td>
    <td  align="right">     None </td>
    <td  align="right">     None </td>
    </tr>
    <tr>
    <td align="left"> Number of points </td>
    <td  align="right">        5 </td>
    <td  align="right">        4 </td>
    </tr>
    """

    statsdict_species = {}
    output_string = ''
    # Include header as a html5 level 3 header
    if header is not None:
        string = '<h3 style="text-align:center;">{}</h>\n'
        output_string += string.format(header)
        output_string += '<br>\n'
    #Generate output strings
    string = '<table style="width: 50%"; border="1"; table align="center";>\n'
    output_string += string
    output_string += '<tr> \n'
    output_string += '  <th>  </th>\n'
    output_string += (len(allmodels)*'  <th align="center";> %9s </th>\n' %
                      tuple(allmodels))
    output_string += '</tr>\n'

    # This is a three nested loops to access to the nested dictionaries in
    #  statsdict_all that contains all the statistics for each species
    #  and each model.
    # To access to one statistic we previously need the species key, and the
    #  model key, e.g.:  statsdict_all[species][model][stat_name]
    # The final output_string contains stats sorted by species, statistics
    #  and model in this order.

    # Loop through species keys
    for species in allspecies:
        output_string += '<tr> \n'
        string = '  <th colspan="%s"; align="center";> %10s </th>\n'
        output_string += string % (len(allmodels) + 1, species)
        output_string += '</tr>\n'

        # Loop through stat_name
        for stat_name in statsnames:
            # Get stat name and unit
            stat_name_unit = _get_stat_name(species, stat_name, statsdict_all)

            output_string += '<tr>\n'
            output_string += '<td align="left"> %s </td>\n' % stat_name_unit

            # loop through all models
            for model in allmodels:
                statsdict_species = statsdict_all[species]

                # Obtain stats value, checking for existence in the
                # stats dictionary or writing NAN value if not
                value = _get_stat_value(statsdict_species, stat_name, model)

                output_string += ("  <td  align=\"right\"> %8s </td>\n"% value)
            output_string += '</tr>\n'

    output_string += '</table>'


    return output_string


def _stats_string_latex(statsdict_all, allmodels, allspecies,
                        statsnames, header):
    #Don't let pylint complain about ' Anomalous backslash in string: '\ '.
    #String constant might be missing an r prefix.
    #(anomalous-backslash-in-string)
    #pylint: disable=W1401
    """
    Given a sorted dictionary of species/model name dictionaries create an
    output string to write stats in latex format (one table for each species).

    :param statsdict_all: dictionary that stores statistics for each species
                          and each model.
    :param allmodels: list containing all models evaluated for all species.
    :param allspecies: list containing all species evaluated.
    :param statsnames: list containing all statistics abbreviation.
    :param header: string passed in to be used as a header.

    :returns: output_string: string with all statistics, to be written
                             in latex.

    Example:

    >>> import config
    >>> statsfile = config.CODE_DIR + 'adaqdocs/figures/stats.tex'
    >>> with open(statsfile, "r") as fin:
    ...     for iline, line in enumerate(fin):
    ...         print(line.strip())
    ...         if iline > 10: break
    \\begin{table}
    \\begin{center}
    \\begin{tabular}{|l|c|c|} \hline
    & \\textbf{model1}  & \\textbf{model2} \\\ \hline
    Number of sites                               &      None &      None \
\\\  \hline
    Number of points                              &         5 &         4 \
\\\  \hline
    Correlation                                   &     -0.21 &     -0.62 \
\\\  \hline
    Bias ($\mu g\ m^{-3}$)                        &     -1.31 &      1.98 \
\\\  \hline
    Normalised Mean Bias                          &     -0.01 &      0.02 \
\\\  \hline
    Modified Normalised Mean Bias                 & -8.47e-03 &      0.02 \
\\\  \hline
    Normalised Mean Gross Error                   &      0.23 &      0.20 \
\\\  \hline
    Fractional Gross Error                        &      0.22 &      0.19 \
\\\  \hline
    """

    statsdict_species = {}
    #Generate output strings
    output_string = ''

    # This is a three nested loops to access to the nested dictionaries in
    #  statsdict_all that contains all the statistics for each species
    #  and each model.
    # To access to one statistic we previously need the species key, and
    #  the model key, e.g.:  statsdict_all[species][model][stat_name]
    # The final output_string contains stats sorted by species, statistics
    #  and model in this order.

    # Loop through species keys
    for species in allspecies:

        output_string += '\\begin{table}\n'
        output_string += '  \\begin{center}\n'
        output_string += ('    \\begin{tabular}{|l|' + \
                          'c|'*len(allmodels) + '} \\hline\n')
        output_string += len(allmodels)*" & \\textbf{%s} " % tuple(allmodels)
        output_string += '\\\\ \\hline \n'

        # Loop through stat_name
        for stat_name in statsnames:

            # Get stat name and unit
            stat_name_unit = _get_stat_name(species, stat_name, statsdict_all)

            # Change any latex's not allowed character before add stat name
            # to the output string.
            for i, j in plotting_functions.LATEX.items():
                stat_name_unit = stat_name_unit.replace(i, j)
            output_string += (" %-45s " % stat_name_unit)

            # loop through all models
            for model in allmodels:
                statsdict_species = statsdict_all[species]

                # Obtain stats value, checking for existence in the
                # stats dictionary or writing NAN value if not
                value = _get_stat_value(statsdict_species, stat_name, model)
                if value is not None:
                    for i, j in plotting_functions.LATEX.items():
                        value = value.replace(i, j)
                    output_string += "& %9s " % value
                else:
                    output_string += "& %9s " % value

            output_string += '\\\\  \\hline\n'
        output_string += '    \\end{tabular}\n'

        # Write out latex table caption using header.
        # Firstly filter incompatible latex characters in species names
        for i, j in plotting_functions.LATEX.items():
            species = species.replace(i, j)

        if header is not None:
            output_string += ('    \\caption{%s: %s}\n' % (species, header))

        output_string += '    \\label{Stats%s}\n' % species
        output_string += '  \\end{center}\n'
        output_string += '\\end{table}\n\n'

    return output_string
    #pylint: enable=W1401

def stats_to_nc(tsstats_list, filename_prefix, nc_append=False):
    """
    Saves to netcdf format.

    :param tsstats_list: list of :class:`timeseries_stats.TimeSeriesStats`
                         objects.
    :param filename_prefix: name of the output file without path or file
                            extension. If set to None, the default value is
                            stats without the directory name.
    :param nc_append: If writing to netcdf file, then if True, append to
                      file of the same name if it already exists. If False,
                      then overwrites the file instead.

    Method:
        This is essentially saving to a netcdf file.
        However if nc_append, then reads in the pre-existing file
        and does a concatenate with any matching cubes.
        However extra code is required to future proof this code against:

          * If the list of statistics coordinates changes
          * If the list of short_names saved in the file changes.
          * If any times are missing then these are added in with nan data
            values.
    """

    #Set up dictionary, whose keys will be filenames, and whose
    #values will be a list of cubes which should be written to that
    #filename.
    file_cubelists = {}

    for stats in tsstats_list:
        #Generate filename using the cube label
        if stats.statscube is None:
            warnings.warn('Cannot output to netcdf file: no statistics cube')
            continue
        model = stats.statscube.attributes['label']
        obs = stats.statscube.attributes['obs']
        label = 'Obs-' + obs + '_' + 'Mod-' + model
        label = label.replace(' ', '')
        label = label.replace(',', '_')
        filename = filename_prefix + '_' + label + '.nc'

        if filename not in file_cubelists:
            file_cubelists[filename] = iris.cube.CubeList()

        cube = stats.statscube

        if nc_append and os.path.isfile(filename):
            #First read file in containing previously saved data
            prev_cubes = iris.load(
                filename,
                iris.AttributeConstraint(
                    short_name=cube.attributes['short_name']))
            if len(prev_cubes) == 1:
                prev_cube = prev_cubes[0]

                if prev_cube.coord('time') == cube.coord('time'):
                    #Don't use the previous cube if they have the same time in.
                    print('Overwriting existing data - same time coordinate.')

                else:

                    #If the prev_cube includes the time trying to be added
                    #(NB should only be one time in the current cube, but
                    # add this check in case not)
                    #then this time needs to be removed from the prev_cube
                    #so it can be overwritten by the current cube.
                    if len(cube.coord('time').points) == 1 and \
                       cube.coord('time').points[0] in \
                       prev_cube.coord('time').points:

                        indices = np.where(prev_cube.coord('time').points != \
                                           cube.coord('time').points[0])
                        #Extract just the indices from the time dimension
                        #(2nd dim)
                        assert prev_cube.coord_dims('time') == (1,)
                        prev_cube = prev_cube[:, indices[0]]

                    #Load data into memory now.
                    #Otherwise when trying to save back to same filename...
                    #data is only accessed from file as it tries to save, but
                    #by that point the file has been already destroyed.
                    #Pylint ignore "Statement has no effect error"
                    #pylint: disable=W0104
                    prev_cube.data
                    #pylint: enable=W0104

                    #Remove Conventions attribute as this is only added on save
                    #to netcdf so won't match the current cube
                    if 'Conventions' in prev_cube.attributes:
                        del prev_cube.attributes['Conventions']

                    #Ensure fill values are set the same for both, otherwise
                    #wont merge
                    ma.set_fill_value(prev_cube.data, 1e20)
                    ma.set_fill_value(cube.data, 1e20)

                    #Get list of integer statistic numbers for both cubes and
                    #the overall list between them.
                    prev_istats = prev_cube.coord('istatistic').points
                    cube_istats = cube.coord('istatistic').points
                    istats = np.unique(list(prev_istats) + list(cube_istats))

                    #Check if the istats from each cube are the same. If they
                    #are not, then need to add the extra statistics in.
                    add_extra_stats = False
                    if len(cube_istats) != len(prev_istats):
                        add_extra_stats = True
                    elif not (cube_istats == prev_istats).all():
                        add_extra_stats = True

                    if add_extra_stats:
                        #Develop a cubelist containing slices for each cube,
                        #with each slice having just a single statistic.
                        prev_cube_slices = iris.cube.CubeList()
                        for prev_cube_slice in \
                            prev_cube.slices_over('istatistic'):
                            prev_cube_slices.append(prev_cube_slice)
                        cube_slices = iris.cube.CubeList()
                        for cube_slice in cube.slices_over('istatistic'):
                            cube_slices.append(cube_slice)


                        for istat in istats:
                            #Find any statistics which don't already exist in
                            #each cube
                            if istat not in prev_istats:
                                #Need to statistic add to prev_cube...
                                #Extract the equivalent slice from the other
                                #cube
                                cube_slice = cube.extract(
                                    iris.Constraint(istatistic=istat))
                                #Extract a sample slice from the current cube
                                prev_cube_slice = prev_cube_slices[0].copy()
                                #Modify the coordinates to take their points
                                #from the other cube (apart from time)
                                for coord in prev_cube_slice.coords():
                                    coord_name = coord.name()
                                    if coord_name != 'time':
                                        coord.points = cube_slice.coords(
                                            coord_name)[0].points
                                #Overwrite data with nans
                                prev_cube_slice.data[:] = np.nan
                                #Add to list of slices
                                prev_cube_slices.append(prev_cube_slice)
                            if istat not in cube_istats:
                                #Need to add to cube
                                #Extract the equivalent slice from the other
                                #cube
                                prev_slice = prev_cube.extract(
                                    iris.Constraint(istatistic=istat))
                                #Extract a sample slice from the current cube
                                cube_slice = cube_slices[0].copy()
                                #Modify the coordinates to take their points
                                #from the other cube (apart from time)
                                for coord in cube_slice.coords():
                                    coord_name = coord.name()
                                    if coord_name != 'time':
                                        coord.points = prev_slice.coords(
                                            coord_name)[0].points
                                #Overwrite data with nans
                                cube_slice.data[:] = np.nan
                                #Add to list of slices
                                cube_slices.append(cube_slice)

                        #Now merge all the slices back into their original cube,
                        #but should now contain the extra statistics
                        prev_cube = prev_cube_slices.merge_cube()
                        cube = cube_slices.merge_cube()

                    #Now finally the current cube and the previous cube should
                    #match sufficiently that they should be able to concatenate
                    #their time-axes together.
                    try:
                        #Try concatenating first as this is much quicker than
                        #bit under except...
                        cube = iris.cube.CubeList(
                            [prev_cube, cube]).concatenate_cube()
                    except iris.exceptions.ConcatenateError:
                        #Usual reason for not being able to concatenate
                        #is that the prev_cube contains times both before and
                        #after the time in cube, which means that they can not
                        #be directly appended to each other. Instead:
                        #Split each cube up into multiple cubes, each of which
                        #only contain a single time
                        cube_list = iris.cube.CubeList(
                            [c for c in cube.slices_over('time')])
                        prev_cube_list = iris.cube.CubeList(
                            [c for c in prev_cube.slices_over('time')])
                        #Now merge the cubes back in the right order
                        cube = (prev_cube_list + cube_list).merge_cube()
                        #Now check and make sure that istatistic is first coord,
                        #and time is second coord - fix the other way around
                        if cube.coord_dims('istatistic') == (1,) and \
                           cube.coord_dims('time') == (0,):
                            #Fix as not as expected
                            cube.transpose([1, 0])

                    #Add missing times in so does not have to be done in code
                    #everytime otherwise
                    if len(cube.coord('time').points) > 2:
                        cube = cube_time.cube_add_missing_times(cube)



            elif not prev_cubes:
                print('No matching cubes found')
            else:
                print('Multiple cubes loaded: ', prev_cubes)

        file_cubelists[filename].append(cube)

    if nc_append:
        #Check if any other cubes in the files which might otherwise be lost
        for filename, cubelist in file_cubelists.items():
            if not os.path.isfile(filename):
                #No file to read in
                continue
            file_cubelist = iris.load(filename)
            #Get list of shortnames available in file
            file_shortnames = [cube.attributes['short_name']
                               for cube in file_cubelist]
            #Get list of shortnames already in use
            cubelist_shortnames = [cube.attributes['short_name']
                                   for cube in cubelist]
            for short_name in file_shortnames:
                if short_name not in cubelist_shortnames:
                    #Also add this cube to cubelist
                    cube = file_cubelist.extract(iris.AttributeConstraint(
                        short_name=short_name), strict=True)

                    #Access data and coords so can overwrite file
                    #Pylint ignore "Statement has no effect error"
                    #Pylint ignore "Instance of 'CubeList' has no 'coords'
                    #member (no-member)" - using strict on the extract above
                    #should guaruntee that a single cube, rather than a cubelist
                    #is returned.
                    #pylint: disable=W0104, E1101
                    for coord in cube.coords():
                        coord.points
                    cube.data
                    #pylint: enable=W0104, E1101
                    #Add to cubelist
                    cubelist.append(cube)


    #Save the cubes
    #(ensure always output in sorted order - cannot guarantee order if just
    # using iteritems on dictionary which can confuse doctests otherwise)
    for filename, cubelist in sorted(file_cubelists.items()):
        iris.save(cubelist, filename)
        print('Saved to ', filename)




def save_stats(tsstats_list, percentile_value=95, filename_prefix=None,
               directory=None, format_list=None,
               header=None, nc_append=False):
    '''
    Given a list of tsstats objects, saves information to a fixed-format,
    e.g. csv file of the format:

    .. literalinclude:: ../adaqdocs/figures/stats.csv

    :param tsstats_list: list of :class:`timeseries_stats.TimeSeriesStats`
     objects.
    :param percentile_value: The percentile number requested by the user
    :param filename_prefix: name of the output file without path or file
                            extension. If set to None, the default value is
                            stats without the directory name.
    :param directory: path pointing to the directory to save the stats file.
    :param format_list: list of file formats (csv, wiki, html, tex, nc).
    :param header: string passed in to be used as a header in ascii files.
    :param nc_append: If writing to netcdf file, then if True, append to
                      file of the same name if it already exists. If False,
                      then overwrites the file instead.

    Example of saving statistics, by making use of label and
    phenomenon keywords in TimeSeriesStats.
    Firstly, set up three sets of TimeSeriesStats objects.

    >>> obslist    = [None,    99.80, 133.56, 151.24, 107.61, 99.9 ]
    >>> modellist1 = [131.12, 131.65,  94.16, 121.12, 139.73, 98.9 ]
    >>> modellist2 = [130.,     132.,  94.16,   None, 122.13, 100.5]
    >>> stats1 = TimeSeriesStats(obslist,modellist1,mdi=None,threshold=100.,
    ... label='model1',phenomenon='O3')
    >>> statsdict1 = stats1.get_basic_stats()
    >>> statsdict1['units'] = 'ug/m3'
    >>> stats2 = TimeSeriesStats(obslist,modellist2,mdi=None,threshold=100.,
    ... label='model2',phenomenon='O3')
    >>> statsdict2 = stats2.get_basic_stats()
    >>> statsdict2['units'] = 'ug/m3'
    >>> stats3 = TimeSeriesStats(obslist,modellist2,mdi=None,threshold=100.,
    ... label='model2',phenomenon='CO')
    >>> statsdict3 = stats3.get_basic_stats()
    >>> statsdict3['units'] = 'ug/m3'

    Can now pass these objects as a list into save_stats

    >>> import config
    >>> header ='Statistics for Tue 25-Mar-2014 23:00 - Thu 27-Mar-2014 00:00'
    >>> save_stats([stats1,stats2,stats3],
    ...            filename_prefix='stats',
    ...            directory=config.CODE_DIR+'/adaqdocs/figures/',
    ...            format_list=['csv', 'wiki-colour', 'html', 'tex'],
    ...            header = header) # doctest: +ELLIPSIS
    Statistics saved to  .../stats.csv
    Statistics saved to  .../stats.wiki-colour
    Statistics saved to  .../stats.html
    Statistics saved to  .../stats.tex

    In order to produce netcdf output, the input data to TimeSeriesStats must
    have been cubes:

    >>> import adaq_functions
    >>> ini_dict, sites_data, od, md_list = adaq_functions.get_exampledata(
    ... exampletype='full') # doctest: +ELLIPSIS
    Reading inifile .../example_data_5days.ini
    Number of sites:  5
    >>> obs_cube = od.extract(short_name='O3', singlecube=True)
    >>> mod_cube = md_list[0].extract(short_name='O3', singlecube=True)
    >>> mod_cube2 = md_list[1].extract(short_name='O3', singlecube=True)
    >>> stats = TimeSeriesStats(obs_cube, mod_cube)
    >>> bias = stats.bias()
    >>> cube = stats.convert_to_cube()
    >>> save_stats([stats],
    ...           filename_prefix='stats',
    ...           directory=config.TEST_DIR,
    ...           format_list=['nc']) # doctest: +ELLIPSIS
    Saved to  .../stats_Obs-AURN_Mod-aqum_oper.nc

    Note that if there were multiple models then these would be saved in
    different netcdf files. However multiple short_names for the same model
    are saved in the same file.

    Using the nc_append attribute, running this code with some later data, the
    output netcdf file can be appended to:

    >>> ini_dict, sites_data, od, md_list = adaq_functions.get_exampledata()
    ... # doctest: +ELLIPSIS
    Reading inifile .../example_data_1days.ini
    Number of sites:  5
    >>> obs_cube = od.extract(short_name='O3', singlecube=True)
    >>> mod_cube = md_list[0].extract(short_name='O3', singlecube=True)
    >>> stats = TimeSeriesStats(obs_cube, mod_cube)
    >>> bias = stats.bias()
    >>> cube = stats.convert_to_cube()
    >>> save_stats([stats],
    ...           filename_prefix='stats',
    ...           directory=config.TEST_DIR,
    ...           format_list=['nc'], nc_append=True) # doctest: +ELLIPSIS
    Saved to  .../stats_Obs-AURN_Mod-aqum_oper.nc

    Load the cube(s) back in to check:

    >>> cubes = iris.load(
    ... config.TEST_DIR+'/stats_Obs-AURN_Mod-aqum_oper.nc')
    >>> print(cubes[0])  # doctest: +ELLIPSIS
    mass_concentration_of_ozone_in_air / (1) (istatistic: 5; time: 2)
         Dimension coordinates:
              istatistic                                x        -
              time                                      -        x
         Auxiliary coordinates:
              statistic                                 x        -
              statistic_long_name                       x        -
              statistic_units                           x        -
         Attributes:
              ...
              label: aqum_oper
              obs: AURN
              short_name: O3

    '''

    if format_list is None:
        format_list = ['csv']

    if filename_prefix is None:
        filename_prefix = 'stats'

    if directory is not None:
        if not directory.endswith('/'):
            #Ensure directory has '/' on the end
            directory += '/'

    if filename_prefix is None:
        if directory is None:
            raise IOError('Filename and directory not fixed!')
        else:
            filename_prefix = directory + 'stats'
    else:
        filename_prefix = directory + filename_prefix

    #Check output directory exists and create if not
    if directory: #Checks not current directory
        if not os.path.isdir(directory):
            os.makedirs(directory)

    #Preferred order of stats
    statsnames = ['nsites', 'npts', 'correlation', 'bias', 'nmb', 'mnmb',
                  'nmge', 'fge', 'rmse', 'fac2', 'ioa', 'threshold',
                  'odds_ratio', 'orss', 'hitrate', 'falsealarmrate',
                  'falsealarmratio', 'o>=t_m>=t', 'o<t_m>=t', 'o>=t_m<t',
                  'o<t_m<t', 'maxobs', 'maxmod', 'pc' + str(percentile_value) + 'obs',
                  'pc' + str(percentile_value) + 'mod', 'meanobs', 'meanmod', 'sdobs',
                  'sdmod', 'perc_correct', 'perc_over', 'perc_under', 'units']

    if isinstance(tsstats_list, dict):
        #Convert to a list
        tsstats_list = list(tsstats_list)


    if 'nc' in format_list:
        #Write out to netcdf file
        stats_to_nc(tsstats_list, filename_prefix, nc_append)

    ascii = False
    for format_type in format_list:
        if format_type in ['csv', 'wiki', 'wiki-colour', 'html', 'tex']:
            ascii = True

    if ascii:
        #Sort out the dictionaries - put into species/model name dictionaries
        #Get list of species and model names
        allspecies = []
        allmodels = []
        statsdict_all = {}
        for idict, stats in enumerate(tsstats_list):
            species = stats.phenomenon
            if species is None:
                species = 'Unknown'
            if species not in allspecies:
                allspecies.append(species)
            if stats.label is not None:
                modelname = stats.label
            else:
                modelname = 'Model'+str(idict)
            if modelname not in allmodels:
                allmodels.append(modelname)
            if species not in statsdict_all:
                statsdict_all[species] = {}
            statsdict_all[species][modelname] = stats.statsdict


        #Generate output strings
        for format_type in format_list:
            filename = filename_prefix+'.'+format_type
            if format_type == 'csv':
                output_string = _stats_string_csv(statsdict_all, allmodels,
                                                  allspecies, statsnames,
                                                  header)
            elif format_type == 'wiki':
                output_string = _stats_string_wiki(statsdict_all, allmodels,
                                                   allspecies, statsnames,
                                                   header)
            elif format_type == 'wiki-colour':
                output_string = _stats_string_wiki(statsdict_all, allmodels,
                                                   allspecies, statsnames,
                                                   header, colour=True)
            elif format_type == 'html':
                output_string = _stats_string_html(statsdict_all, allmodels,
                                                   allspecies, statsnames,
                                                   header)
            elif format_type == 'tex':
                output_string = _stats_string_latex(statsdict_all, allmodels,
                                                    allspecies, statsnames,
                                                    header)
            else:
                #Unknown format, or netcdf
                continue

            #Print data to file
            with open(filename, "w") as fout:
                fout.write(output_string)

            print('Statistics saved to ', filename)


if __name__ == "__main__":

    import doctest
    doctest.testmod()
