"""
Contains statistical functions performed on cubes.
"""
from __future__ import print_function

from six.moves.builtins import str
import warnings
import datetime
from pyproj import Geod

import numpy as np
import iris
import iris.coord_categorisation as coord_cat

import array_statistics
import cube_time


#: Dictionary of possible aggregators to use on a cube.
#: Only common aggregators are currently included but this dictionary
#: could be expanded to include all aggregators in iris.analysis
CUBE_AGGREGATORS = {'NANMEAN' : iris.analysis.Aggregator('nanmean', np.nanmean),
                    'NANMAX'  : iris.analysis.Aggregator('nanmax', np.nanmax),
                    'NANMEDIAN' : iris.analysis.Aggregator('nanmedian', np.nanmedian),
                    'NANMIN'  : iris.analysis.Aggregator('nanmin', np.nanmin),
                    'NANSUM'  : iris.analysis.Aggregator('nansum', np.nansum),
                    'NANSTD'  : iris.analysis.Aggregator('nanstd', np.nanstd),
                    'NANVAR'  : iris.analysis.Aggregator('nanvar', np.nanvar),
                    'NANPERCENTILE' : iris.analysis.Aggregator('nanpercentile',
                                                               np.nanpercentile),
                    'MEAN'    : iris.analysis.MEAN,
                    'MAX'     : iris.analysis.MAX,
                    'MEDIAN'  : iris.analysis.MEDIAN,
                    'MIN'     : iris.analysis.MIN,
                    'COUNT'   : iris.analysis.COUNT,
                    'STD_DEV' : iris.analysis.STD_DEV,
                    'SUM'     : iris.analysis.SUM,
                    'VARIANCE': iris.analysis.VARIANCE,
                    'PERCENTILE': iris.analysis.PERCENTILE
                   }


def aggregate_time(cube, period='hour', aggregator_name='NANMEAN',
                   collapsed=False):
    """
    Aggregate cube according to a specific time component. The time
    co-ordinate is replaced with a co-ordinate for the given period.

    For example, when aggregating over 'hour', all points at the same
    time of day are combined, regardless of which day they are in.
    The resulting 'hour' co-ordinate therefore contains at most 24 points -
    one for each hour.

    .. note:: Not to be confused with :any:`periodic_stat`,
              which is used to (for example) aggregate all hours within
              the same day to a single point.

    :param cube: iris.cube.Cube
    :param period: period to aggregate over. Currently available options
                   are:

                    * 'hour' (gives a diurnal cube),
                    * 'monthly' (gives a monthly cube), also sets up extra coord
                      'monthname' containing three letter month name, eg 'Feb'.
                      Over three years, this will end up with a coord 12 points
                      long.
                    * 'yearmonth', aggregates into months but retains the year
                      component, so over three years will end up with a coord
                      3x12 points long. Also sets up extra coord 'yearmonthname'
                      containing three letter month name and year as a string,
                      eg 'Feb 2014'.

    :param aggregator_name: string name of aggregator to use. By default this
                             is 'NANMEAN' which will take a mean over all
                             values, ignoring NaNs. Other available alternatives
                             include 'NANMAX' (max, ignoring NaNs), 'MEAN' and
                             'MAX', which are the mean
                             and max values, taking NaNs into account (if any
                             value in array is NaN, returns NaN).
                             Uses cube_statistics.CUBE_AGGREGATORS dictionary.
                             Alternatively any aggregator in iris.analysis
                             could be used.
    :param collapsed: logical to determine whether to collapse over all
                      dimensions. If False (default), time dimension is removed
                      and replaced with hour, but all other dimensions kept.
                      If True, the same aggregator is used to collapse all other
                      dimensions. In this case,
                      the only dimension in the returned cube is time.


    Method:

    Works on a cube to return a cube aggregated by required period.
    When the cube is aggregated the time coordinate is no longer monotonic.
    The time coordinate and any other coordinates on this dimension
    (for example forecast period) are removed, and the new period coordinate is
    added instead as a dimcoord.
    To ensure the resulting cube can be plotted in the expected time order
    (eg 0-23Z for period='hour'), the aggregated cube is sliced by period
    and the indiviudal subcubes are then concatenated in the correct order.

    Example:

    >>> import config
    >>> import adaq_data
    >>> sample_data_path = config.SAMPLE_DATADIR+'sites_cube_list/'
    >>> md = adaq_data.ADAQData()
    >>> scl = md.load_ts(sample_data_path+'aqum_oper_5days.nc')

    >>> mod_cube = md.extract(short_name='O3', singlecube=True)
    >>> print(mod_cube) # doctest: +ELLIPSIS
    mass_concentration_of_ozone_in_air / (ug/m3) (site_id: 5; time: 121)
         Dimension coordinates:
              site_id                                    x        -
              time                                       -        x
         Auxiliary coordinates:
              abbrev                                     x        -
              grid_latitude                              x        -
              grid_longitude                             x        -
              latitude                                   x        -
              longitude                                  x        -
              site_altitude                              x        -
              site_name                                  x        -
              site_type                                  x        -
              surface_altitude                           x        -
              forecast_period                            -        x
         Scalar coordinates:
              forecast_day: 1.0 Days
              level_height: 20.00... m, bound=(0.0, 49.99...) m
              model_level_number: 1
              sigma: 0.99..., bound=(1.0, 0.99...)
         Attributes:
              Conventions: CF-1.5
              STASH: m01s34i001
              label: aqum_oper
              short_name: O3
              source: Data from Met Office Unified Model
         Cell methods:
              mean: time (1 hour)

    Calculate diurnal variation, but without collapsing the site_id coordinate,
    use period='hour':

    >>> diurnal_cube = aggregate_time(mod_cube, period='hour', collapsed=False)
    >>> print(diurnal_cube) # doctest: +ELLIPSIS
    mass_concentration_of_ozone_in_air / (ug/m3) (hour: 24; site_id: 5)
         Dimension coordinates:
              hour                                    x            -
              site_id                                 -            x
         Auxiliary coordinates:
              abbrev                                  -            x
              grid_latitude                           -            x
              grid_longitude                          -            x
              latitude                                -            x
              longitude                               -            x
              site_altitude                           -            x
              site_name                               -            x
              site_type                               -            x
              surface_altitude                        -            x
         Scalar coordinates:
              forecast_day: 1.0 Days
              level_height: 20.00... m, bound=(0.0, 49.99...) m
              model_level_number: 1
              sigma: 0.99..., bound=(1.0, 0.99...)
         Attributes:
              Conventions: CF-1.5
              STASH: m01s34i001
              label: aqum_oper
              short_name: O3
              source: Data from Met Office Unified Model
         Cell methods:
              mean: time (1 hour)
              nanmean: hour

    And now collapse the site_id coordinate
    (the nanmean is also taken across site_id):

    >>> diurnal_cube = aggregate_time(mod_cube, period='hour', collapsed=True)
    >>> print(diurnal_cube) # doctest: +ELLIPSIS
    mass_concentration_of_ozone_in_air / (ug/m3) (hour: 24)
         Dimension coordinates:
              hour                                    x
         Scalar coordinates:
              abbrev: YW|ACTH|AH|ABD|HAR
              forecast_day: 1.0 Days
              grid_latitude: 1.381... degrees, \
bound=(-1.895..., 4.658...) degrees
              grid_longitude: -0.021... degrees, \
bound=(-0.772..., 0.729...) degrees
              latitude: 53.877... degrees, \
bound=(50.597..., 57.157...) degrees
              level_height: 20.000... m, bound=(0.0, 49.998...) m
              longitude: -2.521... degrees, \
bound=(-3.716..., -1.326...) degrees
              model_level_number: 1
              sigma: 0.997..., bound=(1.0, 0.994...)
              site_altitude: 195... m, bound=(20, 370) m
              site_id: 35747847.141..., bound=(35628361.140..., 35867333.141...)
              site_name: Yarner_Wood|Auchencorth_Moss|Aston_Hill|\
Aberdeen|Harwell
              site_type: RURAL|RURAL|RURAL|URBAN_BACKGROUND|RURAL
              surface_altitude: 137.04... m, bound=(38.888..., 235.209...) m
         Attributes:
              Conventions: CF-1.5
              STASH: m01s34i001
              label: aqum_oper
              short_name: O3
              source: Data from Met Office Unified Model
         Cell methods:
              mean: time (1 hour)
              nanmean: hour
              nanmean: site_id

    This returned cube has only one dimension - hour, which is the required
    period.

    As an example of using period='month' to get a monthly cube:

    >>> scl = md.load_ts(sample_data_path+'aqum_oper_10days.nc')
    >>> mod_cube = md.extract(short_name='O3', singlecube=True)
    >>> cube = aggregate_time(mod_cube, period='month')
    >>> print(cube.coord('month').points)
    [3 4]
    >>> print(cube.coord('monthname').points)
    ['Mar' 'Apr']

    And finally using period='yearmonth':

    >>> cube = aggregate_time(mod_cube, period='yearmonth')
    >>> print(cube.coord('yearmonth').points)
    [201403 201404]
    >>> print(cube.coord('yearmonthname').points)
    ['Mar 2014' 'Apr 2014']

    """

    #If period already exists as a coordinate, remove first
    if cube.coords(period):
        warnings.warn('Removing pre-existing ' + period + ' coordinate')
        cube.remove_coord(period)

    timeunits = period
    if period == 'hour':
        #Add hour as a categorised co-ordinate
        coord_cat.add_categorised_coord(cube, 'hour',
                                        cube.coord('time'),
                                        cube_time.hour_from_time,
                                        units=timeunits)
    elif period == 'month':
        #Add month as a categorised co-ordinate
        coord_cat.add_categorised_coord(cube, 'month',
                                        cube.coord('time'),
                                        cube_time.month_from_time,
                                        units=timeunits)
    elif period == 'yearmonth':
        #Add month as a categorised co-ordinate
        timeunits = '1'
        coord_cat.add_categorised_coord(cube, 'yearmonth',
                                        cube.coord('time'),
                                        cube_time.yearmonth_from_time,
                                        units=timeunits)

    #Remove coordinates on original time-axis as these will
    #be incorrect/misleading and won't let it aggregate over hour otherwise.
    cube_tmp = cube.copy() #Don't modify input cube
    time_dim = cube_tmp.coord_dims('time')
    time_axis_coords = [coord.name() for coord
                        in cube_tmp.coords(dimensions=time_dim)]
    for coord_name in time_axis_coords:
        if coord_name != period:
            warnings.warn('Removing coordinate on time axis: ' + coord_name)
            cube_tmp.remove_coord(coord_name)

    aggregator = CUBE_AGGREGATORS[aggregator_name]

    #Aggregate cube by period using the aggregator
    aggregated_cube = cube_tmp.aggregated_by(period, aggregator)

    #Collapse cube across all non-time/period dim coords
    #to give a mean cube over all dimensions
    if collapsed:
        dim_coord_names = [coord.name() for coord
                           in aggregated_cube.dim_coords]
        for coord_name in dim_coord_names:
            if coord_name != 'time' and coord_name != period:
                aggregated_cube = aggregated_cube.collapsed(coord_name,
                                                            aggregator)

    # Generate cube in expected time order
    cube_list = iris.cube.CubeList()
    for subcube in aggregated_cube.slices_over(period):
        new_coord = subcube.coord(period)
        subcube.remove_coord(period)
        subcube.add_aux_coord(iris.coords.DimCoord(new_coord.points,
                                                   long_name=period,
                                                   units=timeunits))
        subcube = iris.util.new_axis(subcube, period)
        cube_list.append(subcube)


    new_cube = cube_list.concatenate_cube()


    if period == 'month':
        #Also add name of month as an extra coord
        #Eg 'Feb'
        coord_cat.add_categorised_coord(new_cube, 'monthname',
                                        new_cube.coord('month'),
                                        cube_time.monthname_from_month)
    if period == 'yearmonth':
        #Also add name of month and year as an extra coord
        #Eg 'Feb 2014'
        coord_cat.add_categorised_coord(new_cube, 'yearmonthname',
                                        new_cube.coord('yearmonth'),
                                        cube_time.yearmonthname_from_yearmonth)

        #Also add datetime of month and year as an extra coord
        #Convert to datetime format
        yearmonthdt = np.array(
            [datetime.datetime(int(str(point)[:4]), #year
                               int(str(point)[4:6]), #month
                               1) #first day of month
             for point in new_cube.coord('yearmonth').points])
        #Then convert this to numbers
        units = cube.coord('time').units
        #units = cf_units.Unit('hours since epoch', calendar='gregorian')
        yearmonth_num = np.array([units.date2num(t) for t in yearmonthdt])
        #Can now set up new coord
        coord = iris.coords.AuxCoord(yearmonth_num,
                                     long_name='yearmonthdt',
                                     units=units)
        #Guess bounds, where point is at the beginning of the month
        #so bound covers entire month
        if len(yearmonth_num) > 1:
            coord.guess_bounds(0.0)
        new_cube.add_aux_coord(coord, 0) #Add this coord to new_cube

    return new_cube


def periodic_stat(cube, stat='max', period='day', min_periods=1, aqdates=None):
    """
    Calculate daily, monthly, or yearly statistics.

    .. note:: Not to be confused with :any:`aggregate_time`,
              which is used to (for example) aggregate all points at
              the same hour of different days to a single point.

    :param cube: input iris cube
    :param stat: string, statistic to calculate. Available options are
                 'max', 'mean', 'min', 'sum', 'std' (standard deviation),
                 'var' (variance)
    :param period: period to aggregate over. Currently available options
                   are 'day', 'month', and 'year'.
    :param min_periods: Minimum number of data values required to calculate
                        a value - if less than this number present,
                        then value is set to NaN.
    :param aqdates: If True, then dates are set according to air quality
               standards, ie 00Z is actually 24Z from the previous day.

               If not specified (ie None), the value True or False will
               be guessed from the bounds of the first point.

    :returns: cube containing data for each `period`.

    Example:

    >>> import config
    >>> import adaq_data
    >>> sample_data_path = config.SAMPLE_DATADIR+'sites_cube_list/'
    >>> md = adaq_data.ADAQData()
    >>> scl = md.load_ts(sample_data_path+'aqum_oper_5days.nc')

    >>> mod_cube = md.extract(short_name='O3', singlecube=True)
    >>> print(mod_cube) # doctest: +ELLIPSIS
    mass_concentration_of_ozone_in_air / (ug/m3) (site_id: 5; time: 121)
         Dimension coordinates:
              site_id                                    x        -
              time                                       -        x
         Auxiliary coordinates:
              abbrev                                     x        -
              grid_latitude                              x        -
              grid_longitude                             x        -
              latitude                                   x        -
              longitude                                  x        -
              site_altitude                              x        -
              site_name                                  x        -
              site_type                                  x        -
              surface_altitude                           x        -
              forecast_period                            -        x
         Scalar coordinates:
              forecast_day: 1.0 Days
              level_height: 20.00... m, bound=(0.0, 49.99...) m
              model_level_number: 1
              sigma: 0.99..., bound=(1.0, 0.99...)
         Attributes:
              Conventions: CF-1.5
              STASH: m01s34i001
              label: aqum_oper
              short_name: O3
              source: Data from Met Office Unified Model
         Cell methods:
              mean: time (1 hour)

    Calculate the maximum for each day:

    >>> stat_cube = periodic_stat(mod_cube, "max", "day")
    >>> print(stat_cube)
    mass_concentration_of_ozone_in_air / (ug/m3) (site_id: 5; time: 6)
         Dimension coordinates:
              site_id                                    x        -
              time                                       -        x
         Auxiliary coordinates:
              abbrev                                     x        -
              grid_latitude                              x        -
              grid_longitude                             x        -
              latitude                                   x        -
              longitude                                  x        -
              site_altitude                              x        -
              site_name                                  x        -
              site_type                                  x        -
              surface_altitude                           x        -
              date                                       -        x
              forecast_period                            -        x
         Scalar coordinates:
              forecast_day: 1.0 Days
              level_height: 20.000338 m, bound=(0.0, 49.998882) m
              model_level_number: 1
              sigma: 0.9977165, bound=(1.0, 0.99429625)
         Attributes:
              Conventions: CF-1.5
              STASH: m01s34i001
              label: aqum_oper
              short_name: O3
              source: Data from Met Office Unified Model
         Cell methods:
              mean: time (1 hour)
              nanmax: date

    The new co-ordinate 'date' contains a string representation for each
    period - in this case each day:

    >>> print(stat_cube.coord("date").points)
    ['2014-03-25' '2014-03-26' '2014-03-27' '2014-03-28' '2014-03-29'
     '2014-03-30']

    Note that the time co-ordinate bounds have been merged,
    and that the time points now represent the midpoint:

    >>> print(stat_cube.coord("time")) # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    DimCoord([2014-03-25 23:30:00, 2014-03-26 12:00:00, 2014-03-27 12:00:00,
           2014-03-28 12:00:00, 2014-03-29 12:00:00, 2014-03-30 12:00:00],
           bounds=[[2014-03-25 23:00:00, 2014-03-26 00:00:00],
           [2014-03-26 00:00:00, 2014-03-27 00:00:00],
           [2014-03-27 00:00:00, 2014-03-28 00:00:00],
           [2014-03-28 00:00:00, 2014-03-29 00:00:00],
           [2014-03-29 00:00:00, 2014-03-30 00:00:00],
           [2014-03-30 00:00:00, 2014-03-31 00:00:00]],
           ...)
    """

    #Check cube has a time-coordinate
    if not cube.coords('time'):
        raise ValueError('cube does not have a time coordinate')

    #Define an appropriate aggregator
    aggregator_name = "nan"+stat
    if min_periods > 1:
        aggregator_name += "_min{}periods".format(min_periods)
    aggregator = iris.analysis.Aggregator(
        aggregator_name,
        array_statistics.nanstat_minperiod,
        stat=stat,
        min_periods=min_periods)

    #Choose a date format that groups times in the same period
    if period == 'day':
        date_fmt = '%Y-%m-%d'
    elif period == 'month':
        date_fmt = '%Y-%m'
    elif period == 'year':
        date_fmt = '%Y'
    else:
        raise ValueError("unrecognised aggregation period: {}".format(period))

    #Work on a copy of the input cube
    stat_cube = cube.copy()
    coord = stat_cube.coord('time')

    #Determine whether a point at 00Z represents data for the same day
    # or the previous day, by checking bounds. In fact only check
    # whether the first point is at the end of its bound, and assume
    # that the same applies to all points.
    if aqdates is None:
        if coord.has_bounds():
            aqdates = coord.points[0] == coord.bounds[0, 1]
        else:
            warnings.warn("aqdates could not be guessed but assumed False")
            aqdates = False

    #Add a new coordinate to aggregate by
    if aqdates:
        func = lambda coord, point: cube_time.date_from_time_aq(coord, point,
                                                                fmt=date_fmt)
    else:
        func = lambda coord, point: cube_time.date_from_time(coord, point,
                                                             fmt=date_fmt)
    coord_cat.add_categorised_coord(stat_cube, 'date', coord, func,
                                    units='no_unit')

    #Perform the aggregation
    stat_cube = stat_cube.aggregated_by('date', aggregator)

    return stat_cube


def daily_stat(cube, stat='max', min_periods=1, aqdates=False):
    """
    Calculate daily maximum cube.

    :param cube: input iris cube
    :param stat: string, statistic to calculate. Available options are
                 'max', 'mean', 'min', 'sum', 'std' (standard deviation),
                 'var' (variance)
    :param min_periods: Minimum number of data values required to calculate
                        a value - if less than this number present,
                        then value is set to NaN.
    :param aqdates: If True, then dates are set according to air quality
               standards, ie 00Z is actually 24Z from the previous day.

    >>> import config
    >>> import adaq_data
    >>> samplepath = config.SAMPLE_DATADIR+'sites_cube_list/'
    >>> od = adaq_data.ADAQData()
    >>> obs_scl = od.load_ts(samplepath+'aurn_5days.nc')
    >>> cube = od.extract(short_name='O3', singlecube=True)
    >>> maxcube = daily_stat(cube, stat='max', min_periods=18)
    >>> print(maxcube)
    mass_concentration_of_ozone_in_air / (ug/m3) (site_id: 5; time: 6)
         Dimension coordinates:
              site_id                                    x        -
              time                                       -        x
         Auxiliary coordinates:
              abbrev                                     x        -
              latitude                                   x        -
              longitude                                  x        -
              site_altitude                              x        -
              site_name                                  x        -
              site_type                                  x        -
              date                                       -        x
         Attributes:
              Conventions: CF-1.5
              label: Obs
              short_name: O3
              source: AURN
         Cell methods:
              mean: time (1 hour)
              nanmax_min18periods: date
    >>> print(cube.data[0,:24].max())
    90.0
    >>> np.set_printoptions(formatter={'float':lambda x: '{:5.2f}'.format(x)})
    >>> print(maxcube.data[0,:])
    [90.00 95.00 86.00 88.00 91.00   nan]
    >>> np.set_printoptions()
    """


    #Check cube has a time-coordinate
    try:
        cube.coord('time')
    except:
        raise ValueError('cube does not have a time coordinate')

    cube_tmp = cube.copy()

    #Firstly add the date coordinate
    if aqdates:
        coord_cat.add_categorised_coord(cube_tmp, 'date',
                                        cube_tmp.coord('time'),
                                        cube_time.date_from_time_aq,
                                        units='')
    else:
        coord_cat.add_categorised_coord(cube_tmp, 'date',
                                        cube_tmp.coord('time'),
                                        cube_time.date_from_time,
                                        units='')

    #Then can aggregate to get maximum value for each date
    aggregator_name = "nan"+stat
    if min_periods > 1:
        aggregator_name += "_min%dperiods" % (min_periods)
    max_aggregator = iris.analysis.Aggregator(
        aggregator_name,
        array_statistics.nanstat_minperiod, stat=stat,
        min_periods=min_periods)

    day_cube = cube_tmp.aggregated_by('date', max_aggregator)

    return day_cube

def distance(cube, reflong, reflat):
    """
    Calculate distance of each grid point from the point defined
    by a given longitude and latitude

    Works on a cube with in PlateCarree (WGS84) projection
    Returns a cube which is a copy of the original cube but includes
    an extra coordinate called Great Circle Distance which provides the
    distance in metres from the reference latitude and longitude
    calculated along great circles.

    :param cube: cube on which distance calculation is carried out
    :param reflong: reference longitude
    :param reflat: reference latitude

    Use iris to load some example data

    >>> import iris
    >>> import config
    >>> sample_data_path = config.SAMPLE_DATADIR+'name/'
    >>> cube = iris.load_cube(sample_data_path+'NAME_output.txt')

    Specify a reference latitude and longitude

    >>> reflat = 63.63
    >>> reflong = -19.62

    Compute great circle distances from the reference latitude and longitude

    >>> new_cube = distance(cube, reflong, reflat)

    Extract the great circle distances and print the minimum distance

    >>> distance_points = new_cube.coord('great circle distance').points
    >>> print('{:.3f}'.format(distance_points.min()))
    6979.647
    """

    # Initial checking
    if not isinstance(cube, iris.cube.Cube):
        raise ValueError('code requires a single cube')

    try:
        lons = cube.coord('longitude')
    except:
        raise ValueError('cube does not have a longitude coordinate')
    try:
        lats = cube.coord('latitude')
    except:
        raise ValueError('cube does not have a latitude coordinate')

    xpt, ypt = np.meshgrid(lons.points, lats.points)

    arr_reflong = np.ones_like(xpt) * reflong
    arr_reflat = np.ones_like(ypt) * reflat

    # Calculate great circle distance between points (using our ellipse)
    globe = lons.coord_system.as_cartopy_globe()
    geodetic = Geod(a=globe.semimajor_axis, b=globe.semiminor_axis)
    _, _, dist = geodetic.inv(xpt, ypt, arr_reflong, arr_reflat)

    # Associate with the cube
    dist_coord = iris.coords.AuxCoord(
        dist, long_name='great circle distance',
        units='m', attributes={'origin':
                               'great circle distance calculated using pyproj',
                               'distance_to':
                                   '({},{})'.format(reflong, reflat)})
    cube.add_aux_coord(dist_coord,
                       data_dims=[cube.coord_dims(lats)[0],
                                  cube.coord_dims(lons)[0]])

    return cube

def diurnal(cube, aggregator_name='NANMEAN', collapsed=False):
    """
    Calculate the diurnal variation (typically means) of a cube.
    This function is a wrapper around the :func:`aggregate_time` function.

    :param aggregator_name: string name of aggregator to use. By default this
                             is 'NANMEAN' which will take a mean over all
                             values, ignoring NaNs. Other available alternatives
                             include 'NANMAX' (max, ignoring NaNs), 'MEAN' and
                             'MAX', which are the mean
                             and max values, taking NaNs into account (if any
                             value in array is NaN, returns NaN).
                             Uses cube_statistics.CUBE_AGGREGATORS dictionary.
    :param collapsed: logical to determine whether to collapse over all
                      dimensions. If False (default), time dimension is removed
                      and replaced with hour, but all other dimensions kept.
                      If True, the same aggregator is used to collapse all other
                      dimensions. In this case,
                      the only dimension in the returned cube is time.


    Example:
    >>> import config
    >>> import adaq_data
    >>> sample_data_path = config.SAMPLE_DATADIR+'sites_cube_list/'
    >>> md = adaq_data.ADAQData()
    >>> scl = md.load_ts(sample_data_path+'aqum_oper_5days.nc')

    >>> mod_cube = md.extract(short_name = 'O3', singlecube = True)
    >>> print(mod_cube.summary(shorten=True))
    mass_concentration_of_ozone_in_air / (ug/m3) (site_id: 5; time: 121)

    Calculate diurnal variation, but without collapsing the site_id coordinate:

    >>> diurnal_cube = diurnal(mod_cube, collapsed=False)
    >>> print(diurnal_cube.summary(shorten=True))
    mass_concentration_of_ozone_in_air / (ug/m3) (hour: 24; site_id: 5)

    And now collapse the site_id coordinate
    (the nanmean is also taken across site_id):

    >>> diurnal_cube = diurnal(mod_cube, collapsed=True)
    >>> print(diurnal_cube.summary(shorten=True))
    mass_concentration_of_ozone_in_air / (ug/m3) (hour: 24)

    This returned cube has only one dimension - hour.

    """

    new_cube = aggregate_time(cube,
                              period='hour',
                              aggregator_name=aggregator_name,
                              collapsed=collapsed)

    return new_cube




def match_cubes(cubes):
    """
    Ensure times and sites (if exist) are matching
    - very important for calculating accurate statistics.
    Where-ever any cube has nan data, set the other cubes to nan

    :param cubes: List of cubes to match between each other

    :returns: List of cubes that have been matched

    >>> import config
    >>> import adaq_data
    >>> samplepath = config.SAMPLE_DATADIR+'sites_cube_list/'
    >>> od = adaq_data.ADAQData()
    >>> obs_scl = od.load_ts(samplepath+'aurn_5days.nc')
    >>> obs = od.extract(short_name='NO2', singlecube=True)

    >>> md = adaq_data.ADAQData()
    >>> mod_scl = md.load_ts(samplepath+'aqum_oper_5days.nc')
    >>> mod = md.extract(short_name='NO2', singlecube=True)

    For this example, modify mod such that some of the times at
    the end are missing:

    >>> mod = mod[:,:25]

    Have a look at the number of time points, and the data for
    the first time:

    >>> np.set_printoptions(formatter={'float':lambda x: '{:5.2f}'.format(x)})
    >>> print(len(obs.coord('time').points), obs.data[:,0])
    121 [ 8.20   nan  7.40 83.00 33.70]
    >>> print(len(mod.coord('time').points), mod.data[:,0])
    25 [ 3.39 30.56 13.56 10.31 60.98]

    Note mod has less time points than obs, plus obs has some
    missing data (the second site).

    >>> obs_new, mod_new = match_cubes([obs, mod])
    >>> print(len(obs_new.coord('time').points), obs_new.data[:,0])
    25 [ 8.20   nan  7.40 83.00 33.70]
    >>> print(len(mod_new.coord('time').points), mod_new.data[:,0])
    25 [ 3.39   nan 13.56 10.31 60.98]
    >>> np.set_printoptions()
    """

    if len(cubes) == 1:
        #Only 1 cube, so don't need to match - return input cube
        return cubes

    cube1_coord_names = [coord.name() for coord in cubes[0].coords()]

    #Match times
    if 'time' in cube1_coord_names:
        cubes = cube_time.intersect_cubetime(cubes)

    #Match site ids
    if 'site_id' in cube1_coord_names:
        site_ids = []
        extract = False #Logical to test if site ids are actually different
                        #Don't need to use cube extract (slow)
                        # if site ids are the same.
        for i, cube in enumerate(cubes):
            if i == 0:
                site_ids = set(cube.coord('site_id').points)
            else:
                if site_ids != set(cube.coord('site_id').points):
                    extract = True
                    site_ids = (site_ids & set(cube.coord('site_id').points))
        if extract:
            site_constraint = iris.Constraint(site_id=site_ids)
            newcubes = []
            for cube in cubes:
                newcubes.append(cube.extract(site_constraint))
            cubes = newcubes

    #Check all cubes have the same size/shape data array as the first cube
    for cube in cubes[1:]:
        assert cube.data.shape == cubes[0].data.shape

    #Match nan points
    #First find all indices where we have any nan data
    for i, cube in enumerate(cubes):
        if i == 0:
            nans = np.isnan(cube.data)
        else:
            nans = nans | np.isnan(cube.data)
    indices = np.where(nans)
    #Then set all these points to nan
    for cube in cubes:
        cube.data[indices] = np.nan

    return cubes

def maximum_daily_rolling_8hr_mean(cube, aqdates=True):
    """
    Calculate maximum daily rolling 8 hour mean, as required by Defra's DAQI.
    An 8 hour mean is defined as the average of that hour and
    the 7 previous hours,
    but only defined if at least 75% of the data (ie 6 hours) is available.
    The maximum of all the 8 hour means (defined as 1Z-24Z if aq=True,
    0Z-23Z otherwise) is then calculated, but only if at least
    75% (ie 18hours) of data available.
    Returns a cube with one time per day (instead of per hour).

    >>> import config
    >>> import adaq_data
    >>> samplepath = config.SAMPLE_DATADIR+'sites_cube_list/'
    >>> od = adaq_data.ADAQData()
    >>> obs_scl = od.load_ts(samplepath+'aurn_5days.nc')
    >>> cube = od.extract(short_name='O3', singlecube=True)
    >>> print(cube)
    mass_concentration_of_ozone_in_air / (ug/m3) (site_id: 5; time: 121)
         Dimension coordinates:
              site_id                                    x        -
              time                                       -        x
         Auxiliary coordinates:
              abbrev                                     x        -
              latitude                                   x        -
              longitude                                  x        -
              site_altitude                              x        -
              site_name                                  x        -
              site_type                                  x        -
         Attributes:
              Conventions: CF-1.5
              label: Obs
              short_name: O3
              source: AURN
         Cell methods:
              mean: time (1 hour)
    >>> np.set_printoptions(formatter={'float':lambda x: '{:5.2f}'.format(x)})
    >>> print(cube.data[0,:25])
    [49.00 48.00 49.00 45.00 44.00 44.00 45.00 44.00 42.00 49.00 54.00 60.00
     65.00 73.00 84.00 88.00 90.00 87.00 87.00 82.00 67.00 63.00 63.00 49.00
     54.00]

    >>> cube_md8m = maximum_daily_rolling_8hr_mean(cube, aqdates=True)
    >>> print(cube_md8m)
    mass_concentration_of_ozone_in_air / (ug/m3) (site_id: 5; time: 6)
         Dimension coordinates:
              site_id                                    x        -
              time                                       -        x
         Auxiliary coordinates:
              abbrev                                     x        -
              latitude                                   x        -
              longitude                                  x        -
              site_altitude                              x        -
              site_name                                  x        -
              site_type                                  x        -
              date                                       -        x
         Attributes:
              Conventions: CF-1.5
              label: Obs
              short_name: O3
              source: AURN
         Cell methods:
              mean: time (1 hour)
              nanmax_min18periods: date
    >>> print(cube_md8m.data[0,:2])
    [  nan 82.25]
    >>> np.set_printoptions()

    Note the first value is nan. This is due to the first value in the input
    cube having a time of 00Z - this corresponds to the previous day
    when aq=True, hence this date is included in the output cube,
    but the data value is nan as only one time is included
    (not the minumum of 18 points):

    >>> tunits = cube.coord('time').units
    >>> np.set_printoptions(linewidth=70)
    >>> print(tunits.num2date(cube.coord('time')[:2].points))  # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    [...datetime(2014, 3, 26, 0, 0)
    ...datetime(2014, 3, 26, 1, 0)]
    >>> np.set_printoptions(linewidth=75)
    >>> print(cube_md8m.coord('date')[:2].points)
    ['2014-03-25' '2014-03-26']
    """

    #Check cube has a time-coordinate
    try:
        cube.coord('time')
    except:
        raise ValueError('cube does not have a time coordinate')

    #Check hourly data
    dtpts = cube.coord('time').units.num2date(cube.coord('time').points)
    delta = min(dtpts[1:] - dtpts[:-1])
    assert delta == datetime.timedelta(hours=1)

    #Check no missing times and add in if needed (can be slow)
    #cube = add_missing_times(cube)

    #Calculate rolling 8 hour mean.
    #Requires a minimum of 6 hours of data for each 8hr period.

    #Get a rolling window - this contains a view to the data
    # which has an extra axis of the length of the window, containing
    # the data for that time and the previous times
    taxis = cube.coord_dims('time')[0]
    rw = array_statistics.rolling_window_extended(cube.data, 8, axis=taxis)
    #Then use an aggregator to calculate the means from each of these slices
    mean_aggregator = iris.analysis.Aggregator(
        "nanmean_min%dperiods"%(8),
        array_statistics.nanstat_minperiod,
        stat='mean', min_periods=6)
    data = mean_aggregator.aggregate(rw, axis=taxis+1)
    #Take copy of cube to ensure data not overwritten,
    #  then overwrite data array with this meaned data
    cube_mean = cube.copy()
    cube_mean.data = data

    #Now calculate daily maximum of these means
    #Requires a minimum of 18hrs of means for each day.
    day_cube = daily_stat(cube_mean, stat='max',
                          min_periods=18, aqdates=aqdates)

    return day_cube


if __name__ == '__main__':

    import doctest
    doctest.testmod()
