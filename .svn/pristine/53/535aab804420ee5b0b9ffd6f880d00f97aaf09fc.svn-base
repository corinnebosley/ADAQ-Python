"""
Contains generic functions relating to time which work on any cube.
"""
from six.moves.builtins import zip
from six.moves.builtins import str
from six.moves.builtins import range

import datetime
import warnings

import numpy as np
import numpy.ma as ma
import iris

from distutils.version import LooseVersion
if LooseVersion(iris.__version__) <= '1.13.0':
    iris.FUTURE.cell_datetime_objects = True

def cube_add_missing_times(input_cube):
    """
    For use with cubes with regular time data, but where some
    times are completely missing from the cube.
    This routine will add the extra times, by assuming a regular
    delta derived from the minimum distance between two consecutive times.
    Any data associated with these extra times are set to NaN.
    Note if bounds are added back with guess_bounds, so if
    they are not regular in input cube, the output may not be correct.
    Also note, 360 day calendar not supported.

    Get a simple cube in (although will work also with multiple dimensioned
    cubes)

    >>> import config
    >>> samplepath = config.SAMPLE_DATADIR+'sites_cube_list/'
    >>> constraint = (iris.AttributeConstraint(short_name="O3") &
    ...     iris.Constraint(abbrev="HAR"))
    >>> cube = iris.load_cube(samplepath+"aurn_1days.nc",constraint)
    >>> cube = iris.cube.CubeList([cube[:4],cube[8:12]]).concatenate()[0]
    >>> np.set_printoptions(formatter={'float':lambda x: '{:5.2f}'.format(x)})
    >>> print(cube.data)
    [16.00 20.00 37.00 39.00 31.00 40.00 37.00 23.00]
    >>> print(cube.coord('time')) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
    DimCoord([2014-04-02 00:00:00, 2014-04-02 01:00:00, 2014-04-02 02:00:00,
           2014-04-02 03:00:00, 2014-04-02 08:00:00, 2014-04-02 09:00:00,
           2014-04-02 10:00:00, 2014-04-02 11:00:00], \
bounds=[[2014-04-01 23:00:00, 2014-04-02 00:00:00],
           [2014-04-02 00:00:00, 2014-04-02 01:00:00],
           [2014-04-02 01:00:00, 2014-04-02 02:00:00],
           [2014-04-02 02:00:00, 2014-04-02 03:00:00],
           [2014-04-02 07:00:00, 2014-04-02 08:00:00],
           [2014-04-02 08:00:00, 2014-04-02 09:00:00],
           [2014-04-02 09:00:00, 2014-04-02 10:00:00],
           [2014-04-02 10:00:00, 2014-04-02 11:00:00]], \
standard_name=...'time', calendar=...'gregorian', var_name='time')

    Notice that the hours 4:00 - 7:00 are missing.

    >>> newcube = cube_add_missing_times(cube)
    >>> print(newcube.coord('time'))  # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
    DimCoord([2014-04-02 00:00:00, 2014-04-02 01:00:00, 2014-04-02 02:00:00,
           2014-04-02 03:00:00, 2014-04-02 04:00:00, 2014-04-02 05:00:00,
           2014-04-02 06:00:00, 2014-04-02 07:00:00, 2014-04-02 08:00:00,
           2014-04-02 09:00:00, 2014-04-02 10:00:00, 2014-04-02 11:00:00], \
bounds=[[2014-04-01 23:00:00, 2014-04-02 00:00:00],
           [2014-04-02 00:00:00, 2014-04-02 01:00:00],
           [2014-04-02 01:00:00, 2014-04-02 02:00:00],
           [2014-04-02 02:00:00, 2014-04-02 03:00:00],
           [2014-04-02 03:00:00, 2014-04-02 04:00:00],
           [2014-04-02 04:00:00, 2014-04-02 05:00:00],
           [2014-04-02 05:00:00, 2014-04-02 06:00:00],
           [2014-04-02 06:00:00, 2014-04-02 07:00:00],
           [2014-04-02 07:00:00, 2014-04-02 08:00:00],
           [2014-04-02 08:00:00, 2014-04-02 09:00:00],
           [2014-04-02 09:00:00, 2014-04-02 10:00:00],
           [2014-04-02 10:00:00, 2014-04-02 11:00:00]], \
standard_name=...'time', calendar=...'gregorian', var_name='time')

    >>> print(newcube.data)
    [16.00 20.00 37.00 39.00   nan   nan   nan   nan 31.00 40.00 37.00 23.00]
    >>> np.set_printoptions()

    Notice now the extra hours 4:00 - 7:00 now exist, but their data values are
    set to nan.
    """
    #Notes:
    # Written with support from AVD!
    # This code takes cubes with regular time data, but where some times are
    # completely missing from the cube.
    # It will expand the cube to add these extra times in,
    # setting their data values to nan.
    # This is a non-trivial problem as the time dimension could be any
    # dimension, not always the first, or last dimension for example.
    # It also has to take into account that there may be other coordinates
    # on this dimension, for example forecast period.

    #Method:
    # 1. Figure out which dimension number time is on.
    # 2. Convert the points from the input cube time coordinate to date-time
    #    format
    # 3. Calculate the required time points for the output cube, based on the
    #    start and end points of the input cube and the smallest time delta
    #    between input points.
    # 4. If the output points are exactly the same as the input points, then
    #    there is nothing to do, so can escape from this routine early,
    #    returning the input cube.
    # 5. Generate a boolean numpy array containing True where the required
    #    output points exist already in the input array.
    #    Set to False where the output points don't exist.
    # 6. Generate a list of indices for all output time points to indicate which
    #    indice this corresponds to in the input time points.
    #    Where the point does not exist in the input times, set the indice to
    #    zero to temporarily copy from from the first (zero'th) element.
    # 7. Expand the input cube.
    #    a. Set up an indexing array, where the indices in the time dimension
    #       match the required indice from the input cube.
    #    b. Set up the output cube by expanding the input cube - using values
    #       from the input cube where they exist, or pointing to the first index
    #       of the input cube otherwise.
    #    c. Turn the indexing array into a True/False array, setting to True
    #       indices which were not present in the input array, or
    #       False if they were present in the input array
    #    d. Now can set the data values whose indices are True to np.nan (SLOW)
    #       This therefore sets those points which were missing from the input
    #       array to nan.
    # 8. Fix any other coordinates that are on the same dimension as time as the
    #    newly added points match the zero'th element from the input cube but
    #    should be reset to be regular.
    # 9. Ensure time is still a dimension coordinate.

    #---
    #  1. Figure out which dimension number time is on.

    #Get time coordinate
    time_coord = input_cube.coord('time')
    #Get time dimension
    time_dims = input_cube.coord_dims(time_coord)
    assert len(time_dims) == 1
    i_time_dim = time_dims[0]

    #Access shape of data (fails otherwise!)
    #Pylint ignore "Statement has no effect error"
    #pylint: disable=W0104
    input_cube.data.shape
    #pylint: enable=W0104

    #---
    # 2. Convert the points from the input cube time coordinate to date-time
    #    format
    time_units = time_coord.units
    input_cube_tpts = time_coord.points
    input_cube_dtpts = time_units.num2date(input_cube_tpts)

    #---
    # 3. Calculate the required time points for the output cube, based on the
    #    start and end points of the input cube and the smallest time delta
    #    between input points.
    time_step = np.min(np.diff(input_cube_tpts))
    n_tpts = (input_cube_tpts[-1] - input_cube_tpts[0]) / time_step
    if not np.isclose(n_tpts, int(n_tpts)):
        #ie the smallest delta does not divide the full range
        #Continuing would likely raise an error, so escape with a warning
        warnings.warn("skipping cube with irregular time steps")
        return input_cube
    n_tpts = 1 + int(np.round(n_tpts))
    output_cube_tpts = np.linspace(input_cube_tpts[0],
                                   input_cube_tpts[-1],
                                   n_tpts)
    output_cube_dtpts = time_units.num2date(output_cube_tpts)
    #NB need to use datetime point values instead of raw point values:
    # To ensure equality in time points, eg for setting up mask correctly
    # To ensure linear interpolation gives exactly the same data
    #  values back, which is not the case when using newcube_tpts instead
    #  as 15. changes to eg 14.99999999

    #---
    # 4. If the output points are exactly the same as the input points, then
    #    there is nothing to do, so can escape from this routine early,
    #    returning the input cube.
    if len(output_cube_dtpts) == len(input_cube_dtpts):
        if all(output_cube_dtpts == input_cube_dtpts):
            return input_cube

    #---
    # 5. Generate a boolean numpy array containing True where the required
    #    output points exist already in the input array.
    #    Set to False where the output points don't exist.
    tpts_valid = np.array([dtpt in input_cube_dtpts
                           for dtpt in output_cube_dtpts])

    #---
    # 6. Generate a list of indices for all output time points to indicate which
    #    indice this corresponds to in the input time points.
    #    Where the point does not exist in the input times, set the indice to
    #    zero to temporarily copy from from the first (zero'th) element.
    #    Nb the data at this indice will later be set to nan
    #    This will be used for expanding the input data array.
    input_time_indices = []
    index = 0
    for tpt_valid in tpts_valid:
        if tpt_valid:
            #In input cube
            input_time_indices.append(index)
            index += 1 #Increment index
        else:
            #Not in input cube: set to zero
            input_time_indices.append(0)

    #---
    # 7. Expand the input cube.

    # Note on slices:
    #   array[slice(start,stop,step)] == array[start:stop:step]
    #   array[slice(None)] == array[slice(None,None,None)] == array[:]
    #   array[slice(None), slice(None) ] == array[:,:]
    #   array[slice(None), (0,1,2,3,4) ] == array[:,0:5]
    #   if array[1] = [3,4,5], then array[slice(None), (0,1,0,2,0)]
    #     gives array[1] = [3,4,3,5,3] which has been expanded.

    # 7a. Set up an indexing array, where the indices in the time dimension
    #     match the required indice from the input cube.
    # Start with a "slice(None)" (equivalent to :) for each key
    # up to the time dimension.
    cube_all_indices = [slice(None)] * (i_time_dim + 1)
    # Replace the time key with our sequence of values (it must be a tuple).
    cube_all_indices[i_time_dim] = tuple(input_time_indices)
    # 7b. Set up the output cube by expanding the input cube - using values
    #     from the input cube where they exist, or pointing to the first index
    #     of the input cube otherwise.
    # Where values exist in the input cube,
    # they are correctly copied to the output cube
    # Where values don't exist,
    # temporarily set to the first value from the input cube
    output_cube = input_cube[tuple(cube_all_indices)]
    # 7c. Turn the indexing array into a True/False array, setting to True
    #     indices which were not present in the input array, or
    #     False if they were present in the input array
    cube_all_indices[i_time_dim] = ~tpts_valid
    # 7d. Now can set the data values whose indices are True to np.nan (SLOW)
    #     This therefore sets those points which were missing from the input
    #     array to nan.
    #!!! NB this is the slowest line of the code as finally accessing the large
    #    data array from the cube...
    output_cube.data[tuple(cube_all_indices)] = np.nan

    #---
    # 8. Fix any other coordinates that are on the same dimension as time as the
    #    newly added points match the zero'th element from the input cube but
    #    should be reset to be regular and as expected.
    for coord in output_cube.coords(contains_dimension=i_time_dim):
        coord_name = coord.name()
        coord.points = np.linspace(coord.points[0],
                                   coord.points[-1],
                                   coord.points.size)
        if input_cube.coord(coord_name).has_bounds():
            bounds = np.empty(coord.bounds.shape)
            bounds[:, 0] = np.linspace(coord.bounds[0, 0],
                                       coord.bounds[-1, 0],
                                       coord.bounds[:, 0].size)
            bounds[:, 1] = np.linspace(coord.bounds[0, 1],
                                       coord.bounds[-1, 1],
                                       coord.bounds[:, 1].size)
            coord.bounds = bounds

    #---
    # 9. Ensure time is still a dimension coordinate.

    #Finally, make time-coordinates back into dimensions coords
    #  (if dim coord in input_cube)
    for coord in input_cube.coords(contains_dimension=i_time_dim,
                                   dim_coords=True):
        iris.util.promote_aux_coord_to_dim_coord(output_cube, coord.name())

    return output_cube


def cube_tpoints_dt(cube):
    """
    Works on any cube which has a time coordinate.
    Returns points from the time-coordinate in datetime format

    Get a cube and examine standard output for points:

    >>> import config
    >>> sample_data_path = config.SAMPLE_DATADIR+'ukhires/'
    >>> cubes = iris.load(sample_data_path+'uk_hires.pp')
    >>> np.set_printoptions(formatter={'float':lambda x: '{:5.2f}'.format(x)})
    >>> print(cubes[0].coord('time').points)
    [349618.00 349619.00 349620.00]
    >>> np.set_printoptions()

    Now use cube_tpoints_dt to print in nice date-time format:

    >>> dtpoints = cube_tpoints_dt(cubes[0])
    >>> print(dtpoints)  # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    [...datetime(2009, 11, 19, 10, 0)
     ...datetime(2009, 11, 19, 11, 0)
     ...datetime(2009, 11, 19, 12, 0)]
    """

    #Check cube has time-coordinate
    try:
        cube.coords('time')
        tcoord = cube.coord('time')
    except:
        raise ValueError('cube does not have a time coordinate')
    tunits = tcoord.units
    dtpoints = tunits.num2date(tcoord.points)

    return dtpoints


def date_from_time(coord, point, fmt='%Y-%m-%d'):
    """
    Calculate date string from the coordinate point value.
    """
    return coord.units.num2date(point).strftime(fmt)

def date_from_time_aq(coord, point, fmt='%Y-%m-%d'):
    """
    Calculate date string.
    As this is for use for air quality purposes, 24Z should be treated
    as the same date as 23Z.
    """

    dt = coord.units.num2date(point)
    if dt.hour == 0 and dt.minute == 0:
        #24Z => Date should be from previous day (same as 23Z)
        dt -= datetime.timedelta(hours=1)
    return dt.strftime(fmt)

def hour_from_time(coord, point):
    """
    Category function to calculate hour given time, for use in
    add_categorised_coord.
    """
    return coord.units.num2date(point).hour

def month_from_time(coord, point):
    """
    Category function to calculate month given time, for use in
    add_categorised_coord.
    """
    return coord.units.num2date(point).month

def monthname_from_month(coord, point):
    """
    Category function to calculate month name given integer month, for use in
    add_categorised_coord. Returns three letter month name, eg 'Feb'
    """
    #Use arbitary year and date within month as these are not used
    return datetime.datetime(2000, point, 1).strftime('%b')

def yearmonth_from_time(coord, point):
    """
    Category function to calculate year and month given time, for use in
    add_categorised_coord. Returns integer in format YYYYmm, eg 201402.
    """
    return int(coord.units.num2date(point).strftime('%Y%m'))

def yearmonthname_from_yearmonth(coord, point):
    """
    Category function to calculate year and month name (month Year eg
    'Feb 2014') given yearmonth (integer in format YYYYmm, eg 201402),
    for use in add_categorised_coord.
    """
    year = int(str(point)[:4])
    month = int(str(point)[4:6])
    #Use arbitary date (1st) within month as these are not used

    return datetime.datetime(year, month, 1).strftime('%b %Y')


def extract_all_forecast_days(cubelist):
    """
    Extract all forecast days, putting each forecast day into a separate cube
    within the output cubelist.

    :param cubelist: Input iris.cube.CubeList (can also take a single cube which
                     is immediately converted to a cubelist)

    :returns: cubelist - iris.cube.CubeList with >= number of input cubes,
              each cube only containing a single forecast_day. Time is promoted
              to a dimension coordinate.

    .. note:: It is not guaranteed that the returned order for a single cube
              name will be in forecast_day order.

    Get some example data from pp files, which have a forecast_day coordinate:

    >>> import config
    >>> import pp_data
    >>> pp = pp_data.PPData()
    >>> start_dt = datetime.datetime(2014, 4, 3, 00)
    >>> end_dt = datetime.datetime(2014, 4, 7, 00)
    >>> filenames=config.SAMPLE_DATADIR+'aqum_output/oper_forecast/*201404*.pp'
    >>> gcl = pp.readdata(filenames=filenames, short_name_list=['O3'],
    ...    start_datetime=start_dt, end_datetime=end_dt, forecast_day=None)
    >>> print(gcl)
    0: mass_fraction_of_ozone_in_air / (kg kg-1) \
(-- : 267; grid_latitude: 182; grid_longitude: 146)
    >>> np.set_printoptions(formatter={'float':lambda x: '{:5.2f}'.format(x)})
    >>> print(gcl[0].coord('forecast_day').points)
    [ 1.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00
      2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00
      2.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00
      3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00
      3.00  4.00  4.00  4.00  4.00  4.00  4.00  4.00  4.00  4.00  4.00  4.00
      4.00  4.00  4.00  4.00  4.00  4.00  4.00  4.00  4.00  4.00  4.00  4.00
      4.00  5.00  5.00  5.00  5.00  5.00  5.00  5.00  5.00  5.00  5.00  5.00
      5.00  5.00  5.00  5.00  5.00  5.00  5.00  5.00  5.00  5.00  5.00  5.00
      5.00  0.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00
      1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00
      1.00  1.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00
      2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00
      2.00  2.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00
      3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00
      3.00  3.00  4.00  4.00  4.00  4.00  4.00  4.00  4.00  4.00  4.00  4.00
      4.00  4.00  4.00  4.00  4.00  4.00  4.00  4.00  4.00  4.00  4.00  4.00
      4.00  4.00  0.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00
      1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00
      1.00  1.00  1.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00
      2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00
      2.00  2.00  2.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00
      3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00
      3.00  3.00  3.00]

    >>> allcubes = extract_all_forecast_days(gcl)
    >>> for cube in allcubes:
    ...    print(cube.summary(True))
    ...    np.set_printoptions(formatter={'float':lambda x: '{:5.2f}'.format(x)})
    ...    print(np.unique(cube.coord('forecast_day').points))
    mass_fraction_of_ozone_in_air / (kg kg-1) \
(time: 2; grid_latitude: 182; grid_longitude: 146)
    [ 0.00]
    mass_fraction_of_ozone_in_air / (kg kg-1) \
(time: 49; grid_latitude: 182; grid_longitude: 146)
    [ 1.00]
    mass_fraction_of_ozone_in_air / (kg kg-1) \
(time: 72; grid_latitude: 182; grid_longitude: 146)
    [ 2.00]
    mass_fraction_of_ozone_in_air / (kg kg-1) \
(time: 72; grid_latitude: 182; grid_longitude: 146)
    [ 3.00]
    mass_fraction_of_ozone_in_air / (kg kg-1) \
(time: 48; grid_latitude: 182; grid_longitude: 146)
    [ 4.00]
    mass_fraction_of_ozone_in_air / (kg kg-1) \
(time: 24; grid_latitude: 182; grid_longitude: 146)
    [ 5.00]
    >>> np.set_printoptions()
    """

    if isinstance(cubelist, iris.cube.Cube):
        cubelist = iris.cube.CubeList([cubelist])
    assert isinstance(cubelist, iris.cube.CubeList)

    output_cubelist = iris.cube.CubeList()
    for cube in cubelist:
        #This routine requires forecast_day coordinate, so check
        #that these are available before continuing.
        cube_coord_names = [coord.name() for coord in cube.coords()]
        assert 'time' in cube_coord_names
        assert 'forecast_day' in cube_coord_names
        forecast_days = np.unique(cube.coord('forecast_day').points)
        for fd in forecast_days:
            fd_cube = cube.extract(iris.Constraint(forecast_day=fd))
            #Check that time coordinate is monotonic
            if not fd_cube.coord('time').is_monotonic():
                #Make monotonic:
                #By breaking each time into a separate cube, the
                #iris merge puts them back together in a monotonic order.
                time_cubelist = iris.cube.CubeList()
                for cube_slice in fd_cube.slices_over('time'):
                    time_cubelist.append(cube_slice)
                fd_cube = time_cubelist.merge_cube()
            #Check that is now monotonic:
            assert fd_cube.coord('time').is_monotonic()
            #If more than a single time point, then then promote to ensure the
            #axis is named 'time' rather than being anonymous
            if len(fd_cube.coord('time').points) > 1:
                iris.util.promote_aux_coord_to_dim_coord(fd_cube, 'time')
            output_cubelist.append(fd_cube)

    return output_cubelist

def extract_latest_forecast_days(cube, forecast_day='latest', start_dt=None):
    """
    Extract a cube with a montonic time dimension which depending on the
    setting of the 'forecast_day' parameter may contain a stretch of
    day 1 forecasts, plus then a forecast from a single forecast run at
    the end.

    :param cube: Input cube to extract data from. Must have a time and a
                 forecast_day coordinate.
    :param forecast_day: * 'latest' generates a cube which has day 1 forecasts
                           where possible, followed by a multi-day forecast from
                           a single forecast run
                         * 'forecast' generates a cube which only has a
                           multi-day forecat from a single forecast run.
    :param start_dt: datetime formatted date-time to start from. Only used if
                     forecast_day='forecast'. This then represents the start
                     day of the first full day of the forecast.
                     If this is instead set to None, then using
                     forecast_day='forecast' gets data from the final forecast
                     that is available from the input cube.

    Note this routine is tested in doctests of pp_data.py and maccens_data.py

    Get some example data from pp files, which have a forecast_day coordinate:

    >>> import config
    >>> import pp_data
    >>> pp = pp_data.PPData()
    >>> start_dt = datetime.datetime(2014, 4, 3, 00)
    >>> end_dt = datetime.datetime(2014, 4, 7, 00)
    >>> filenames=config.SAMPLE_DATADIR+'aqum_output/oper_forecast/*201404*.pp'
    >>> gcl = pp.readdata(filenames=filenames, short_name_list=['O3'],
    ...    start_datetime=start_dt, end_datetime=end_dt, forecast_day=None)
    >>> cube = gcl[0]
    >>> np.set_printoptions(formatter={'float':lambda x: '{:5.2f}'.format(x)})
    >>> print(cube.coord('forecast_day').points)
    [ 1.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00
      2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00
      2.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00
      3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00
      3.00  4.00  4.00  4.00  4.00  4.00  4.00  4.00  4.00  4.00  4.00  4.00
      4.00  4.00  4.00  4.00  4.00  4.00  4.00  4.00  4.00  4.00  4.00  4.00
      4.00  5.00  5.00  5.00  5.00  5.00  5.00  5.00  5.00  5.00  5.00  5.00
      5.00  5.00  5.00  5.00  5.00  5.00  5.00  5.00  5.00  5.00  5.00  5.00
      5.00  0.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00
      1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00
      1.00  1.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00
      2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00
      2.00  2.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00
      3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00
      3.00  3.00  4.00  4.00  4.00  4.00  4.00  4.00  4.00  4.00  4.00  4.00
      4.00  4.00  4.00  4.00  4.00  4.00  4.00  4.00  4.00  4.00  4.00  4.00
      4.00  4.00  0.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00
      1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00
      1.00  1.00  1.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00
      2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00
      2.00  2.00  2.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00
      3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00
      3.00  3.00  3.00]

    Now try extracting the 'latest' available forecast days:

    >>> latestcube = extract_latest_forecast_days(cube, 'latest')
    >>> print(latestcube.coord('forecast_day').points)
    [ 1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00
      1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00
      1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00
      1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00
      1.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00
      2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00
      2.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00
      3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00
      3.00]

    And then a cube with just a single forecast run in:

    >>> fcstcube = extract_latest_forecast_days(cube, 'forecast',
    ... start_dt=start_dt)
    >>> print(fcstcube.coord('forecast_day').points)
    [ 0.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00
      1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00
      1.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00
      2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00  2.00
      2.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00
      3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00
      3.00  4.00  4.00  4.00  4.00  4.00  4.00  4.00  4.00  4.00  4.00  4.00
      4.00  4.00  4.00  4.00  4.00  4.00  4.00  4.00  4.00  4.00  4.00  4.00
      4.00]
    >>> np.set_printoptions()

    """

    if forecast_day not in ['latest', 'forecast']:
        return cube

    #This routine requires both time and forecast_day coordinates, so check
    #that these are available before continuing.
    cube_coord_names = [coord.name() for coord in cube.coords()]
    assert 'time' in cube_coord_names
    assert 'forecast_day' in cube_coord_names

    tunits = cube.coord('time').units

    #Set up cubelist to contain partial cubes which will be
    #joined to become single cube at the end
    cubelist = iris.cube.CubeList()

    #Take day 1 where possible
    day1_cube = cube.extract(iris.Constraint(forecast_day=1))

    if forecast_day == 'latest':
        #Convert this day1_cube into a cubelist,
        #to ensure time can be made monotonic
        for cube_slice in day1_cube.slices_over('time'):
            cubelist.append(cube_slice)

    #Get times after final day 1
    #First find the final time from the day 1 cube
    day1_maxt = max(day1_cube.coord('time').points)
    #And convert this to datetime format so can be used in iris constraints
    day1_maxt_dt = tunits.num2date(day1_maxt)

    if forecast_day == 'latest':
        #Extract another cube, containing all times after the final time
        #from the day 1 cube, only for forecast_days > 1 as =1 is in day1_cube
        other_days_cube = cube.extract(
            iris.Constraint(forecast_day=lambda fd: fd > 1) &
            iris.Constraint(time=lambda t: t.point > day1_maxt_dt))
    elif forecast_day == 'forecast':
        #Extract a cube which contains all times after the time at the
        #end of the last day1 - 1 day.
        if start_dt is not None:
            #time at end of first required day
            day1_maxt_dt = start_dt + datetime.timedelta(days=1)

        #Don't limit by forecast day here as we want all forecasts days from
        #single forecast run, including possibly day 0 (if available)
        other_days_cube = cube.extract(
            iris.Constraint(time=lambda t:
                            t.point >= day1_maxt_dt-datetime.timedelta(days=1)))


    #If a cube is found:
    if other_days_cube is not None:
        #Check if multiple points available for each time
        if not iris.util.monotonic(other_days_cube.coord('time').points,
                                   strict=True):

            #Multiple points available for each time,
            # eg from different forecast runs
            #So now loop through each available forecast_day and pick out
            # points in cube for the requested forecast_day and with the time
            # which corresponds to the expected time in relation to the last
            # time from the day 1 cube.
            fd_max = other_days_cube.coord('forecast_day').points.max()
            if forecast_day == 'latest':
                fd_range = range(2, int(fd_max+1))
            elif forecast_day == 'forecast':
                fd_range = range(0, int(fd_max+1))
            for fd in fd_range:
                fd_cube = other_days_cube.extract(
                    iris.Constraint(forecast_day=fd) & \
                    iris.Constraint(time=lambda t:
                                    day1_maxt_dt+datetime.timedelta(days=fd-2) \
                                    < t.point <= day1_maxt_dt + \
                                    datetime.timedelta(days=fd-1)))
                if fd_cube is not None:
                    #Add individual time slices to cubelist to ensure
                    #when merged time will be montonic
                    for cube_slice in fd_cube.slices_over('time'):
                        cubelist.append(cube_slice)
        else:
            #Only a single point is available for each time,
            # so use this cube in its entirety
            if not cubelist:
                cubelist.append(other_days_cube)
            else:
                #Add individual time slices to cubelist to ensure
                #when merged with day1_cube time will be montonic
                for cube_slice in other_days_cube.slices_over('time'):
                    cubelist.append(cube_slice)

    #Merge cubelist back into a single cube
    if len(cubelist) == 1:
        #Only one cube, so don't need to merge, just take first cube.
        newcube = cubelist[0]
    elif len(cubelist) > 1:
        newcube = cubelist.merge_cube()
    else:
        newcube = None

    return newcube


def get_startenddt(cube):
    """
    Works on any cube to return start and end times in datetime-format.
    Note, uses start of first bound, and end of last bound if it has bounds,
    otherwise uses points.

    >>> import config
    >>> sample_data_path = config.SAMPLE_DATADIR+'ukhires/'
    >>> cubes = iris.load(sample_data_path+'uk_hires.pp')
    >>> cube = cubes[0]
    >>> cube.coord('time').guess_bounds()
    >>> np.set_printoptions(formatter={'float':lambda x: '{:5.2f}'.format(x)})
    >>> print(cube.coord('time').bounds)
    [[349617.50 349618.50]
     [349618.50 349619.50]
     [349619.50 349620.50]]
    >>> np.set_printoptions()
    >>> startdt, enddt = get_startenddt(cube)
    >>> print(startdt)
    2009-11-19 09:30:00
    >>> print(enddt)
    2009-11-19 12:30:00

    """
    try:
        cube.coords('time')
        tcoord = cube.coord('time')
    except:
        raise ValueError('cube does not have a time coordinate')
    tunits = tcoord.units
    if tcoord.has_bounds():
        startdt = tunits.num2date(tcoord.bounds[0][0])
        enddt = tunits.num2date(tcoord.bounds[-1][1])
    else:
        startdt = tunits.num2date(tcoord.points[0])
        enddt = tunits.num2date(tcoord.points[-1])
    return startdt, enddt


def intersect_cubetime(cubes):
    """
    Works on any cube which has a time coordinate, including a sitescube.
    Returns cubes which contain the same times by using the intersection of the
    times from all cubes.
    :param cubes: list of cubes which require intersecting.

    Load cubes with different (but intersecting) time coordinates:

    >>> import config
    >>> sample_data_path = config.SAMPLE_DATADIR+'ukhires/'
    >>> cubes = iris.load(sample_data_path+'uk_hires.pp')
    >>> cube1 = cubes[0]
    >>> cube2 = cubes[1]
    >>> print(cube_tpoints_dt(cube1))  # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    [...datetime(2009, 11, 19, 10, 0)
     ...datetime(2009, 11, 19, 11, 0)
     ...datetime(2009, 11, 19, 12, 0)]
    >>> print(cube_tpoints_dt(cube2))  # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    [...datetime(2009, 11, 19, 10, 0)]

    Now make use of intersect_cubetime to return the same cubes where
    the time-coordinate is the same for both.

    >>> newcube1, newcube2 = intersect_cubetime([cube1, cube2])
    >>> print(cube_tpoints_dt(newcube1))  # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    [...datetime(2009, 11, 19, 10, 0)]
    >>> print(cube_tpoints_dt(newcube2))  # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    [...datetime(2009, 11, 19, 10, 0)]


    """

    if len(cubes) == 1:
        #Only 1 cube, so nothing to intersect
        return cubes

    extract = False
    for i, cube in enumerate(cubes):
        try:
            cube.coords('time')
            time = cube.coord('time')
        except:
            raise ValueError('cube does not have a time coordinate')

        if i == 0:
            units = time.units
            times = units.num2date(time.points)
        else:
            #Check time units are the same
            if units == time.units:
                pass
            elif units.calendar == 'gregorian' and \
                 time.units.calendar == 'standard':
                #Gregorian and Standard are equivalent
                pass
            elif units.calendar == 'standard' and \
                 time.units.calendar == 'gregorian':
                pass
            else:
                raise ValueError('Time Units do not match')

            if set(times) != set(units.num2date(time.points)):
                #Will need to extract newcubes as sets are different
                extract = True

                #Find intersection
                times = (set(times) & set(units.num2date(time.points)))

    if extract:

        #Sort into numerical order
        times = sorted(times)

        #Set up time constraint
        tconstraint = iris.Constraint(time=lambda cell:
                                      cell.point in times)

        #Extract using time constraint
        newcubes = []
        for cube in cubes:
            newcubes.append(cube.extract(tconstraint))

    else:
        #Don't need to extract - times are already the same
        newcubes = cubes

    return newcubes


def plume_arrival_time(cube, threshold):
    '''
    Compute the arrival time of a plume at all locations using a threshold
    value to determine when the plume has arrived.
    (written with support from AVD)

    Arrival time is determined based on the first time step in the cube,
    not the release time and a new cube is created where the data points
    are hours (or fractions of hours) between the reference time and the
    arrival time.

    .. note:: Currently this subroutine can only be used on cubes with
              three dimensions

    .. note:: This subroutine assumes that the time dimension in the
              cube is in hours and will fail with an error where this
              is not true

    First import name_data and config

    >>> import name_data
    >>> import config

    Use the sample data path to locate data

    >>> sample_data_path = config.SAMPLE_DATADIR+'name/'

    Read in the data

    >>> name = name_data.NAMEData()
    >>> name.readdata(sample_data_path + 'Fields_grid1*',
    ... field_attributes = {'Species': 'CAESIUM-137'})
    [<iris 'Cube' of CAESIUM-137_AIR_CONCENTRATION / (Bq / m^3) \
(time: 9; latitude: 90; longitude: 180)>]

    Set a threshold

    >>> threshold = 1.0e-7

    Compute arrival times

    >>> time_cube = plume_arrival_time(name.gridded_cube_list[0], threshold)
    >>> print(time_cube.data[67, 98])
    96.0

    '''

    shape = cube.data.shape
    if len(shape) != 3:
        raise ValueError("Cube must have 3-dimensions")

    # Locate time coordinate
    try:
        cube.coords('time')
        t_coord = cube.coord('time')
        t_dim = cube.coord_dims('time')[0]
    except:
        raise ValueError('cube does not have a time coordinate')

    # Check that time coordinate units contains hours
    if 'hour' not in str(t_coord.units):
        raise ValueError('time coordinate units need to be hours')

    # Determine location/ name of remaining coordinates
    dims_list = []
    for nind, coord in enumerate(cube.coords()):
        if coord.name() != 'time':
            dims_list.append(nind)

    # Create a masked array
    time_array = ma.array(np.zeros(shape[1:]))
    time_array[:] = ma.masked

    # Loop over values which are greater than the threshold
    coordinates = np.where(cube.data >= threshold)
    ts = coordinates[t_dim]
    coord1 = coordinates[dims_list[0]]
    coord2 = coordinates[dims_list[1]]

    # Extract the first time step
    first_timestep = t_coord.points[0]

    # In the new array add data where the threshold is exceeded.
    for tind, c1, c2 in zip(ts, coord1, coord2):
        if (time_array[c1, c2] is ma.masked) or \
                (time_array[c1, c2] > t_coord.points[tind] - first_timestep):
            time_array[c1, c2] = t_coord.points[tind] - first_timestep

    # Add coordinates and units and create cube
    threshold_coord = iris.coords.AuxCoord(threshold,
                                           long_name='Threshold',
                                           units=cube.units)
    a_coord = cube.coords()[dims_list[0]]
    b_coord = cube.coords()[dims_list[1]]

    start_time = t_coord.units.num2date(first_timestep).isoformat(' ')

    if 'Species' in cube.attributes:
        long_name = 'Time of arrival of {}'.format(cube.attributes['Species'])
    else:
        long_name = 'Time of arrival'

    time_cube = iris.cube.Cube(time_array,
                               long_name=long_name,
                               units='hours since '+start_time,
                               attributes=cube.attributes,
                               dim_coords_and_dims=[(a_coord, 0),
                                                    (b_coord, 1)],
                               aux_coords_and_dims=[(threshold_coord, None)])

    # Add in source term coordinates if present
    names = [coord.name() for coord in cube.coords()]
    if 'source_latitude' in names:
        new_coord = cube.coord('source_latitude')
        time_cube.add_aux_coord(new_coord)
        new_coord = cube.coord('source_longitude')
        time_cube.add_aux_coord(new_coord)

    # Add short name
    time_cube.attributes['short_name'] = 'TimeOfArrival'

    return time_cube


if __name__ == '__main__':

    import doctest
    doctest.testmod()
