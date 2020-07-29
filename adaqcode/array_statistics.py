"""
Functions that calculate statistics on numpy arrays (not cubes).
"""

from six.moves.builtins import range
import iris
import numpy as np


def calc_forecast_day(forecast_period, runtime, day_start_hour=None):
    """
    Basic function to return the forecast day given the forecast_period
    (leadtime) and model runtime. Note Day 1 is the first full day of a
    forecast.

    :param forecast_period: Forecast range(s) for this time (leadtime, number of
                           hours since forecast started)
    :param runtime: Hour of day at which forecast run starts
                    (usually 18, 0, or 12)
    :param day_start_hour: Hour of the day to refer to as the start of a day.
                           By default this is set to 1 (1Z, so a day is
                           01Z - 24Z), or 12 (for 12Z runtimes,
                           so a day is 13Z - 12Z).

    For example a forecast_period of 7 hours from an 18Z model run (valid
    therefore at 01Z):

    >>> print(calc_forecast_day(7,18))
    1

    Or a forecast_period of 24 hours from a 00Z model run (valid therefore
    at midnight):

    >>> print(calc_forecast_day(24,0))
    1

    Or a forecast_period of 12 hours from a 12Z model run (valid at midnight):

    >>> print(calc_forecast_day(12,12))
    1

    Can also pass in a list of forecast periods to convert:

    >>> print(calc_forecast_day([3,7,24,36], 18))
    [0, 1, 1, 2]

    """

    #Set up hour of day to refer to as the start of a day
    #Note for hourly means, the convention is for the coordinate point to be
    #placed at the end of the meaning period, hence by default day_start_hour
    #is set to 1, which is equivalent to 0Z-1Z mean period.
    #However if runtime is 12Z, then for historical reasons, the first 'Day'
    #is counted as 12-13Z -> 11-12Z.
    #Note these defaults can be overridden with the keyword day_start_hour.
    if day_start_hour is None:
        day_start_hour = 1
        if runtime == 12:
            day_start_hour = 13

    #Number of hours from runtime to first day_start_hour:
    nhrs_to_start_of_day = (day_start_hour - runtime)%24
    #Number of hours after first day_start_hour:
    nhrs_after_start_of_first_day = np.array(forecast_period) - \
                                    np.array(nhrs_to_start_of_day)
    #Number of days after first day_start_hour:
    ndays_after_start_of_first_day = nhrs_after_start_of_first_day // 24
    #But day_start_hour should be referred to as 'Day1', so add an extra day on
    forecast_day = ndays_after_start_of_first_day + 1

    #Convert back to same type as input:
    if isinstance(forecast_period, int):
        forecast_day = int(forecast_day)
    elif isinstance(forecast_period, list):
        forecast_day = list(forecast_day)
    elif isinstance(forecast_day, np.ndarray):
        #Ensure all integers
        forecast_day = forecast_day.astype(int)

    return forecast_day


def calc_forecast_periods(forecast_day, runtime, day_start_hour=None):
    """
    Basic function to return a list of forecast_periods (leadtimes) for each
    forecast day given the model runtime. Note Day 1 is the first full day of a
    forecast.

    :param forecast_day: Forecast day required
    :param runtime: Hour of day at which forecast run starts
                    (usually 18, 0, or 12)
    :param day_start_hour: Hour of the day to refer to as the start of a day.
                           By default this is set to 1 (1Z, so a day is
                           01Z - 24Z), or 12 (for 12Z runtimes,
                           so a day is 13Z - 12Z).

    For example the forecast ranges for Day 1 from an 18Z model:

    >>> print(calc_forecast_periods(1, 18))
    [7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, \
26, 27, 28, 29, 30]

    """
    #Set up hour of day to refer to as the start of a day
    if day_start_hour is None:
        day_start_hour = 1
        if runtime == 12:
            day_start_hour = 13

    #Number of days since first day_start_hour:
    ndays_after_start_of_first_day = int(forecast_day) - 1
    #Number of hours since first day_start_hour:
    nhrs_after_start_of_first_day = ndays_after_start_of_first_day*24
    #Number of hours from runtime to first day_start_hour:
    nhrs_to_start_of_day = (day_start_hour - runtime)%24
    #Number of hours since runtime at start of forecast_day:
    nhrs_to_start_of_fcst_day = nhrs_to_start_of_day + \
                                nhrs_after_start_of_first_day
    #A day is 24 hours:
    forecast_periods = list(range(nhrs_to_start_of_fcst_day,
                                  nhrs_to_start_of_fcst_day+24))

    return forecast_periods

def nanstat_minperiod(array, stat='max', min_periods=1, axis=None):
    """
    Calculates a statistic, ignoring any nan values.
    If there is less than the required min_periods points
    of data, then this maximum value is set to nan.

    :param array: numpy array to calculate statistic on
    :param stat: string, statistic to calculate. Available options are
                 'max', 'mean', 'min', 'sum', 'std' (standard deviation),
                 'var' (variance)
    :param min_periods: minimum number of non-nan points required in
                        a period to calculate a non-nan answer
    :param axis: axis to calculate statistic over.

    >>> data = np.arange(10.)
    >>> data[4:8] = np.nan
    >>> np.set_printoptions(formatter={'float':lambda x: '{:5.0f}'.format(x)})
    >>> print(data)
    [    0     1     2     3   nan   nan   nan   nan     8     9]
    >>> np.set_printoptions()
    >>> print(nanstat_minperiod(data, stat='max'))
    9.0
    >>> print(nanstat_minperiod(data, stat='max', min_periods=7))
    nan

    Calculate mean instead:

    >>> print('{:.3f}'.format(nanstat_minperiod(data, stat='mean')))
    3.833

    Also works on multiple-axis arrays:

    >>> data2d = data.reshape(2,5)
    >>> np.set_printoptions(formatter={'float':lambda x: '{:5.0f}'.format(x)})
    >>> print(data2d)
    [[    0     1     2     3   nan]
     [  nan   nan   nan     8     9]]
    >>> print(nanstat_minperiod(data2d, stat='max', min_periods=2))
    9.0
    >>> print(nanstat_minperiod(data2d, stat='max', min_periods=2, axis=0))
    [  nan   nan   nan     8   nan]
    >>> print(nanstat_minperiod(data2d, stat='max', min_periods=2, axis=1))
    [    3     9]
    >>> np.set_printoptions()
    """

    #Calculate statistic for entire array
    if stat == 'max':
        output_array = np.nanmax(array, axis=axis)
    elif stat == 'mean':
        output_array = np.nanmean(array, axis=axis)
    elif stat == 'min':
        output_array = np.nanmin(array, axis=axis)
    elif stat == 'sum':
        output_array = np.nansum(array, axis=axis)
    elif stat == 'std':
        output_array = np.nanstd(array, axis=axis)
    elif stat == 'var':
        output_array = np.nanvar(array, axis=axis)
    else:
        raise ValueError("Unknown input value for stat " + stat)

    #Convert to floats - int(np.nan) is invalid
    output_array = output_array.astype(float)
    if array.ndim > 1:
        #multi-dimensional
        indices = np.where(~(np.nansum(~np.isnan(array), axis=axis)
                             >= min_periods))
        if indices[0].size:
            output_array[indices] = np.float64(np.nan)
            #note: nan must be an np.float64, not a standard float
            #this will allow it to have a dtype
            #otherwise, fails when used in aggregation if set to
            #standard np.nan
    else:
        if np.nansum(~np.isnan(array)) < min_periods:
            output_array = np.float64(np.nan)

    return output_array


def rolling_window_extended(array, window=1, step=1, axis=0):
    """
    Similar to `iris.util.rolling_window`, but extended at the
    beginning to add nan values where data is missing.

    :param array: numpy array containing input data to which rolling
                  window will be added to
    :param window: integer, size of rolling window
    :param step: integer, size of step between rolling windows
    :param axis: integer, axis to take rolling window over

    :returns: rolling window - an array that is a view of the original
              array, with an added dimension of the size of the given
              window at axis+1.

    For example:

    >>> d = np.arange(5.)
    >>> rwi = iris.util.rolling_window(d, window=3)
    >>> np.set_printoptions(formatter={'float':lambda x: '{:5.0f}'.format(x)})
    >>> print(rwi)
    [[    0     1     2]
     [    1     2     3]
     [    2     3     4]]
    >>> rwe = rolling_window_extended(d, window=3)
    >>> print(rwe)
    [[  nan   nan     0]
     [  nan     0     1]
     [    0     1     2]
     [    1     2     3]
     [    2     3     4]]
    >>> np.set_printoptions()

    Also works on arrays with multiple axes:

    >>> d = np.arange(30.).reshape(5,6)
    >>> rwe = rolling_window_extended(d,window=3, axis=1)
    >>> print(rwe.shape)
    (5, 6, 3)

    """
    #Based on iris.util.rolling_window

    #Check not integers - convert to floats if needed
    #(otherwise can't add NaNs in)
    if array.dtype == "int64" or array.dtype == "int32":
        array = array.astype(float)

    #Calculate standard rolling window - only starts
    # once including first window values
    rw = iris.util.rolling_window(array, window, step, axis=axis)

    #Now extend to include first parts of windows
    for win in np.arange(window-1, 0, -1):
        #Get the first array from the rolling window
        rw_extra = iris.util.rolling_window(array, win, step, axis=axis)
        rw_extra0 = np.take(rw_extra, 0, axis=axis)
        #And get a shortened version up to fill in gaps with NaNs
        rw_nan = iris.util.rolling_window(array, window-win, step, axis=axis)
        rw_nan0 = np.take(rw_nan, 0, axis=axis)
        rw_nan0[:] = np.nan
        #Add these two arrays together to get an array the same shape
        # as the required first array from the main rolling window
        values = np.concatenate((rw_nan0, rw_extra0), axis=axis)
        #Add this arry to the beginning of the main rolling window
        rw = np.insert(rw, 0, values, axis=axis)

    return rw


if __name__ == '__main__':

    import doctest
    doctest.testmod()
