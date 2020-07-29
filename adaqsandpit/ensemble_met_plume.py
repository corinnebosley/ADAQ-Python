"""
ensemble_met_plume.py

Experimental code for plotting ensemble meteorological data
The code plots time series of the ensemble mean with
shaded regions showing the 25-75 percentile range and the
10-90 percentile range. Can be used on any ensemble output of
met data from NAME. Currently meteorological variable and
file name are hardwired.

Susan, 19 April 2017
"""

import iris
import iris.coords as icoords
import matplotlib.pyplot as plt
import iris.plot as iplt
import pandas
import datetime
import matplotlib.dates as mpldates


_MET_SHORT_NAME = {
    '_TEMPERATURE_(C)': 'temperature',
    '_RELATIVE_HUMIDITY_(%)': 'relative_humidity',
    '_PRECIPITATION_RATE_(MM/HR)': 'precipitation',
    '_CLOUD_AMOUNT_(OKTAS)': 'cloud',
    '_WIND_DIRECTION_(DEGREES)': 'wind_direction',
    '_WIND_SPEED': 'wind_speed',
    '_BOUNDARY_LAYER_DEPTH': 'boundary_layer_depth',
    '_PASQUILL_STABILITY': 'pasquill_stability',
    '_PRESSURE_(PA)': 'pressure',
    '_POTENTIAL_TEMPERATURE_(K)': 'potential_temperature',
    '_SEA_LEVEL_PRESSURE_(PA)': 'sea_level_pressure'
    }

def set_xaxis_date_fmt(ax):
    """
    (Borrowed from line_plot)
    Set the dates on an x-axis nicely, spread over 2 lines
    If x-axis does not include dates, then has no effect.

    :param ax: matplotlib axis instance

    """

    xax = ax.get_xaxis() # get the x-axis

    #Get hold of existing locator and formatter objects
    major_locator = xax.get_major_locator()
    major_formatter = xax.get_major_formatter()

    if not isinstance(major_formatter,
                      pandas.tseries.converter.PandasAutoDateFormatter):
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

def my_callback(cube, field, filename):

    """
    Function to:
     * Add an ensemble coordinate
     * Rename cubes
     * Remove unwanted coordinates
    """

    unwanted_keys = ['Number of field cols',
                     'Number of preliminary cols',
                     'Run time',
                     'Met data']
    for key in unwanted_keys:
        if key in cube.attributes:
            del cube.attributes[key]

    if cube.long_name in _MET_SHORT_NAME:
        cube.attributes['short_name'] = _MET_SHORT_NAME[cube.long_name]
        cube.rename(cube.attributes['short_name'])

    if not cube.coords('realization'):
        ensemble_number = filename.strip('.txt').split('_')[-1]
        realization = ensemble_number[1:]
        ensemble_coord = icoords.AuxCoord(realization,
                                          standard_name='realization')
        cube.add_aux_coord(ensemble_coord)


def plot_met_plume(cube):
    '''
    Main plotting code.
    '''

    mean_cube = cube.collapsed(['realization'], iris.analysis.MEAN)
    iqr_cube = cube.collapsed(['realization'],
                              iris.analysis.PERCENTILE,
                              percent=[0, 10, 25, 75, 90, 100])
    t_coord = iqr_cube.coord('time')
    times = t_coord.units.num2date(t_coord.points)

    iplt.plot(mean_cube, 'k', linewidth=2)
    plt.fill_between(times, iqr_cube.data[0, :], iqr_cube.data[1, :],
                     edgecolor='none', color='#cce6ff')
    plt.fill_between(times, iqr_cube.data[1, :], iqr_cube.data[2, :],
                     edgecolor='none', color='#66b3ff')
    plt.fill_between(times, iqr_cube.data[2, :], iqr_cube.data[3, :],
                     edgecolor='none', color='#0080ff')
    plt.fill_between(times, iqr_cube.data[3, :], iqr_cube.data[4, :],
                     edgecolor='none', color='#66b3ff')
    plt.fill_between(times, iqr_cube.data[4, :], iqr_cube.data[5, :],
                     edgecolor='none', color='#cce6ff')
    ax = plt.gca()
    ax = set_xaxis_date_fmt(ax)
    ylabel = '{} [{}]'.format(cube.name().replace('_', ' ').title(),
                              cube.units)

    plt.ylabel(ylabel)

    plt.show()

def main():
    '''
    Main code. Change workdir to plot a different dataset and
    wind_speed to plot a different variable.
    '''

    workdir = '/data/users/apdg/python_sample_data/name_ensemble/'
    cube = iris.load_cube(workdir + 'Met_Data_*.txt',
                          'wind_speed',
                          callback=my_callback)
    plot_met_plume(cube)

if __name__ == '__main__':
    main()
