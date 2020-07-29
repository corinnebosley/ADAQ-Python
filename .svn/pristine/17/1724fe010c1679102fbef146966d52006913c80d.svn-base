# Crown copyright 2014
"""
Functions relating to meteorology
"""
from __future__ import division
from __future__ import print_function
import warnings
import iris
import numpy as np
import cf_units
import constants

def relative_humidity(cubelist):
    """
    Generates cube of relative humidity. Expression for water saturation vapour
    pressure comes from UMDP029-The Large Scale Cloud Scheme and Saturated
    Specific Humidity.

    :param cubelist: cubelist containing cubes of pressure, temperature and
               specific humidity.
    :return cubelist: cubelist containing same cubes as input list but with a
                relative humidity cube added.
    """

    #Check input list contains required cubes and extract
    try:
        pressure = cubelist.extract('air_pressure', strict=True)
    except:
        raise ValueError('Missing pressure data required for relative \
                             humidity generation')
    try:
        temperature = cubelist.extract('air_temperature', strict=True)
    except:
        raise ValueError('Missing temperature data required for relative \
                             humidity generation')
    try:
        spc_humid = cubelist.extract('specific_humidity', strict=True)
    except:
        raise ValueError('Missing specific_humidity data required for \
                             relative humidity generation')


    #Make relative humidity cube
    rel_humid = spc_humid.copy()
    rel_humid.rename('relative_humidity')
    rel_humid.units = '%'
    rel_humid.attributes['short_name'] = 'rh'

    #See UMDP reference in docstring for source of these numbers used in the
    #calculation of the water saturation vapour pressure
    a1 = 10.79574
    a2 = -5.028
    a3 = 1.50475E-04
    a4 = 0.42873E-03
    temp_trp = 273.16
    e1 = -8.2969
    e2 = 4.76955

    t1 = a1*(1.0-temp_trp/temperature.data)
    t2 = a2*np.log10(temperature.data/temp_trp)
    t3a = 10.0**(e1*((temperature.data/temp_trp) -1))
    t3 = a3*(1.0-t3a)
    t4 = a4*(10.0**(e2*(1.0-temp_trp/temperature.data)))
    t5 = 0.78614
    t6 = 2.0

    p_wv_sat = 10.0**(t1+t2+t3+t4+t5+t6)

    p_wv = pressure.data*spc_humid.data/(spc_humid.data*(1.0- \
                                constants.MMASS_RATIO)+constants.MMASS_RATIO)
    rel_humid.data = 100.0*p_wv/p_wv_sat

    cubelist.append(rel_humid)

    return cubelist


def relative_humidity_list(md_list):
    """
    Modifies a cubelist of ADAQData objects by adding a relative humidity
    cube to each list item.     This function takes as input a list
    of ADAQData objects which must contain cubes of pressure,
    temperature and specific humidity.

    :param md_list: list of model data ADAQData objects in the form of
               :class:`adaq_data.ADAQData` objects.

    :return md_list: list of model data ADAQData objects with relative humidity
             cube added to both sites and gridded cubelists.

    Get some sample data

    >>> import adaq_functions
    >>> ini_dict, sites_data, od, md_list = adaq_functions.get_exampledata(
    ...    sites_cube_list=True, exampletype='3d',
    ...    gridded_cube_list=True) # doctest: +ELLIPSIS
    Reading inifile /.../example_data_3d.ini
    Number of sites:  5
    >>> print(md_list[0]) # doctest: +ELLIPSIS
    <class '....ADAQData'> - Subclass of ADAQData - Contains:
    sites_cube_list:
    0: air_pressure / (Pa)                 (site_id: 5; time: 2; \
model_level_number: 38)
    1: air_temperature / (K)               (site_id: 5; time: 2; \
model_level_number: 38)
    2: mass_fraction_of_ozone_in_air / (kg kg-1) (site_id: 5; time: 2; \
model_level_number: 38)
    3: specific_humidity / (kg kg-1)       (site_id: 5; time: 2; \
model_level_number: 38)
    gridded_cube_list:
    0: air_pressure / (Pa)                 (time: 2; model_level_number: 38; \
grid_latitude: 182; grid_longitude: 146)
    1: air_temperature / (K)               (time: 2; model_level_number: 38; \
grid_latitude: 182; grid_longitude: 146)
    2: mass_fraction_of_ozone_in_air / (kg kg-1) (time: 2; \
model_level_number: 38; grid_latitude: 182; grid_longitude: 146)
    3: specific_humidity / (kg kg-1)       (time: 2; model_level_number: 38; \
grid_latitude: 182; grid_longitude: 146)
    trajectory_cube_list:
    < No cubes >

    Call routine to generate and add relative humidity cube

    >>> md_list = relative_humidity_list(md_list)
    >>> md = md_list[0]
    >>> print(md)  # doctest: +ELLIPSIS
    <class '....ADAQData'> - Subclass of ADAQData - Contains:
    sites_cube_list:
    0: air_pressure / (Pa)                 (site_id: 5; time: 2; \
model_level_number: 38)
    1: air_temperature / (K)               (site_id: 5; time: 2; \
model_level_number: 38)
    2: mass_fraction_of_ozone_in_air / (kg kg-1) (site_id: 5; time: 2; \
model_level_number: 38)
    3: specific_humidity / (kg kg-1)       (site_id: 5; time: 2; \
model_level_number: 38)
    4: relative_humidity / (%)             (site_id: 5; time: 2; \
model_level_number: 38)
    gridded_cube_list:
    0: air_pressure / (Pa)                 (time: 2; \
model_level_number: 38; grid_latitude: 182; grid_longitude: 146)
    1: air_temperature / (K)               (time: 2; \
model_level_number: 38; grid_latitude: 182; grid_longitude: 146)
    2: mass_fraction_of_ozone_in_air / (kg kg-1) (time: 2; \
model_level_number: 38; grid_latitude: 182; grid_longitude: 146)
    3: specific_humidity / (kg kg-1)       (time: 2; \
model_level_number: 38; grid_latitude: 182; grid_longitude: 146)
    4: relative_humidity / (%)             (time: 2; \
model_level_number: 38; grid_latitude: 182; grid_longitude: 146)
    trajectory_cube_list:
    < No cubes >

    Extract a rel humid cube and check mean value

    >>> rh = md.extract(short_name='rh', singlecube=True, gridded=True)
    >>> mean_rh_field = rh.data[0,0,:,:].mean()
    >>> print('Mean rh field value = {:.3f} {}'.format(mean_rh_field, rh.units))
    Mean rh field value = 80.107 %

    """
    for md in md_list:
        if md is not None:
            if md.sites_cube_list:
                md.sites_cube_list = relative_humidity(md.sites_cube_list)
            if md.gridded_cube_list:
                md.gridded_cube_list = relative_humidity(md.gridded_cube_list)

    return md_list

def wind_speed(cubelist, short_name='ws'):
    """
    Calculate wind speed from x and y wind components (u and v)
    which are cubes within cubelist.
    Appends new cube to existing cubelist.

    :param cubelist: iris cubelist, must contain cubes named 'x_wind' and
                     'y_wind'
    :param short_name: required short_name to be added to output cube
    :returns: cubelist with extra wind speed cube included

    >>> import config
    >>> sample_datadir = config.SAMPLE_DATADIR+'/aqum_output/oper/met/'
    >>> cubelist = iris.load(sample_datadir+'prods_op_aqum_20170602_18.000.pp',
    ... ['x_wind','y_wind'])
    >>> print(cubelist) # doctest: +NORMALIZE_WHITESPACE
    0: x_wind / (m s-1)  (time: 6; grid_latitude: 183; grid_longitude: 146)
    1: y_wind / (m s-1)  (time: 6; grid_latitude: 183; grid_longitude: 146)

    >>> cubelist = wind_speed(cubelist)
    >>> print(cubelist)
    0: x_wind / (m s-1)                    (time: 6; grid_latitude: 183; \
grid_longitude: 146)
    1: y_wind / (m s-1)                    (time: 6; grid_latitude: 183; \
grid_longitude: 146)
    2: wind_speed / (m s-1)                (time: 6; grid_latitude: 183; \
grid_longitude: 146)

    >>> ws = cubelist.extract('wind_speed', strict=True)
    >>> print(ws) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
    wind_speed / (m s-1)                (time: 6; grid_latitude: 183; \
grid_longitude: 146)
         Dimension coordinates:
              time                      x                 -                    -
              grid_latitude             -                 x                    -
              grid_longitude            -                 -                    x
         Auxiliary coordinates:
              forecast_period           x                 -                    -
         Scalar coordinates:
              forecast_reference_time: 2017-06-02 18:00:00
              height: 10.0 m
         Attributes:
              short_name: ws
    """

    try:
        x_wind = cubelist.extract('x_wind', strict=True)
    except iris.exceptions.ConstraintMismatchError:
        warnings.warn('wind speed not calculated: x_wind missing')
        return cubelist
    try:
        y_wind = cubelist.extract('y_wind', strict=True)
    except iris.exceptions.ConstraintMismatchError:
        warnings.warn('wind speed not calculated: y_wind missing')
        return cubelist

    #Check they are on same grid (U&V often output on different grids)
    assert x_wind.coord('grid_longitude') == y_wind.coord('grid_longitude')
    assert x_wind.coord('grid_latitude') == y_wind.coord('grid_latitude')

    ws = (x_wind**2 + y_wind**2)**0.5
    ws.convert_units('m s-1')
    ws.rename('wind_speed')
    if 'label' in x_wind.attributes:
        ws.attributes['label'] = x_wind.attributes['label']
    ws.attributes['short_name'] = short_name
    cubelist.append(ws)

    return cubelist

def wind_direction(cubelist, short_name='wind_dir'):
    """
    Calculate the direction wind is going to,
    using the x and y wind components (u and v)
    which are within the cubelist.
    Also rotates the winds to a regular lat-long grid if originally on a rotated grid.
    Appends new cube to existing cubelist

    :param cubelist: iris cubelist, must contain cubes named
                     'x_wind' and'y_wind'
    :param short_name: required short_name to be added to output cube
    :returns: cubelist with extra wind direction cube included

    >>> import config
    >>> sample_datadir = config.SAMPLE_DATADIR+'/aqum_output/oper/met/'
    >>> cubelist = iris.load(sample_datadir+'prods_op_aqum_20170602_18.000.pp',
    ... ['x_wind','y_wind'])
    >>> print(cubelist) # doctest: +NORMALIZE_WHITESPACE
    0: x_wind / (m s-1)  (time: 6; grid_latitude: 183; grid_longitude: 146)
    1: y_wind / (m s-1)  (time: 6; grid_latitude: 183; grid_longitude: 146)

    >>> cubelist = wind_direction(cubelist)
    >>> print(cubelist)
    0: x_wind / (m s-1)                    (time: 6; grid_latitude: 183; grid_longitude: 146)
    1: y_wind / (m s-1)                    (time: 6; grid_latitude: 183; grid_longitude: 146)
    2: wind_to_direction / (degrees)       (time: 6; grid_latitude: 183; grid_longitude: 146)

    >>> wind_dir = cubelist.extract('wind_to_direction', strict=True)
    >>> print(wind_dir)
    wind_to_direction / (degrees)       (time: 6; grid_latitude: 183; grid_longitude: 146)
         Dimension coordinates:
              time                           x                 -                    -
              grid_latitude                  -                 x                    -
              grid_longitude                 -                 -                    x
         Auxiliary coordinates:
              forecast_period                x                 -                    -
         Scalar coordinates:
              forecast_reference_time: 2017-06-02 18:00:00
              height: 10.0 m
         Attributes:
              STASH: m01s03i225
              short_name: wind_dir
              source: Data from Met Office Unified Model
              um_version: 10.8

    >>> x_wind = cubelist.extract('x_wind', strict=True)
    >>> y_wind = cubelist.extract('y_wind', strict=True)
    >>> print(x_wind.data[0,0,0], y_wind.data[0,0,0], wind_dir.data[0,0,0])
    2.25 -3.875 142.07626
    """

    #Getting the x and y wind components
    try:
        x_wind = cubelist.extract('x_wind', strict=True)
    except:
        raise ValueError('wind direction not calculated, x_wind missing')
    try:
        y_wind = cubelist.extract('y_wind', strict=True)
    except:
        raise ValueError('wind direction not calculated, y_wind missing')

    #Checking that the wind components are on the same grid
    assert x_wind.coord('grid_longitude') == y_wind.coord('grid_longitude')
    assert x_wind.coord('grid_latitude') == y_wind.coord('grid_latitude')
    wind_dir = x_wind.copy()

    #Rotating the wind compnonets, from the rotated model, to real world
    lat_lon_coord_system = iris.coord_systems.GeogCS(
        semi_major_axis=iris.fileformats.pp.EARTH_RADIUS)
    u_wind, v_wind = iris.analysis.cartography.rotate_winds(x_wind, y_wind, lat_lon_coord_system)

    #Getting the wind direction
    wind_dir.data = np.rad2deg(np.arctan2(u_wind.data, v_wind.data)) % 360

    #Changing units and names
    wind_dir.rename('wind_to_direction')
    wind_dir.attributes['short_name'] = short_name
    wind_dir.units = cf_units.Unit('degrees')
    cubelist.append(wind_dir)

    return cubelist

if __name__ == "__main__":

    import doctest
    doctest.testmod()
