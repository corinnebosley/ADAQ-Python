#pylint: disable=line-too-long
"""
Module to hold functions, constants and colour bars related to air quality indices.
Currently included:

* **DAQI**: UK Daily Air Quality Index.
  http://uk-air.defra.gov.uk/assets/documents/reports/cat14/1304251155_Update_on_Implementation_of_the_DAQI_April_2013_Final.pdf
* **CAMS**: Colour scale and index levels used for
  Copernicus Atmospheric Monitoring Service. Eg see
  http://macc-raq.copernicus-atmosphere.eu/index.php?category=forecasts

Note when this module is loaded, extra colour bars are registered,
so they can be called directly by their names:

>>> import matplotlib.pyplot as plt
>>> cmap = plt.get_cmap('DAQI')
>>> print(cmap) # doctest: +ELLIPSIS
<matplotlib.colors.ListedColormap object at ...>
>>> print(cmap.name)
DAQI
>>> cmap = plt.get_cmap('CAMS')
>>> print(cmap) # doctest: +ELLIPSIS
<matplotlib.colors.ListedColormap object at ...>

"""
#pylint: disable=line-too-long

from __future__ import division
from __future__ import absolute_import
from __future__ import print_function

import warnings
import datetime

import numpy as np
from matplotlib import colors, cm
import iris

import cube_statistics
import chemistry

# List of species that map to PM10 index levels
_PM10_SPECIES = ['PM10', 'PM10_NH42SO4', 'PM10_BC', 'PM10_BB', 'PM10_OCFF',
                 'PM10_DUST', 'PM10_NH4NO3', 'PM10_SOA']
# List of species that map to PM2p5 index levels
_PM2p5_SPECIES = ['PM2p5', 'PM2p5_NH42SO4', 'PM2p5_BC', 'PM2p5_BB', 'PM2p5_OCFF',
                  'PM2p5_DUST', 'PM2p5_NH4NO3', 'PM2p5_SOA']


#-------- DAQI ----------

#Dictionary of Daily Air Quality Index levels in ug/m3.
#This is a private method to contain the levels only once
#- DAQI_LEVELS should be used to enable checking against multiple short_names
# Note, pylint is disabled for this constant to make it much easier
# to match levels to comments
# pylint: disable=line-too-long
#                 [   1  |   2  |  3   |  4   |  5   |  6   |  7   |  8   | 9   | 10-     ]
#                 [       Low          |    Moderate        |        High       |VeryHigh ]
# pylint: disable=bad-whitespace
_DAQI_LEVS = \
    {
        'DAQI':
        np.array(([1,      2,     3,     4,     5,     6,     7,     8,     9,     10,     11])),
        'O3':
        np.array(([0.,  33.5,  66.5, 100.5, 120.5, 140.5, 160.5, 187.5, 213.5,  240.5,  500.0])),
        'NO2':
        np.array(([0.,  67.5, 134.5, 200.5, 267.5, 334.5, 400.5, 467.5, 534.5,  600.5, 1000.0])),
        'SO2':
        np.array(([0.,  88.5, 177.5, 266.5, 354.5, 443.5, 532.5, 710.5, 887.5, 1064.5, 2000.0])),
        'PM2p5':
        np.array(([0.,  11.5,  23.5,  35.5,  41.5,  47.5,  53.5,  58.5,  64.5,   70.5,  140.0])),
        'PM10':
        np.array(([0.,  16.5,  33.5,  50.5,  58.5,  66.5,  75.5,  83.5,  91.5,  100.5,  200.0]))
    }
# pylint: enable=line-too-long
# pylint: enable=bad-whitespace


#: Dictionary of Daily Air Quality Index (DAQI) levels, with a key of short_name
#: in ug/m3, as defined by
# pylint: disable=line-too-long
#: http://uk-air.defra.gov.uk/assets/documents/reports/cat14/1304251155_Update_on_Implementation_of_the_DAQI_April_2013_Final.pdf
# pylint: enable=line-too-long
#:
#: DAQI_LEVELS['short_name'] = Numpy array of maximum values, such that:
#:  * array[i] gives the maximum value for DAQI i.
#:  * array[0] = 0. which is a sensible min value for contours
#:  * array[10] is a sensible max value for contours, approx 2x DAQI-10
#:
#: For short_names not directly included in the DAQI,
#: levels are also set, but such that they correspond to DAQI levels
#: for appropriate species, for example NO = 10% of NO2 levels
#: and PM10_DUST = PM10 levels.
#:
#: .. note:: This assumes all data is in ug/m3 as per definition of the DAQI.
#:           No explict checking of units are done here.
#:
DAQI_LEVELS = {'DAQI' : _DAQI_LEVS['DAQI'],
               'O3' : _DAQI_LEVS['O3'],
               'O3_DAQI' : _DAQI_LEVS['DAQI'],
               'NO2' : _DAQI_LEVS['NO2'],
               'NO2_DAQI' : _DAQI_LEVS['DAQI'],
               'NO' : _DAQI_LEVS['NO2'] * 0.1,
               'SO2' : _DAQI_LEVS['SO2'],
               'SO2_DAQI' : _DAQI_LEVS['DAQI'],
               # CO: Based on previous version of DAQI
               'CO' : 1000.*np.array(([0., 3.9, 7.7, 11.6, 13.5, 15.5,
                                       17.4, 19.3, 21.3, 23.2])),
               'CO_DAQI' : _DAQI_LEVS['DAQI']
              }
#Also add all the PM10/PM2p5 species
for species in _PM10_SPECIES:
    DAQI_LEVELS[species] = _DAQI_LEVS['PM10']
    DAQI_LEVELS[species+'_DAQI'] = _DAQI_LEVS['DAQI']
for species in _PM2p5_SPECIES:

    DAQI_LEVELS[species] = _DAQI_LEVS['PM2p5']
    DAQI_LEVELS[species+'_DAQI'] = _DAQI_LEVS['DAQI']


def daqi_applylevels(cube):
    """
    Apply DAQI levels to a cube.

    :param cube: input cube
    :returns: cube with value converted to DAQI-based levels if possible,
              otherwise None is returned.

    This routine will apply the DAQI level values to any cube that has a
    short_name in DAQI_LEVELS and with units of ug/m3
    (for any others, None is returned)

    Converts values in ug/m3 to their appropriate DAQI level (1-10).
    Renames cube to daily_air_quality_index_of_<molecule-name>
    and sets units to 1. It also appends _DAQI to the short name.
    Note this routine does NOT apply any time-meaning to the cube, it just
    converts whatever values of data it is given. The cube returned may
    therefore not be a true DAQI.

    >>> import config
    >>> samplepath = config.SAMPLE_DATADIR+'sites_cube_list/'
    >>> constraint = (iris.AttributeConstraint(short_name="O3") &
    ...     iris.Constraint(abbrev="YW"))
    >>> o3_ugm3 = iris.load_cube(samplepath+"aurn_1days.nc", constraint)
    >>> print(o3_ugm3.summary(shorten=True))
    mass_concentration_of_ozone_in_air / (ug/m3) (time: 25)

    >>> o3_daqilevs = daqi_applylevels(o3_ugm3)
    >>> print(o3_daqilevs.summary(shorten=True))
    daily_air_quality_index_of_ozone / (1) (time: 25)
    >>> print(o3_daqilevs.attributes['short_name'])
    O3_DAQI
    >>> print(o3_ugm3.data.max(), o3_daqilevs.data.max())
    59.0 2.0
    """

    short_name = cube.attributes['short_name']
    if short_name in DAQI_LEVELS and cube.units == 'ug/m3':
        daqi_levels = DAQI_LEVELS[short_name]
        daqi_cube = cube.copy()
        molecule_name = chemistry.get_molecule_name(cube.name())
        daqi_cube.rename('daily_air_quality_index_of_'+molecule_name)
        daqi_cube.attributes['short_name'] += '_DAQI'
        daqi_cube.units = '1'

        for ilev, level in enumerate(daqi_levels):
            #Loop through each (increasing) level and overwrite data values
            #with DAQI level if greater/equal to minimum required for this level
            indices = np.where(cube.data >= level)
            daqi_cube.data[indices] = ilev+1
    else:
        daqi_cube = None

    return daqi_cube

def calc_species_daqi(cube):
    """
    Calculate DAQI for a single cube.

    :param cube: Input cube
    :returns: cube, converted to DAQI if possible. Otherwise returns None.

    This routine does appropriate time-averaging to convert a cube to a daily
    cube and then applys the DAQI levels. This therefore results in the DAQI
    for this particular species.

    >>> import config
    >>> samplepath = config.SAMPLE_DATADIR+'sites_cube_list/'
    >>> constraint = (iris.AttributeConstraint(short_name="PM2p5") &
    ...     iris.Constraint(abbrev="HAR"))
    >>> pm_ugm3 = iris.load_cube(samplepath+"aurn_5days.nc", constraint)
    >>> print(pm_ugm3.summary(shorten=True))
    mass_concentration_of_pm2p5_ambient_aerosol_in_air / (ug/m3) (time: 121)
    >>> print(pm_ugm3.data.max())
    69.0

    >>> pm_daqi = calc_species_daqi(pm_ugm3)
    >>> print(pm_daqi) # doctest: +ELLIPSIS
    daily_air_quality_index_of_pm2p5_ambient_aerosol / (1) (time: 6)
         Dimension coordinates:
              time                                              x
         Auxiliary coordinates:
              date                                              x
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
              short_name: PM2p5_DAQI
              source: AURN
         Cell methods:
              mean: time (1 hour)
              nanmean_min18periods: date
    >>> np.set_printoptions(formatter={'float':lambda x: '{:5.0f}'.format(x)})
    >>> print(pm_daqi.data)
    [  nan     2     2     3     4     3]
    >>> np.set_printoptions()

    """

    assert isinstance(cube, iris.cube.Cube)

    #By default, return None
    #Note this might be because already converted to DAQI
    # (short_name would end _DAQI)
    #Or because it is something that can not be converted,
    # eg a met variable.
    species_daqi = None

    if cube.units == 'ug/m3':
        #Note can't convert to DAQI if not in correct units!

        daily_cube = None

        #Check working with hourly data (can't calculate DAQI otherwise)
        dtpts = cube.coord('time').units.num2date(cube.coord('time').points)
        delta = min(dtpts[1:] - dtpts[:-1])
        if delta == datetime.timedelta(hours=1):

            short_name = cube.attributes['short_name']
            if short_name in ['NO2', 'SO2', 'NO']:
                daily_cube = cube_statistics.daily_stat(
                    cube, stat='max', min_periods=18, aqdates=True)
            elif short_name in _PM10_SPECIES or short_name in _PM2p5_SPECIES:
                daily_cube = cube_statistics.daily_stat(
                    cube, stat='mean', min_periods=18, aqdates=True)
            elif short_name in ['O3', 'CO']:
                daily_cube = cube_statistics.maximum_daily_rolling_8hr_mean(
                    cube, aqdates=True)

            if daily_cube is not None:
                #Has matched one of the above species and converted to a daily
                #cube, therefore can apply DAQI levels.
                species_daqi = daqi_applylevels(daily_cube)

    return species_daqi


def calc_daqi(cubelist):
    """
    Given a list of cubes, calculate a total DAQI from these species.
    Note the list of cubes should contain O3, NO2, SO2, PM10, and PM2p5.
    If any of these are missing, then the DAQI will still be calculated
    (assuming at least one species available), but it will not be a
    true representation of the DAQI.

    :param cubelist: list of iris cubes. Should contain required species,
                     given by the short_names:

                     * 'O3' or 'O3_DAQI'
                     * 'NO2' or 'NO2_DAQI'
                     * 'SO2' or 'SO2_DAQI'
                     * 'PM10' or 'PM10_DAQI'
                     * 'PM2p5' or 'PM2p5_DAQI'

    :returns: cube - containing DAQI values in data, and correct
              name, short_name and units.

    >>> import config
    >>> samplepath = config.SAMPLE_DATADIR+'sites_cube_list/'
    >>> cubelist = iris.load(samplepath+"aqum_oper_5days.nc")

    Sort cubes alphabetically (only required for doctests)

    >>> cubelist = iris.cube.CubeList(sorted(list(
    ... cubelist), key=lambda cube: cube.name()))
    >>> print(cubelist)
    0: mass_concentration_of_nitrogen_dioxide_in_air / \
(ug/m3) (site_id: 5; time: 121)
    1: mass_concentration_of_ozone_in_air / (ug/m3) \
(site_id: 5; time: 121)
    2: mass_concentration_of_pm10_dry_aerosol_in_air / \
(ug/m3) (site_id: 5; time: 121)
    3: mass_concentration_of_pm2p5_dry_aerosol_in_air / \
(ug/m3) (site_id: 5; time: 121)
    4: mass_concentration_of_sulfur_dioxide_in_air / \
(ug/m3) (site_id: 5; time: 121)

    >>> daqi = calc_daqi(cubelist)
    >>> print(daqi) # doctest: +ELLIPSIS
    daily_air_quality_index / (1)       (site_id: 5; time: 6)
         Dimension coordinates:
              site_id                           x        -
              time                              -        x
         Auxiliary coordinates:
              abbrev                            x        -
              grid_latitude                     x        -
              grid_longitude                    x        -
              latitude                          x        -
              longitude                         x        -
              site_altitude                     x        -
              site_name                         x        -
              site_type                         x        -
              surface_altitude                  x        -
              date                              -        x
              forecast_period                   -        x
         Scalar coordinates:
              forecast_day: 1.0 Days
              level_height: 20.00... m, bound=(0.0, 49.99...) m
              model_level_number: 1
              sigma: 0.99..., bound=(1.0, 0.99...)
         Attributes:
              Conventions: CF-1.5
              STASH: m01s34i004
              label: aqum_oper
              short_name: DAQI
              source: Data from Met Office Unified Model
         Cell methods:
              mean: time (1 hour)
              nanmax_min18periods: date
    >>> np.set_printoptions(formatter={'float':lambda x: '{:5.0f}'.format(x)})
    >>> print(daqi.data)
    [[  nan     3     3     3     3     4]
     [  nan     3     3     3     3     2]
     [  nan     3     3     3     3     3]
     [  nan     3     3     3     3     3]
     [  nan     2     2     3     3     3]]
    >>> np.set_printoptions()

    """

    available_short_names = [cube.attributes['short_name']
                             for cube in cubelist]

    #Set up list to store individual species converted to DAQI in.
    daqi_cube_list = iris.cube.CubeList()

    #List of species that are required to calculate DAQI
    required_short_names = ['NO2', 'SO2', 'PM10', 'PM2p5', 'O3']

    #Loop through to find cube from input cube list for each of these species
    for short_name in required_short_names:

        if short_name + '_DAQI' in available_short_names:
            #Already converted to DAQI levels
            #Extract species and add to list of DAQI cubes
            species_daqi = cubelist.extract(
                iris.AttributeConstraint(short_name=short_name + '_DAQI'),
                strict=True)
            daqi_cube_list.append(species_daqi)

        elif short_name in available_short_names:
            #Basic cube available - extract this
            species_cube = cubelist.extract(
                iris.AttributeConstraint(short_name=short_name),
                strict=True)
            #Convert to DAQI levels
            species_daqi = calc_species_daqi(species_cube)
            if species_daqi is not None:
                #And add to list of DAQI cubes
                daqi_cube_list.append(species_daqi)

        else:
            #Required species not found
            warnings.warn(short_name +
                          ' not found - can not be used in DAQI calculation')

    if not daqi_cube_list:
        daqi = None
    else:
        #Calculate maximum across all species.

        #Pick cube with most times to act as the starting template
        nt = [len(cube.coord('time').points) for cube in daqi_cube_list]
        daqi_cube = daqi_cube_list[np.argmax(nt)]
        daqi_slices = iris.cube.CubeList([sl for sl in daqi_cube.slices_over('time')])

        for species_daqi in daqi_cube_list:
            for sp_slice in species_daqi.slices_over('time'):
                #Find matching daqi slice at the same time
                daqi_slice = None
                for d_slice in daqi_slices:
                    if d_slice.coord('date').points[0] == sp_slice.coord('date').points[0]:
                        daqi_slice = d_slice
                        break
                if daqi_slice is None: #No match found
                    print('species slice:', sp_slice)
                    print('date in species slice:', sp_slice.coord('date'))
                    print('dates in daqi cube:', daqi_cube.coord('dates'))
                    raise ValueError('Matching time not found for DAQI calculation')
                #fmax: If a nan exists for a value,
                #      then overwritten by other species if possible.
                daqi_slice.data = np.fmax(daqi_slice.data, sp_slice.data)

        #Merge individual slices back together
        daqi = daqi_slices.merge_cube()
        #Ensure site_id is first dim if available (rearranged in merging)
        if daqi.coords('site_id'):
            if daqi.coord_dims('site_id') != (0,):
                daqi.transpose()
        daqi.rename('daily_air_quality_index')
        daqi.attributes['short_name'] = 'DAQI'
        daqi.units = '1'

    return daqi

def convert_all_daqi(cubelist):
    """
    Convert all cubes to daqi if possible and also add daqi into cubelist
    For all cubes in cubelist, try and convert to daqi
    If they can be converted, daqi version added to cubelist,
    along with original cube

    :param cubelist: list of iris cubes
    :returns: new cubelist - list of iris cubes, includes input cubes,
              plus cubes converted to daqi where possible,
              plus total daqi cube.

    >>> import config
    >>> samplepath = config.SAMPLE_DATADIR+'sites_cube_list/'
    >>> cubelist = iris.load(samplepath+"aqum_oper_5days.nc")

    Sort cubes alphabetically (only required for doctests)

    >>> cubelist = iris.cube.CubeList(sorted(list(
    ... cubelist), key=lambda cube: cube.name()))
    >>> print(cubelist)
    0: mass_concentration_of_nitrogen_dioxide_in_air / (ug/m3) \
(site_id: 5; time: 121)
    1: mass_concentration_of_ozone_in_air / (ug/m3) \
(site_id: 5; time: 121)
    2: mass_concentration_of_pm10_dry_aerosol_in_air / (ug/m3) \
(site_id: 5; time: 121)
    3: mass_concentration_of_pm2p5_dry_aerosol_in_air / (ug/m3) \
(site_id: 5; time: 121)
    4: mass_concentration_of_sulfur_dioxide_in_air / (ug/m3) \
(site_id: 5; time: 121)

    >>> newcubelist = convert_all_daqi(cubelist)
    >>> print(newcubelist)
    0: mass_concentration_of_nitrogen_dioxide_in_air / (ug/m3) \
(site_id: 5; time: 121)
    1: mass_concentration_of_ozone_in_air / (ug/m3) \
(site_id: 5; time: 121)
    2: mass_concentration_of_pm10_dry_aerosol_in_air / (ug/m3) \
(site_id: 5; time: 121)
    3: mass_concentration_of_pm2p5_dry_aerosol_in_air / (ug/m3) \
(site_id: 5; time: 121)
    4: mass_concentration_of_sulfur_dioxide_in_air / (ug/m3) \
(site_id: 5; time: 121)
    5: daily_air_quality_index_of_nitrogen_dioxide / (1) \
(site_id: 5; time: 6)
    6: daily_air_quality_index_of_ozone / (1) \
(site_id: 5; time: 6)
    7: daily_air_quality_index_of_pm10_dry_aerosol / (1) \
(site_id: 5; time: 6)
    8: daily_air_quality_index_of_pm2p5_dry_aerosol / (1) \
(site_id: 5; time: 6)
    9: daily_air_quality_index_of_sulfur_dioxide / (1) \
(site_id: 5; time: 6)
    10: daily_air_quality_index / (1)       \
(site_id: 5; time: 6)

    """

    for cube in cubelist:
        #Convert individual species cubes
        species_daqi = calc_species_daqi(cube)
        if species_daqi is not None:
            cubelist.append(species_daqi)
    #Calculate overall DAQI
    daqi = calc_daqi(cubelist)
    if daqi is not None:
        #Add this DAQI to the end of the cubelist
        cubelist.append(daqi)

    return cubelist



def daqi_cmap():
    '''
    Generate colour map based on DAQI colours

    >>> cmap = daqi_cmap()

    To plot and display colour map:

    >>> import config
    >>> import plotting_functions
    >>> filename = config.CODE_DIR + "/adaqdocs/figures/daqi_cmap.png"
    >>> plotting_functions.plot_cmap(cmap, filename=filename)
    ... # doctest: +ELLIPSIS
    Saved figure .../daqi_cmap.png

    .. image:: ../adaqdocs/figures/daqi_cmap.png
       :scale: 50%
    '''
# pylint: disable=bad-whitespace
    colours = [[156, 255, 156],
               [ 49, 255,   0],
               [ 49, 207,   0],
               [255, 255,   0],
               [255, 207,   0],
               [255, 154,   0],
               [255, 100, 100],
               [255,   0,   0],
               [153,   0,   0],
               [206,  48, 255]]
# pylint: enable=bad-whitespace
    colours = np.array(colours)/255.
    cmap = colors.ListedColormap(colours, name='DAQI')

    #Register colour map, so can be called
    # eg plt.get_cmap('DAQI')
    cm.register_cmap(cmap=cmap)

    return cmap


#-------- CAMS ----------

# Contour levels for CAMS (Copernicus Atmospheric Monitoring Service)
# This is a private method to contain the levels only once
# - CAMS_LEVELS should be used to enable checking against multiple short_names
_CAMS_LEVS = {
    'O3' : np.array(([0., 20., 40., 60., 80., 100., 120.,
                      140., 160., 180., 200., 240., 500.])),
    'NO2': np.array(([0., 2., 5., 10., 20., 30., 40.,
                      50., 75., 100., 150., 200., 300.])),
    'SO2': np.array(([0., 2., 5., 10., 20., 30., 40.,
                      50., 75., 100., 150., 200., 800.])),
    'CO': np.array(([0., 50., 100., 150., 200., 250., 300.,
                     350., 400., 500., 700., 1000., 2000.])),
    'PM2p5': np.array(([0., 2., 5., 10., 20., 30., 40.,
                        50., 75., 100., 150., 200., 500.])),
    'PM10': np.array(([0., 2., 5., 10., 20., 30., 40.,
                       50., 75., 100., 150., 200., 500.])),
    }

#: Dictionary of contour levels for  CAMS
#: (Copernicus Atmospheric Monitoring Service), with a key of short_name
#: in ug/m3.
CAMS_LEVELS = {'O3' : _CAMS_LEVS['O3'],
               'NO2' : _CAMS_LEVS['NO2'],
               'NO' : _CAMS_LEVS['NO2'] * 0.1,
               'SO2' : _CAMS_LEVS['SO2'],
               'CO' : _CAMS_LEVS['CO'],
               'PM10': _CAMS_LEVS['PM10'],
               'PM10_NH42SO4': _CAMS_LEVS['PM10'],
               'PM10_BC': _CAMS_LEVS['PM10'],
               'PM10_BB': _CAMS_LEVS['PM10'],
               'PM10_OCFF': _CAMS_LEVS['PM10'],
               'PM10_DUST': _CAMS_LEVS['PM10'],
               'PM10_NH4NO3': _CAMS_LEVS['PM10'],
               'PM2p5': _CAMS_LEVS['PM2p5'],
               'PM2p5_NH42SO4': _CAMS_LEVS['PM2p5'],
               'PM2p5_BC': _CAMS_LEVS['PM2p5'],
               'PM2p5_BB': _CAMS_LEVS['PM2p5'],
               'PM2p5_OCFF': _CAMS_LEVS['PM2p5'],
               'PM2p5_DUST': _CAMS_LEVS['PM2p5'],
               'PM2p5_NH4NO3': _CAMS_LEVS['PM2p5']}

def cams_cmap():
    '''
    Generate colour map based on DAQI colours

    >>> cmap = cams_cmap()

    To plot and display colour map:

    >>> import config
    >>> import plotting_functions
    >>> filename = config.CODE_DIR + "/adaqdocs/figures/cams_cmap.png"
    >>> plotting_functions.plot_cmap(cmap, filename=filename)
    ... # doctest: +ELLIPSIS
    Saved figure .../cams_cmap.png

    .. image:: ../adaqdocs/figures/cams_cmap.png
       :scale: 50%

    '''

    #rgb colours [ [red, green, blue], [r,g,b],... ]
    colours = [[220, 220, 220],
               [18, 101, 209],
               [38, 129, 240],
               [79, 165, 244],
               [0, 220, 0],
               [80, 240, 80],
               [160, 229, 48],
               [229, 220, 49],
               [229, 175, 45],
               [240, 126, 38],
               [250, 60, 60],
               [240, 0, 132]]
    colours = np.array(colours)/255.
    cmap = colors.ListedColormap(colours, name='CAMS')

    #Register colour map, so can be called
    # eg plt.get_cmap('CAMS')
    cm.register_cmap(cmap=cmap)

    return cmap


#----- Ensure any colour maps are registered on loading aq_indices.py ----

daqi_cmap()
cams_cmap()

if __name__ == '__main__':

    import doctest
    doctest.testmod()
