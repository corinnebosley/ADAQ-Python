"""
Class for reading grib2 data from ECMWF into an ADAQData class.
Currently set up to work with C-IFS composition files and IFS meteorology
files, however it should work for, or be simple to adapt to, any other EC grib2
files.
"""
from __future__ import division
from __future__ import print_function

from six.moves.builtins import str

import os.path
import re
import datetime
import warnings

import numpy as np
import cf_units
import iris

import adaq_data
import array_statistics
import cube_time

#Conversion for CAMS ensemble constituent types to param and units
#(Note conversion from species given in filename for CAMS ensemble to
# short_name is also defined within get_fileinfo: if short_names are changed
# in the dictionary below they many also need to be modified in get_fileinfo)
CAMSENS_CONSTTYPE = {
    0:
    {'cube_name': 'mass_concentration_of_ozone_in_air',
     'units': cf_units.Unit('kg m-3'),
     'short_name': 'O3'},
    4:
    {'cube_name': 'mass_concentration_of_carbon_monoxide_in_air',
     'units': cf_units.Unit('kg m-3'),
     'short_name': 'CO'},
    5:
    {'cube_name': 'mass_concentration_of_nitrogen_dioxide_in_air',
     'units': cf_units.Unit('kg m-3'),
     'short_name': 'NO2'},
    8:
    {'cube_name': 'mass_concentration_of_sulfur_dioxide_in_air',
     'units': cf_units.Unit('kg m-3'),
     'short_name': 'SO2'},
    9:
    {'cube_name': 'mass_concentration_of_ammonia_in_air',
     'units': cf_units.Unit('kg m-3'),
     'short_name': 'NH3'},
    11:
    {'cube_name': 'mass_concentration_of_nitrogen_monoxide_in_air',
     'units': cf_units.Unit('kg m-3'),
     'short_name': 'NO'},
    40008:
    {'cube_name': 'mass_concentration_of_pm10_ambient_aerosol_in_air',
     'units': cf_units.Unit('kg m-3'),
     'short_name': 'PM10'},
    40009:
    {'cube_name': 'mass_concentration_of_pm2p5_ambient_aerosol_in_air',
     'units': cf_units.Unit('kg m-3'),
     'short_name': 'PM2p5'},
    60013:
    {'cube_name': 'mass_concentration_of_nmvoc_expressed_as_carbon_in_air',
     'units': cf_units.Unit('kg m-3'),
     'short_name': 'NMVOC'},
    60018:
    {'cube_name': 'mass_concentration_of_peroxyacetyl_nitrate_in_air',
     'units': cf_units.Unit('kg m-3'),
     'short_name': 'PAN'},
    64000:
    {'cube_name': 'grain_concentration_of_birch_pollen_in_air',
     'units': cf_units.Unit('m-3'),
     'short_name': 'birch_pollen'},
    64001:
    {'cube_name': 'grain_concentration_of_grass_pollen_in_air',
     'units': cf_units.Unit('m-3'),
     'short_name': 'grass_pollen'},
    64002:
    {'cube_name': 'grain_concentration_of_olive_pollen_in_air',
     'units': cf_units.Unit('m-3'),
     'short_name': 'olive_pollen'},
    64003:
    {'cube_name': 'grain_concentration_of_ragweed_pollen_in_air',
     'units': cf_units.Unit('m-3'),
     'short_name': 'ragweed_pollen'}
    }



GRIB_EC_DICT = {
    '31' :
    {'short_name' : 'sea_ice_area_frac'},
    '39' :
    {'cube_name': 'Volumetric_soil_water_layer_1',
     'units'    : 'm3 m-3',
     'short_name' : 'vsml1'},
    '40' :
    {'cube_name': 'Volumetric_soil_water_layer_2',
     'units'    : 'm3 m-3',
     'short_name' : 'vsml2'},
    '41' :
    {'cube_name': 'Volumetric_soil_water_layer_3',
     'units'    : 'm3 m-3',
     'short_name' : 'vsml3'},
    '42' :
    {'cube_name': 'Volumetric_soil_water_layer_4',
     'units'    : 'm3 m-3',
     'short_name' : 'vsml4'},
    '129':
    {'short_name' : 'z'},
    '130':
    {'short_name' : 'T'},
    '131':
    {'short_name' : 'u'},
    '132':
    {'short_name' : 'v'},
    '133':
    {'short_name' : 'q'},
    '134' :
    {'cube_name': 'surface_air_pressure',
     'units'    : 'Pa',
     'short_name' : 'p_msl'},
    '139' :
    {'cube_name': 'Soil_temperature_level_1',
     'units'    : 'K',
     'short_name' : 'soil_T_lev1'},
    '152' :
    {'cube_name': 'logarithm_of_surface_pressure',
     'units'    : '1',
     'short_name' : 'lnsp'},
    '170' :
    {'cube_name': 'Soil_temperature level_2',
     'units'    : 'K',
     'short_name' : 'soil_T_lev2'},
    '172' :
    {'short_name' : 'laf'},
    '175' :
    {'short_name' : 'lsm'},
    '183' :
    {'cube_name': 'Soil_temperature_level_3',
     'units'    : 'K',
     'short_name' : 'soil_T_lev3'},
    '235' :
    {'short_name': 'T_1p5'},
    '236' :
    {'cube_name': 'Soil_temperature level_4',
     'units'    : 'K',
     'short_name' : 'soil_T_lev4'},
    '246' :
    {'cube_name': 'Specific_cloud_liquid_water_content',
     'units'    : 'kg kg-1',
     'short_name' : 'clwc'},
    '247' :
    {'cube_name': 'Specific_cloud_ice_water_content',
     'units'    : 'kg kg-1',
     'short_name' : 'ciwc'},
    '248' :
    {'cube_name': 'Cloud_cover',
     'units'    : '1',
     'short_name' : 'cloud_cover'},
    '3066' :
    {'short_name' : 'snow_depth'},
    '210001' :
    {'cube_name': 'mass_fraction_of_sea_salt_(0.03_-_0.5_um)_aerosol_'
                  'particles_in_air',
     'units'    : 'kg kg-1',
     'short_name' : 'SS_ECbin1'},
    '210002' :
    {'cube_name': 'mass_fraction_of_sea_salt_(0.5_-_5_um)_aerosol_'
                  'particles_in_air',
     'units'    : 'kg kg-1',
     'short_name' : 'SS_ECbin2'},
    '210003' :
    {'cube_name': 'mass_fraction_of_sea_salt_(5_-_20_um)_aerosol_'
                  'particles_in_air',
     'units'    : 'kg kg-1',
     'short_name' : 'SS_ECbin3'},
    '210004' :
    {'cube_name': 'mass_fraction_of_dust_(0.03_-_0.55_um)_dry_aerosol_'
                  'particles_in_air',
     'units'    : 'kg kg-1',
     'short_name' : 'DUST_ECbin1'},
    '210005' :
    {'cube_name': 'mass_fraction_of_dust_(0.55_-_0.9_um)_dry_aerosol_'
                  'particles_in_air',
     'units'    : 'kg kg-1',
     'short_name' : 'DUST_ECbin2'},
    '210006' :
    {'cube_name': 'mass_fraction_of_dust_(0.9_-_20_um)_dry_aerosol_'
                  'particles_in_air',
     'units'    : 'kg kg-1',
     'short_name' : 'DUST_ECbin3'},
    '210007' :
    {'cube_name': 'mass_fraction_of_hydrophobic_organic_matter_'
                  'dry_aerosol_particles_in_air',
     'units' : 'kg kg-1',
     'short_name' : 'hydrophobic_OC'},
    '210008' :
    {'cube_name': 'mass_fraction_of_hydrophilic_organic_matter_'
                  'dry_aerosol_particles_in_air',
     'units' : 'kg kg-1',
     'short_name' : 'hydrophilic_OC'},
    '210009' :
    {'cube_name': 'mass_fraction_of_hydrophobic_black_carbon_dry_'
                  'aerosol_particles_in_air',
     'units' : 'kg kg-1',
     'short_name' : 'hydrophobic_BC'},
    '210010' :
    {'cube_name': 'mass_fraction_of_hydrophilic_black_carbon_dry_'
                  'aerosol_particles_in_air',
     'units' : 'kg kg-1',
     'short_name' : 'hydrophilic_BC'},
    '210011' :
    {'cube_name': 'mass_fraction_of_sulfate_dry_aerosol_particles_'
                  'in_air',
     'units' : 'kg kg-1',
     'short_name' : 'SO4'},
    '210073' :
    {'cube_name': 'mass_concentration_of_pm2p5_ambient_aerosol_in_air',
     'short_name' : 'PM2p5'},
    '210074' :
    {'cube_name': 'mass_concentration_of_pm10_ambient_aerosol_in_air',
     'short_name' : 'PM10'},
    '210122' :
    {'cube_name': 'mass_fraction_of_sulfur_dioxide_in_air',
     'units'    : 'kg kg-1',
     'short_name' : 'SO2'},
    '210121' :
    {'cube_name': 'mass_fraction_of_nitrogen_dioxide_in_air',
     'units'    : 'kg kg-1',
     'short_name' : 'NO2'},
    '210123' :
    {'cube_name': 'mass_fraction_of_carbon_monoxide_in_air',
     'units'    : 'kg kg-1',
     'short_name' : 'CO'},
    '210124' :
    {'cube_name': 'mass_fraction_of_formaldehyde_in_air',
     'units'    : 'kg kg-1',
     'short_name' : 'HCHO'},
    '210203' :
    {'cube_name': 'mass_fraction_of_ozone_in_air',
     'units'    : 'kg kg-1',
     'short_name' : 'O3'},
    '217004' :
    {'cube_name': 'mass_fraction_of_methane_in_air',
     'units'    : 'kg kg-1',
     'short_name' : 'CH4'},
    '217006' :
    {'cube_name': 'mass_fraction_of_nitric_acid_in_air',
     'units'    : 'kg kg-1',
     'short_name' : 'HNO3'},
    '217013' :
    {'cube_name': 'mass_fraction_of_peroxyacetyl_nitrate_in_air',
     'units'    : 'kg kg-1',
     'short_name' : 'PAN'},
    '217016' :
    {'cube_name': 'mass_fraction_of_isoprene_in_air',
     'units'    : 'kg kg-1',
     'short_name' : 'C5H8'},
    '217019' :
    {'cube_name': 'mass_fraction_of_ammonia_in_air',
     'units'    : 'kg kg-1',
     'short_name' : 'NH3'},
    '217027' :
    {'cube_name': 'mass_fraction_of_nitrogen_monoxide_in_air',
     'units'    : 'kg kg-1',
     'short_name' : 'NO'},
    '217045' :
    {'cube_name': 'mass_fraction_of_ethane_in_air',
     'units'    : 'kg kg-1',
     'short_name' : 'C2H6'},
    '217047' :
    {'cube_name': 'mass_fraction_of_propane_in_air',
     'units'    : 'kg kg-1',
     'short_name' : 'C3H8'}
    }


class ECGribData(adaq_data.ADAQData):
    """
    Subclass of ADAQData, which contains extra functionality specific to
    ECMWF grib2 format data.
    Also works with these cubes to interpolate to site locations and build
    sites_cube_list.

    .. note:: All files must be in grib2 format, so any grib1 format data
              (eg from IFS) must first be converted to grib2 format. This can
              be done using "grib_set -s edition=2 grib1file grib2file"

    .. note:: Currently set up to work with C-IFS composition files and IFS
              meteorology files, however the readdata() routine should work for,
              or be simple to adapt to, any other ECMWF grib2 files.

    **Example using grgcops (C-IFS) data:**

    Start by setting up dates and required directory to read from:

    >>> import config
    >>> datadir = config.SAMPLE_DATADIR + '/ECgrib/grgcops/'
    >>> start_dt = datetime.datetime(2016,8,15,0)
    >>> end_dt = datetime.datetime(2016,8,15,6)

    Now initialise class:

    >>> ecdata = ECGribData()

    Generate list of filenames:

    >>> filenames = ecdata.get_filenames(datadir, forecast_day=1,
    ... start_datetime=start_dt, end_datetime=end_dt)
    >>> print(filenames) # doctest: +ELLIPSIS
    [.../20160814T0000Z_grgcops_27.grib2', '...20160814T0000Z_grgcops_30.grib2']

    Now read in these files, just for the surface fields and for ozone and
    air temperature

    Note, by passing in some of the keywords to the get_filenames argument, they
    have now been set up as class attributes, so don't need to be passed in
    again to readdata.

    >>> gcl = ecdata.readdata(surface_only=True, short_name_list=['O3', 'T'])
    >>> print(gcl)
    0: air_temperature / (K)               (time: 2; latitude: 60; \
longitude: 137)
    1: mass_fraction_of_ozone_in_air / (kg kg-1) (time: 2; latitude: 60; \
longitude: 137)
    >>> o3 = gcl.extract(iris.AttributeConstraint(short_name='O3'),
    ... strict=True)
    >>> temperature = gcl.extract(iris.AttributeConstraint(short_name='T'),
    ... strict=True)
    >>> print(o3)  # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    mass_fraction_of_ozone_in_air / (kg kg-1) (time: 2; latitude: 60; \
longitude: 137)
         Dimension coordinates:
              time                                 x            -              -
              latitude                             -            x              -
              longitude                            -            -              x
         Auxiliary coordinates:
              forecast_period                      x            -              -
         Scalar coordinates:
              forecast_day: 1 Days
              level_pressure: 0.0 Pa
              model_level_number: 60
              sigma: 1.0
         Attributes:
              ...
              label: ECGrib
              paramID: 210203
              short_name: O3
    >>> print("{:.3e},{:.3e}".format(o3.data.min(), o3.data.max()))
    6.919e-12,2.993e-07
    >>> print("{:.3f},{:.3f}".format(temperature.data.min(),
    ... temperature.data.max()))
    247.147,315.541

    **Example using IFS meteorology data:**

    Here we just read in a single file to save time when running doctests:

    >>> datadir = config.SAMPLE_DATADIR + '/ECgrib/ifs/'
    >>> start_dt = datetime.datetime(2016,8,20,6)
    >>> end_dt = datetime.datetime(2016,8,20,9)

    For the purposes of the doctest, we assume that todays date is 30th August
    2016 to force the filename generation to correctly assume the file dates
    are for 2016.

    >>> today_dt = datetime.datetime(2016,8,30)

    Generate list of filenames:

    >>> ecdata = ECGribData()
    >>> filenames = ecdata.get_filenames(datadir, forecast_day=1,
    ... start_datetime=start_dt, end_datetime=end_dt, today_dt=today_dt)
    >>> print(filenames) # doctest: +ELLIPSIS
    ['/.../ifs/B1D08181200082006001.grib2']

    Read data at all model levels:

    >>> gcl = ecdata.readdata()
    >>> print(gcl)
    0: Soil_temperature level_4 / (K)      (latitude: 533; longitude: 1226)
    1: Soil_temperature_level_1 / (K)      (latitude: 533; longitude: 1226)
    2: Soil_temperature level_2 / (K)      (latitude: 533; longitude: 1226)
    3: Volumetric_soil_water_layer_3 / (m3 m-3) (latitude: 533; longitude: 1226)
    4: Volumetric_soil_water_layer_2 / (m3 m-3) (latitude: 533; longitude: 1226)
    5: Soil_temperature_level_3 / (K)      (latitude: 533; longitude: 1226)
    6: Volumetric_soil_water_layer_1 / (m3 m-3) (latitude: 533; longitude: 1226)
    7: Volumetric_soil_water_layer_4 / (m3 m-3) (latitude: 533; longitude: 1226)
    8: Specific_cloud_liquid_water_content / (kg kg-1) \
(model_level_number: 137; latitude: 533; longitude: 1226)
    9: Specific_cloud_ice_water_content / (kg kg-1) \
(model_level_number: 137; latitude: 533; longitude: 1226)
    10: Cloud_cover / (1)                   (model_level_number: 137; \
latitude: 533; longitude: 1226)
    11: air_temperature / (K)               (model_level_number: 137; \
latitude: 533; longitude: 1226)
    12: geopotential / (m2 s-2)             (latitude: 533; longitude: 1226)
    13: land_area_fraction / (1)            (latitude: 533; longitude: 1226)
    14: sea_ice_area_fraction / (1)         (latitude: 533; longitude: 1226)
    15: specific_humidity / (kg kg-1)       (model_level_number: 137; \
latitude: 533; longitude: 1226)
    16: surface_air_pressure / (Pa)         (latitude: 533; longitude: 1226)
    17: surface_temperature / (K)           (latitude: 533; longitude: 1226)
    18: thickness_of_snowfall_amount / (m)  (latitude: 533; longitude: 1226)
    19: x_wind / (m s-1)                    (model_level_number: 137; \
latitude: 533; longitude: 1226)
    20: y_wind / (m s-1)                    (model_level_number: 137; \
latitude: 533; longitude: 1226)

    Check the temperature cube:

    >>> temperature = gcl.extract(iris.AttributeConstraint(short_name='T'),
    ... strict=True)
    >>> print(temperature)  # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    air_temperature / (K)               (model_level_number: 137; latitude: 533; longitude: 1226)
         Dimension coordinates:
              model_level_number                           x              -               -
              latitude                                     -              x               -
              longitude                                    -              -               x
         Auxiliary coordinates:
              level_pressure                               x              -               -
              sigma                                        x              -               -
         Scalar coordinates:
              forecast_day: 1 Days
              forecast_period: 42 hours
              time: 2016-08-20 06:00:00
         Attributes:
              ...
              label: ECGrib
              paramID: 130
              short_name: T

    >>> print("{:.3f},{:.3f}".format(temperature.data.min(),
    ... temperature.data.max()))
    183.937,314.662

    This ECGribData class can also be used to read grib-format MACC/CAMS
    ensemble data:

    >>> datadir = config.SAMPLE_DATADIR + '/ECgrib/cams_ens/'
    >>> start_dt = datetime.datetime(2018,1,1,0)
    >>> end_dt = datetime.datetime(2018,1,3,0)
    >>> ecdata = ECGribData()
    >>> filenames = ecdata.get_filenames(datadir, forecast_day=1,
    ... start_datetime=start_dt, end_datetime=end_dt, short_name_list=['O3'])
    >>> gcl = ecdata.readdata()
    >>> print(gcl[0]) # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    mass_concentration_of_ozone_in_air / (kg m-3) (time: 48; latitude: 400; \
longitude: 700)
         Dimension coordinates:
              time                                     x    -    -
              latitude                                 -    x    -
              longitude                                -    -    x
         Auxiliary coordinates:
              forecast_period                          x    -    -
         Scalar coordinates:
              forecast_day: 1 Days
              sigma: 1.0
         Attributes:
              ...
              label: ECGrib
              paramID: 0
              short_name: O3

    >>> print("{:.3e},{:.3e}".format(gcl[0].data.min(), gcl[0].data.max()))
    7.282e-17,1.139e-07

    """

    def __init__(self, label='ECGrib'):
        """
        Initiates a class from EC grib data as a subclass of
        :class:`adaq_data.ADAQData`.
        """

        adaq_data.ADAQData.__init__(self)

        #Now some optional extra attributes specific to this type of data:
        self.label = label #Label
        self.level = None #Model level
        self.short_name_list = None #List of short names
        self.start_datetime = None #Starting datetime
        self.end_datetime = None #End datetime
        self.forecast_day = None #Forecast day
        self.filenames = None #List of raw data filenames
        self.sites_data = None #sites_data object
        self.model_level_number = None #model level

    def readdata(self, filenames=None, short_name_list=None, level=None,
                 start_datetime=None, end_datetime=None, label=None,
                 forecast_day=None, model_level_number=None,
                 surface_only=False):

        """
        Create a gridded_cube_list from filenames specified.

        :param filenames: list of files to read, can include wildcards etc
        :param short_name_list: list of short names (None gets everything)
        :param level: model level height - eg 0 for the surface at 0m.
        :param start_datetime: datetime object start date
        :param end_datetime: datetime object end date
        :param label: label for the cubes
        :param forecast_day: This can be an integer, whose number corresponds
                             to the forecast day, eg 1 if using the first full
                             forecast day from each model run.
                             Alternatively, set to "forecast" to get all
                             leadtimes from a single model run.
                             Or set to "latest" to get the latest available
                             data, which would be day 1 forecasts where
                             available and then the forecast from a single run
                             at the end.
        :param model_level_number: Model level number to read.
                                   Set to None to get all levels
        :param surface_only: Only read data on the surface level.
        """

        #Set class attributes on basis of keywords if set
        if filenames is not None:
            self.filenames = filenames
        if short_name_list is not None:
            self.short_name_list = short_name_list
        if forecast_day is not None:
            self.forecast_day = forecast_day
        if level is not None:
            self.level = level
        if start_datetime is not None:
            self.start_datetime = start_datetime
        if end_datetime is not None:
            self.end_datetime = end_datetime
        if label is not None:
            self.label = label
        if model_level_number is not None:
            self.model_level_number = model_level_number

        if self.filenames is None:
            raise ValueError('No filenames given')

        #Now set up iris constraints:
        constraints = None


        if self.short_name_list is not None:
            sname_constraint = iris.AttributeConstraint(
                short_name=lambda c: c in self.short_name_list)
            constraints = constraints & sname_constraint

        if self.start_datetime is not None:
            time_constraint = iris.Constraint(time=lambda c:
                                              c.point >= self.start_datetime)
            constraints = constraints & time_constraint

        if self.end_datetime is not None:
            time_constraint = iris.Constraint(time=lambda c:
                                              c.point <= self.end_datetime)
            constraints = constraints & time_constraint

        if not (self.forecast_day is None or \
                self.forecast_day == 'forecast' or \
                self.forecast_day == 'latest'):

            day_constraint = iris.Constraint(
                forecast_day=int(self.forecast_day))
            constraints = constraints & day_constraint

        if self.model_level_number is not None:
            level_constraint = iris.Constraint(
                model_level_number=self.model_level_number)
            constraints = constraints & level_constraint

        if surface_only:
            #Sigma should be 1.0 at the surface
            sigma_constraint = iris.Constraint(sigma=1.0)
            constraints = constraints & sigma_constraint

        #Finally read in data using iris load:

        self.gridded_cube_list = iris.load(self.filenames,
                                           constraints,
                                           callback=self.__callback)

        #Extract forecast days required through the setting of forecast_day.
        if self.forecast_day == 'latest' or self.forecast_day == 'forecast':
            for i, cube in enumerate(self.gridded_cube_list):
                self.gridded_cube_list[i] = \
                    cube_time.extract_latest_forecast_days(
                        cube, self.forecast_day, start_dt=self.start_datetime)
        elif self.forecast_day == 'all':
            self.gridded_cube_list = cube_time.extract_all_forecast_days(
                self.gridded_cube_list)

        #Add bounds to time coordinate
        for cube in self.gridded_cube_list:
            if cube.coord('time'):
                if len(cube.coord('time').points) > 1:
                    if not cube.coord('time').has_bounds():
                        #Assume points are given at the end of the
                        #meaning period
                        cube.coord('time').guess_bounds(bound_position=1.0)


        return self.gridded_cube_list


    def __callback(self, cube, field, filename):
        """
        Private method to provide callback to iris.load.
        Adds standard_name, short_name and units if known.
        Adds label attribute.
        Adds grib 'paramID' as an attribute
        Add forecast_day coordinate
        Fixes some issues with EC data model levels starting at top of model.
        """

        # add paramID to the cube attributes so that we can use it
        # for constrained loading etc
        section = field.sections[4]
        try:
            paramid = section['paramId']
        except KeyError:
            paramid = section.get_computed_key('paramId')
        cube.attributes['paramID'] = paramid

        # Indicate that for EC data the model level numbering starts
        # from 1 at the top of the model.
        for ml_coord in cube.coords('model_level_number'):
            ml_coord.attributes['positive'] = 'down'

        # Modify paramid if set to 0 for CAMS-ensemble grib data
        # in order to match numbers in GRIB_EC_DICT
        if paramid == 0 and os.path.basename(filename)[:17] == 'W_fr-meteofrance,':
            table_info = CAMSENS_CONSTTYPE[section['constituentType']]
            if not cube.coords('sigma'):
                #Add sigma coord = 1 to indicate a surface field
                sigma_coord = iris.coords.DimCoord(np.array([1.]),
                                                   units=cf_units.Unit('1'),
                                                   long_name='sigma')
                cube.add_aux_coord(sigma_coord)
            #Ensure lat and lon points are equal for all cubes to the same
            #degree of rounding, otherwise cubes will not concatenate with
            #each other
            cube.coord('longitude').points = np.around(
                cube.coord('longitude').points, 4)
            cube.coord('latitude').points = np.around(
                cube.coord('latitude').points, 4)
        elif str(paramid) in GRIB_EC_DICT:
            table_info = GRIB_EC_DICT[str(paramid)]
        if table_info:
            #Add cube name
            if 'cube_name' in table_info:
                cube.rename(table_info['cube_name'])

            # set units if unknown
            if cube.units == cf_units._UNKNOWN_UNIT_STRING:
                if 'units' in table_info:
                    cube.units = table_info['units']

            #Set short name
            if 'short_name' in table_info:
                cube.attributes['short_name'] = table_info['short_name']

        else:
            print('Unknown paramid:', paramid, cube.name())

        #Sort out sigma due to bug in iris - solution provided by AVD
        #(present at /opt/scitools/environments/default_legacy/2019_02_27)
        #Note this bug causes the lowest model level to have non-one value
        #for sigma - we rely on a value of sigma=1 to identify the surface
        #level.
        if cube.coords('sigma'):
            section = field.sections[4]
            template = section['productDefinitionTemplateNumber']
            if template in (0, 1, 8, 9, 10, 11, 15, 40):
                nv = section['NV']
                typeoffirstfixedsurface = section['typeOfFirstFixedSurface']
                if (nv > 0 and typeoffirstfixedsurface in (105, 119)):
                    # N.B. both 105 and 119 may be used for hybrid pressure
                    #levels
                    pv = section['pv']
                    scaled_value = section['scaledValueOfFirstFixedSurface']

                    offset = scaled_value
                    level_coord = cube.coord('level_pressure')
                    level_coord.points = pv[offset:offset+1]

                    offset = scaled_value + nv // 2
                    sigma_coord = cube.coord('sigma')
                    sigma_coord.points = pv[offset:offset+1]

        #Sort out lnsp (logarithm of surface pressure)
        #This is a surface field, however it has some of the vertical
        #coordinates incorrectly set.
        #For safety, remove the model_level_number
        #rather than resetting it as all other 2D surface cubes don't have
        #this coordinate at all.
        #Reset the value of sigma - this should be 1 for a surface field
        if cube.name() == 'logarithm_of_surface_pressure':
            #cube.remove_coord('level_pressure')
            #Model level number is set to 1, but really should be the maximum
            #model level
            cube.remove_coord('model_level_number')
            #Sigma is set to 0, but for a surface field should be 1.
            cube.coord('sigma').points = 1.0

        #Add forecast day:
        #Only checked with files that have only a single time in it.
        #However it *should* work with multiple-time files.
        #Note if ever any problems here, then see maccens_data.__callback for
        #possible solution.
        tunits = cube.coord('time').units

        ref_time = cube.coord('forecast_reference_time').points
        ref_time_hrs = np.array([t.hour for t in tunits.num2date(ref_time)])

        forecast_periods = cube.coord('forecast_period').points

        forecast_days = array_statistics.calc_forecast_day(forecast_periods,
                                                           ref_time_hrs)

        #However... Use assumption that:
        #For an AQUM forecast whose 'day1' forecast starts at 01Z on
        #2nd of month...
        #This must have used a grgcops forecast from 12Z or 00Z of the
        #previous day...
        #So set up forecast day so they day1 forecast from both represents
        #the same date...
        #This is the equivalent to subtracting 1 from the true forecast_day
        #However for CAMS-ensemble (files beginning 'W_fr-meteofrance,',
        #this is does not need to be done as they are for the correct day
        #already
        if os.path.basename(filename)[:17] != 'W_fr-meteofrance,':
            forecast_days -= 1

        cube.add_aux_coord(iris.coords.DimCoord(forecast_days,
                                                long_name='forecast_day',
                                                units='Days'))

        #Remove forecast_reference_time coordinate
        cube.remove_coord('forecast_reference_time')


        # Set label attribute of cube to value set in call to init
        cube.attributes['label'] = self.label


    def extract_sites(self, sites_data=None):
        """
        Extract site-specific data from gridded_cube_list into
        sites_cube_list, given

        :param sites_data: site information from :class:`sites_info.SitesInfo`

        """
        if sites_data is not None:
            self.sites_data = sites_data
        if self.sites_data is None:
            raise ValueError("ModelCube: Cannot extract sites without self.sites_data set")

        self.sites_cube_list = self.extract_scl_from_gridded(sites_data,
                                                             'longitude',
                                                             'latitude')

        if self.forecast_day == 'all':
            #Sort out so that different forecast_days are put into different
            #cubes (otherwise day 1,3,5 are in one cube, 0,2,4 in another)
            self.sites_cube_list = cube_time.extract_all_forecast_days(
                self.sites_cube_list)


        return self.sites_cube_list



    def get_fileinfo(self, filename, today_dt=None):
        """
        Generate dictionary of information that can be gathered from the
        filename.

        :param filename: string, name of file to gather information from.
        :param today_dt: datetime object containing todays date. This should
                         not normally be required as defaults to todays real
                         date otherwise, but can be useful for doctests.

        :returns: A dictionary of information that can be gathered
                  from the filename:

          * filename - filename
          * data_start - datetime format of the start time of data in file
          * data_end - datetime format of the end timeof data in file
          * forecast_periods - list of expected leadtimes in file (note
            leadtimes correspond to the end of the period for means)
          * forecast_days - list of integer days which data corresponds to
          * day1_date - date of day1 forecast
          * short_name - optional, if a species name is given within the
            filename string

        .. note:: For IFS data without any years defined in filename, the data
                  is assumed to be recent forecast data, so the year is set such
                  that the forecast run time is less than or equal to todays
                  date.

        Example of grgcops (C-IFS) data:

        >>> ecdata = ECGribData()
        >>> filename = '/path/to/data/20160815T0000Z_grgcops_27.grib2'
        >>> fileinfo = ecdata.get_fileinfo(filename)
        >>> fileinfo == {'forecast_days': [1],
        ... 'day1_date': datetime.date(2016, 8, 16),
        ... 'data_end': datetime.datetime(2016, 8, 16, 3, 0),
        ... 'data_start': datetime.datetime(2016, 8, 16, 3, 0),
        ... 'forecast_periods': [27],
        ... 'filename': '20160815T0000Z_grgcops_27.grib2'}
        True

        Example of IFS data retrieved containing year in filename:

        >>> filename = '/path/to/data/B1D201608181200081901001.grib2'
        >>> fileinfo = ecdata.get_fileinfo(filename)
        >>> fileinfo ==  {'forecast_days': [0],
        ... 'day1_date': datetime.date(2016, 8, 19),
        ... 'data_end': datetime.datetime(2016, 8, 19, 1, 0),
        ... 'data_start': datetime.datetime(2016, 8, 19, 1, 0),
        ... 'forecast_periods': [13],
        ... 'filename': 'B1D201608181200081901001.grib2'}
        True

        Example IFS data retrieved using routine dissemination feed, so does not
        contain year in filename. Here we assume the todays date is 20160820:

        >>> filename = '/path/to/data/B1D08181200082001001.grib2'
        >>> fileinfo = ecdata.get_fileinfo(filename,
        ... today_dt = datetime.datetime(2016,8,20))
        >>> fileinfo == {'forecast_days': [1],
        ... 'day1_date': datetime.date(2016, 8, 19),
        ... 'data_end': datetime.datetime(2016, 8, 20, 1, 0),
        ... 'data_start': datetime.datetime(2016, 8, 20, 1, 0),
        ... 'forecast_periods': [37],
        ... 'filename': 'B1D08181200082001001.grib2'}
        True

        Example grib-format MACC ensemble data:

        >>> filename = ('/path/to/data/W_fr-meteofrance,MODEL,ENSEMBLE+FORECAST'
        ... '+SURFACE+O3+0H24H_C_LFPW_20180101000000.grib2')
        >>> fileinfo = ecdata.get_fileinfo(filename)
        >>> fileinfo == {
        ... 'data_start': datetime.datetime(2018, 1, 1, 0, 0),
        ... 'data_end': datetime.datetime(2018, 1, 2, 0, 0),
        ... 'forecast_days': [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        ... 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        ... 'forecast_periods': [ 0,  1,  2,  3,  4,  5,  6,  7,  8,
        ... 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24],
        ... 'filename': ('W_fr-meteofrance,MODEL,ENSEMBLE+FORECAST+SURFACE+'
        ... 'O3+0H24H_C_LFPW_20180101000000.grib2'),
        ... 'short_name': 'O3'}
        True

        """

        fileinfo = {}
        #Remove path if attached
        filename = os.path.basename(filename)

        fileinfo['filename'] = filename
        split_file = filename.split('_')


        if filename[:17] == 'W_fr-meteofrance,':
            #CAMS-ensemble
            split_file_plus = filename.split('+')
            #print(split_file_plus)
            forecast_period_str = split_file_plus[4].split('_')[0]
            leadtime_start, leadtime_end, _ = forecast_period_str.split('H')
            ref_dt_str = split_file_plus[4].split('_')[3].split('.')[0]
            ref_dt = datetime.datetime.strptime(ref_dt_str, '%Y%m%d%H%M%S')

            fileinfo['data_start'] = ref_dt + datetime.timedelta(
                hours=int(leadtime_start))
            fileinfo['data_end'] = ref_dt + datetime.timedelta(
                hours=int(leadtime_end))
            fileinfo['forecast_periods'] = list(np.arange(int(leadtime_start),
                                                          int(leadtime_end)+1))
            fileinfo['forecast_days'] = list(
                array_statistics.calc_forecast_day(fileinfo['forecast_periods'],
                                                   ref_dt.hour,
                                                   day_start_hour=1))
            fileinfo['short_name'] = split_file_plus[3]
            #Convert short name used in file to ADAQpython short name if
            #different
            short_name_dict = {'PM25': 'PM2p5',
                               'PANS': 'PAN',
                               'OLIVEPOLLEN': 'olive_pollen',
                               'GRASSPOLLEN': 'grass_pollen',
                               'BIRCHPOLLEN': 'birch_pollen',
                               'RAGWEEDPOLLEN': 'ragweed_pollen'}
            if fileinfo['short_name'] in short_name_dict:
                fileinfo['short_name'] = short_name_dict[fileinfo['short_name']]


        elif len(split_file) > 1:
            if split_file[1] == 'grgcops':
                #Filename of the format yyyymmddThhmmZ_grgcops_LLL.grib2
                #where LLL is leadtime (upto 3 characters)
                #Each file has a single forecast time in it
                ref_dt_str = split_file[0]
                ref_time_dt = datetime.datetime.strptime(ref_dt_str,
                                                         '%Y%m%dT%H%MZ')

                leadtime = int(split_file[2].split('.')[0])
                leadtime_delta = datetime.timedelta(hours=leadtime)

                fileinfo['data_start'] = ref_time_dt + leadtime_delta
                fileinfo['data_end'] = fileinfo['data_start']
                fileinfo['forecast_periods'] = [leadtime]

                #Calculate forecast day.
                #Use assumption that:
                #For an AQUM forecast whose 'day1' forecast starts at 01Z on
                #2nd of month...
                #This must have used a grgcops forecast from 12Z or 00Z of the
                #previous day...
                #So set up forecast day so the day1 forecast from both
                #represents the same date...
                #Note both a 0Z and a 12Z run represents same day here, so don't
                #allow default day_start_hour to be set.
                #This is the equivalent to subtracting 1 from the true
                #forecast_day
                fileinfo['forecast_days'] = [
                    array_statistics.calc_forecast_day(leadtime,
                                                       ref_time_dt.hour,
                                                       day_start_hour=1)
                    -1]

                #Date of day1
                #Day 1 date is for day after reference time
                fileinfo['day1_date'] = ref_time_dt.date() + \
                                        datetime.timedelta(days=1)

        elif filename[:3] == 'B1D':
            #Filename of format B1DMMDDHHIImmddhhiiE
            #See http://www.ecmwf.int/en/forecasts/documentation-and-support/\
            # data-delivery/manage-transmission-ecpds/real-time-data-file
            # for naming convention information
            #MMDDHHII= month, date, hour, minute for data run time
            #mmddhhii= month, date, hour, minute for data valid time
            #Or of format (user-defined for ecflow suite):
            # B1DYYYYMMDDHHIImmddhhiiE
            #which has extra YYYY- data run year.
            if len(filename) == 26:
                #B1DMMDDHHIImmddhhiiE.grib2

                ref_dt_str = filename[3:11]
                ref_time_dt = datetime.datetime.strptime(ref_dt_str,
                                                         '%m%d%H%M')
                valid_dt_str = filename[11:19]
                valid_time_dt = datetime.datetime.strptime(valid_dt_str,
                                                           '%m%d%H%M')
                #Unknown year- assume it is recent forecast data, so set year
                #accordingly.
                if today_dt is None:
                    today_dt = datetime.datetime.now()
                #Assume the same year as today if ref time is before or
                #equal to today,
                #Otherwise set year to last year
                if ref_time_dt.replace(year=today_dt.year) <= today_dt:
                    year = today_dt.year
                else:
                    year = today_dt.year - 1
                ref_time_dt = ref_time_dt.replace(year=year)
                valid_time_dt = valid_time_dt.replace(year=year)
                #Check that valid_time_dt >= ref_time_dt, otherwise need to
                #increment year in valid_time_dt
                #eg ref time is 31st December 2015,
                #valid time therefore might be 1st Jan 2016.
                if valid_time_dt < ref_time_dt:
                    valid_time_dt = valid_time_dt.replace(year=year+1)


            elif len(filename) == 30:
                #B1DYYYYMMDDHHIImmddhhiiE.grib2

                ref_dt_str = filename[3:15]
                ref_time_dt = datetime.datetime.strptime(ref_dt_str,
                                                         '%Y%m%d%H%M')
                valid_dt_str = filename[15:23]
                valid_time_dt = datetime.datetime.strptime(valid_dt_str,
                                                           '%m%d%H%M')
                #Correct year for valid_time_dt
                valid_time_dt = valid_time_dt.replace(year=ref_time_dt.year)

            else:
                warnings.warn('Unknown filename format:' + filename)
                return fileinfo #Unknown format



            #Only a single hour in file:
            fileinfo['data_start'] = valid_time_dt
            fileinfo['data_end'] = fileinfo['data_start']
            #Lead time in hours
            leadtime = int((valid_time_dt - ref_time_dt).total_seconds()
                           / (60.*60.))
            fileinfo['forecast_periods'] = [leadtime]

            #Calculate forecast day.
            #Use assumption that:
            #For an AQUM forecast whose 'day1' forecast starts at 01Z on
            #2nd of month...
            #This must have used an IFS forecast from 12Z or 00Z of the
            #previous day...
            #So set up forecast day so the day1 forecast from both represents
            #the same date...
            #Note both a 0Z and a 12Z run represents same day here, so don't
            #allow default day_start_hour to be set.
            #This is the equivalent to subtracting 1 from the true forecast_day
            fileinfo['forecast_days'] = [
                array_statistics.calc_forecast_day(leadtime,
                                                   ref_time_dt.hour,
                                                   day_start_hour=1)
                -1]

            #Date of day1
            #Day 1 date is for day after reference time
            fileinfo['day1_date'] = ref_time_dt.date() + \
                                    datetime.timedelta(days=1)


        return fileinfo

    def get_filenames(self, directory, start_datetime=None, end_datetime=None,
                      forecast_day=None, today_dt=None, short_name_list=None):

        """
        Generate the list of filenames from specifically formatted filenames, by
        using the filename to determine the dates/hours that will be contained
        within the file, hence the returned filename list can be limited to only
        the required files.

        :param directory: Directory containing files.
                          Could also be a directory with a wildcard at the end,
                          eg '/my/directory/aqum_2014033*.grib2' to limit files
        :param start_datetime: Datetime formatted start date
        :param end_datetime: Datetime formatted end date
        :param forecast_day: This can be an integer, whose number corresponds
                             to the forecast day, eg 1 if using the first full
                             forecast day from each model run.
                             Alternatively, set to "forecast" to get all
                             leadtimes from a single model run.
        :param today_dt: datetime object containing todays date. This should
                         not normally be required as defaults to todays real
                         date otherwise, but can be useful for doctests.
        :param short_name_list: list of short names which are required which
                                can be used to limit filenames if a species
                                name is given within the filename string.


        All files are expected to be in grib2 format, with the '.grib2'
        extension.

        Currently allowed formatting of filenames:

          * yyyymmddThhmmZ_grgcops_LLL.grib2 - this should contain grgcops
            (C-IFS) composition data
          * B1DMMDDHHIImmddhhiiE.grib2 - IFS meteorology files from routine
            dissemination
          * B1DYYYYMMDDHHIImmddhhiiE.grib2 - IFS meteorology files which have
            been modified to include year within filename.
          * W_fr-meteofrance,MODEL,ENSEMBLE+FORECAST+SURFACE+O3
            +0H24H_C_LFPW_20180101000000.grib2
            - Archived CAMS/MACC ensemble files in grib2 format, see
            www.regional.atmosphere.copernicus.eu/doc/Guide_Numerical_Data_CAMS_new.pdf

        Example using grgcops data for day1 forecast:

        >>> import config
        >>> datadir = config.SAMPLE_DATADIR + '/ECgrib/'

        >>> ecdata = ECGribData()
        >>> filenames = ecdata.get_filenames(datadir + 'grgcops/',
        ... forecast_day=1, start_datetime=datetime.datetime(2016,8,15),
        ... end_datetime=datetime.datetime(2016,8,18))
        >>> print(filenames) # doctest: +ELLIPSIS
        ['.../grgcops/20160814T0000Z_grgcops_27.grib2', \
'.../20160814T0000Z_grgcops_30.grib2', '.../20160814T0000Z_grgcops_33.grib2', \
'.../20160814T0000Z_grgcops_36.grib2', '.../20160814T0000Z_grgcops_39.grib2', \
'.../20160814T0000Z_grgcops_42.grib2', '.../20160814T0000Z_grgcops_45.grib2', \
'.../20160814T0000Z_grgcops_48.grib2', '.../20160815T0000Z_grgcops_27.grib2', \
'.../20160815T0000Z_grgcops_30.grib2', ..., \
'.../20160816T0000Z_grgcops_27.grib2', '.../20160816T0000Z_grgcops_30.grib2', \
..., '.../20160816T0000Z_grgcops_48.grib2']

        For full forecast range instead:

        >>> ecdata = ECGribData()
        >>> filenames = ecdata.get_filenames(datadir + 'grgcops/',
        ... forecast_day='forecast',
        ... start_datetime=datetime.datetime(2016,8,16))
        >>> print(filenames) # doctest: +ELLIPSIS
        ['.../grgcops/20160815T0000Z_grgcops_24.grib2', \
'.../20160815T0000Z_grgcops_27.grib2', \
'.../20160815T0000Z_grgcops_30.grib2', \
..., \
'.../20160815T0000Z_grgcops_117.grib2', \
'.../20160815T0000Z_grgcops_120.grib2']

        And for the latest available data:

        >>> ecdata = ECGribData()
        >>> filenames = ecdata.get_filenames(datadir + 'grgcops/',
        ... forecast_day='latest', start_datetime=datetime.datetime(2016,8,16))
        >>> print(filenames) # doctest: +ELLIPSIS
        ['.../grgcops/20160814T0000Z_grgcops_48.grib2', \
'.../20160815T0000Z_grgcops_27.grib2', '.../20160815T0000Z_grgcops_30.grib2', \
'.../20160815T0000Z_grgcops_33.grib2', ..., \
'.../20160815T0000Z_grgcops_48.grib2', '.../20160816T0000Z_grgcops_18.grib2', \
'.../20160816T0000Z_grgcops_21.grib2', ..., \
'.../20160816T0000Z_grgcops_117.grib2', '.../20160816T0000Z_grgcops_120.grib2']

        Example using IFS data from routine dissemination. As no year is given
        in the files, the year is set such that the forecast run time is less
        than or equal to todays date.
        Here we assume the todays date is 20160820:

        >>> ecdata = ECGribData()
        >>> filenames = ecdata.get_filenames(datadir + 'ifs/',
        ... forecast_day='forecast',
        ... start_datetime=datetime.datetime(2016,8,19),
        ... today_dt=datetime.datetime(2016,8,20))
        >>> print(filenames) # doctest: +ELLIPSIS
        ['/.../B1D08181200081900001.grib2', '/.../B1D08181200081906001.grib2', \
'/.../B1D08181200081912001.grib2', '/.../B1D08181200081918001.grib2', \
'/.../B1D08181200082000001.grib2', '/.../B1D08181200082006001.grib2', \
..., \
'/.../B1D08181200082212001.grib2', '/.../B1D08181200082218001.grib2', \
'/.../B1D08181200082300001.grib2']

        """

        if forecast_day is not None:
            self.forecast_day = forecast_day
        if start_datetime is not None:
            self.start_datetime = start_datetime
        if end_datetime is not None:
            self.end_datetime = end_datetime
        if short_name_list is not None:
            self.short_name_list = short_name_list

        self.filenames = []
        if not directory.endswith('/'):
            directory += '/'
        listdir = os.listdir(directory)

        #Apply natural sorting - so *_102 goes after *_99,
        #rather than before *_18
        #Algorithm from
        #https://blog.codinghorror.com/sorting-for-humans-natural-sort-order/
        convert = lambda text: int(text) if text.isdigit() else text
        alphanum_key = lambda key: [convert(c)
                                    for c in re.split('([0-9]+)', key)]
        listdir.sort(key=alphanum_key)

        #Set up variable to hold the date of the last day 1 in all the files
        #(note as the files are checked in reverse order, this
        # will be in the first file found)
        #This will be used for limiting files for forecast_day='latest'
        last_day1_date = None

        #Loop through files in reverse order,
        # to make finding latest files easier
        for filename in reversed(listdir):

            if not filename.endswith('.grib2'):
                #Not a grib2 format file - discard
                continue

            fileinfo = self.get_fileinfo(filename, today_dt=today_dt)

            #Discard file if all data is before required start time
            if self.start_datetime is not None:
                if fileinfo['data_end'] is not None:
                    if fileinfo['data_end'] < self.start_datetime:
                        continue

            #Discard file if all data is after required end time
            if self.end_datetime is not None:
                if fileinfo['data_start'] is not None:
                    if fileinfo['data_start'] > self.end_datetime:
                        continue

            if self.forecast_day is not None:
                if self.forecast_day == 'forecast':
                    if fileinfo['day1_date'] != self.start_datetime.date():
                        continue
                elif self.forecast_day == 'latest':
                    #Keep all files which contain day 1
                    #Also keep all files which have their day1_date matching
                    # the last day1_date of all files. This only works
                    # if the files are done in reverse order as the last
                    # file in sorted order should have the latest date.
                    if last_day1_date is None:
                        #First time through loop - pick out the date for
                        #the day1 from this file. It should be the last
                        #day1 date as the files are checked in reverse
                        #order.
                        last_day1_date = fileinfo['day1_date']
                    if not (1 in fileinfo['forecast_days'] or \
                       last_day1_date == fileinfo['day1_date']):
                        #If the forecast_day is not day 1, or if the day1 date
                        #of the file is not the last_day1_date, then ignore
                        #this file.
                        continue
                elif self.forecast_day == 'all':
                    #Keep all files
                    pass
                else:
                    #Discard file if file does not contain required forecast_day
                    if fileinfo['forecast_days'] is not None:
                        if int(self.forecast_day) not in \
                           fileinfo['forecast_days']:
                            continue

            #Discard file if it doesn't contain a required short_name
            #(if available)
            if self.short_name_list is not None and 'short_name' in fileinfo:
                if fileinfo['short_name'] not in self.short_name_list:
                    continue

            self.filenames.append(directory+filename)
            #Now reverse back to get in sorted order again
            self.filenames.sort(key=alphanum_key)

        return self.filenames


if __name__ == '__main__':

    import doctest
    doctest.testmod()
