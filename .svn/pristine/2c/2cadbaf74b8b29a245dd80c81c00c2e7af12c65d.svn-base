"""
Class (and example) for working with UM pp files,
to read them into ADAQData class.
"""
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

from six.moves.builtins import str
from six.moves.builtins import range
from six.moves.builtins import object

import os
import datetime
import glob

import numpy as np
import iris

import adaq_data
import cube_time
import array_statistics
import shell_commands
import meteorology
import cube_chemistry

# This is a dictionary only to be used for loading pp data
# Allow the callback to add the correct name and also the units
# This is primarily intended for stashcodes not currently
# present in Iris and should be regularly reviewed to remove
# redundant values. Note ozone left in place to illustrate
# contents of dictionary in case needed for new stash codes
_METAREL_DICT = {
    'm01s34i001' :
        {'cube_name': 'mass_fraction_of_ozone_in_air',
         'units'    : 'kg kg-1'
        },
    'm01s34i069':
        {'cube_name' : 'mass_fraction_of_methanol_in_air',
         'units'     : 'kg kg-1'
         },
    'm01s34i092' :
        {'cube_name' : 'mass_fraction_of_ethene_in_air',
         'units'     : 'kg kg-1'
        },
    'm01s34i094' :
        {'cube_name' : 'mass_fraction_of_butane_in_air',
         'units'     : 'kg kg-1'
        },
    'm01s34i100' :
        {'cube_name' : 'mass_fraction_of_xylene_in_air',
         'units'     : 'kg kg-1'
        },
    'm01s50i228' :
        {'cube_name' : 'photolysis_rate_of_ozone_to_1D_oxygen_atom',
         'units'     : 's-1'
        },
    'm01s50i229' :
        {'cube_name' : 'photolysis_rate_of_nitrogen_dioxide',
         'units'     : 's-1'
        },
    'm01s50i230' :
        {'cube_name': 'mass_concentration_of_nmvoc_expressed_as_carbon_in_air',
         'units'    : 'ug m-3'
        }
}
# a dictionary of short names for pp used to add a short name
# to the cubes on load and to constrain the loading by short name
# These are our choices specific
# to AQ verification. The main restriction
# on this is that anything to be used as a shortname in our plotting and
# verification needs to be given a short name here
_STASH_2_SHORT_NAME = {
    'm01s00i004': 'theta',
    'm01s00i010': 'spec_hum',
    'm01s00i025': 'BL_depth',
    'm01s00i033': 'orog',
    'm01s00i090': 'MURK',
    'm01s00i101': 'SO2',
    'm01s00i103': 'SO4_AIT',
    'm01s00i104': 'SO4_ACC',
    'm01s00i105': 'SO4_DIS',
    'm01s00i107': 'NH3',
    'm01s00i108': 'BC_FRSH',
    'm01s00i109': 'BC_AGD',
    'm01s00i110': 'BC_CLD',
    'm01s00i111': 'BIO_FRSH',
    'm01s00i112': 'BIO_AGD',
    'm01s00i113': 'BIO_CLD',
    'm01s00i114': 'OCFF_FRSH',
    'm01s00i115': 'OCFF_AGD',
    'm01s00i116': 'OCFF_CLD',
    'm01s00i407': 'p_rho',
    'm01s00i408': 'p',
    'm01s03i225': 'u_10m',
    'm01s03i226': 'v_10m',
    'm01s03i236': 'T_1p5',
    'm01s03i237': 'q_1p5',
    'm01s03i245': 'rh_1p5',
    'm01s03i274': 'accum_nitrate_drydepflx',
    'm01s03i300': 'NH3_drydepflx',
    'm01s03i305': 'stableBL',
    'm01s03i307': 'well_mixedBL',
    'm01s05i226': 'total_precip',
    'm01s16i004': 'T',
    'm01s16i222': 'p_msl',
    'm01s17i220': 'PM10',
    'm01s17i221': 'PM2p5',
    'm01s17i222': 'PM10_NH42SO4',
    'm01s17i223': 'PM2p5_NH42SO4',
    'm01s17i224': 'PM10_BC',
    'm01s17i225': 'PM2p5_BC',
    'm01s17i226': 'PM10_BB',
    'm01s17i227': 'PM2p5_BB',
    'm01s17i228': 'PM10_OCFF',
    'm01s17i229': 'PM2p5_OCFF',
    'm01s17i230': 'PM10_SOA',
    'm01s17i231': 'PM2p5_SOA',
    'm01s17i232': 'PM10_SS',
    'm01s17i233': 'PM2p5_SS',
    'm01s17i234': 'PM10_DUST',
    'm01s17i235': 'PM2p5_DUST',
    'm01s17i236': 'PM10_NH4NO3',
    'm01s17i237': 'PM2p5_NH4NO3',
    'm01s34i001': 'O3',
    'm01s34i002': 'NO',
    'm01s34i003': 'NO3',
    'm01s34i004': 'NO2',
    'm01s34i005': 'N2O5',
    'm01s34i007': 'HNO3',
    'm01s34i008': 'H2O2',
    'm01s34i009': 'CH4',
    'm01s34i010': 'CO',
    'm01s34i011': 'HCHO',
    'm01s34i014': 'C2H6',
    'm01s34i016': 'MeCHO',
    'm01s34i017': 'PAN',
    'm01s34i018': 'C3H8',
    'm01s34i022': 'Me2CO',
    'm01s34i027': 'C5H8',
    'm01s34i069': 'CH3OH',
    'm01s34i092': 'C2H4',
    'm01s34i093': 'C3H6',
    'm01s34i094': 'C4H10',
    'm01s34i097': 'TOLUENE',
    'm01s34i100': 'XYLENE',
    'm01s34i991': 'CH3O2',
    'm01s34i993': 'HO2',
    'm01s34i995': 'OH',
    'm01s34i997': 'O1D',
    'm01s50i021': 'O3_drydep',
    'm01s50i079': 'CO_drydep',
    'm01s50i173': 'NO2_drydep',
    'm01s50i174': 'HNO3_drydep',
    'm01s50i179': 'PAN_drydep',
    'm01s50i188': 'N2O5_wetdep',
    'm01s50i190': 'HNO3_wetdep',
    'm01s50i192': 'HCHO_wetdep',
    'm01s50i227': 'NOy',
    'm01s50i228': 'jO1D',
    'm01s50i229': 'jNO2',
    'm01s50i230': 'NMVOC'
    }
#Dictionary for reverse lookups
_SHORT_NAME_2_STASH = {val: key for key, val in _STASH_2_SHORT_NAME.items()}

#Dictionary of derived quantities, referenced by short name.
#'components' of the derived value are given as stash codes.
#These stashcodes should also be included in the _STASH_2_SHORT_NAME dictionary
#(Also used for retrieving required stash codes from mass).
#'function' should be a function that can be called to calculate this short
#  name - it should take in gridded cube list as the first argument with
#  any other keyword arguments passed in using the 'kwargs' dictionary. It
#  should return the gridded cube list with the new cube added.
_DERIVED_SHORT_NAME_DICT = {
    'ws_10m': {
        'components': ['u_10m', 'v_10m'],
        'standard_name': 'wind_speed',
        'function': meteorology.wind_speed,
        'kwargs': {'short_name': 'ws_10m'}
        },
    'wind_dir_10m' : {
        'components': ['u_10m', 'v_10m'],
        'standard_name': 'wind_to_direction',
        'function': meteorology.wind_direction,
        'kwargs': {'short_name': 'wind_dir_10m'}
        },
    'NOx': {
        'components': ['NO', 'NO2'],
        'standard_name': 'mole_fraction_of_nox_in_air',
        'function': cube_chemistry.calculate_nox
        },
    'Ox': {
        'components': ['O3', 'NO2'],
        'standard_name': 'mole_fraction_of_ox_in_air',
        'function': cube_chemistry.calculate_ox
        },
    'PM10_NO3': {
        'components': ['PM10_NH4NO3'],
        'function': cube_chemistry.calculate_no3_aerosol,
        'kwargs': {'pm': 'PM10'}
        },
    'PM2p5_NO3': {
        'components': ['PM2p5_NH4NO3'],
        'function': cube_chemistry.calculate_no3_aerosol,
        'kwargs': {'pm': 'PM2p5'}
        },
    'PM10_NH4': {
        'components': ['PM10_NH4NO3', 'PM10_NH42SO4'],
        'function': cube_chemistry.calculate_nh4_aerosol,
        'kwargs': {'pm': 'PM10'}
        },
    'PM2p5_NH4': {
        'components': ['PM2p5_NH4NO3', 'PM2p5_NH42SO4'],
        'function': cube_chemistry.calculate_nh4_aerosol,
        'kwargs': {'pm': 'PM2p5'}
        },
    'PM10_SO4': {
        'components': ['PM10_NH42SO4'],
        'function': cube_chemistry.calculate_so4_aerosol,
        'kwargs': {'pm': 'PM10'}
        },
    'PM2p5_SO4': {
        'components': ['PM2p5_NH42SO4'],
        'function': cube_chemistry.calculate_so4_aerosol,
        'kwargs': {'pm': 'PM2p5'}
        },
    'Total_VOC': {
        'components': ['HCHO', 'C2H6', 'MeCHO', 'C3H8',
                       'Me2CO', 'C5H8', 'C3H6', 'TOLUENE', 'CH3OH',
                       'C2H4', 'C4H10', 'XYLENE'],
        'function': cube_chemistry.calculate_total_voc
        }
    }


class PPData(adaq_data.ADAQData):

    """
    Subclass of ADAQData, which contains extra functionality specific to
    pp model data. In particular reads in pp files into an iris cubelist
    (one cube per species). Also works with these cubes to interpolate to
    site locations and build sites_cube_list.

    Example:

    Start by setting up dates and required pp files:

    >>> import config
    >>> import datetime
    >>> start_dt = datetime.datetime(2014, 4, 4, 00) #0Z on 4th April
    >>> end_dt = datetime.datetime(2014, 4, 5, 6) #6Z on 5th April
    >>> sample_data_path = config.SAMPLE_DATADIR+'aqum_output/oper/'
    >>> ppfiles = [sample_data_path+'aqum_20140403_T+000_QM18.pp',
    ... sample_data_path+'aqum_20140403_T+012_QM18.pp',
    ... sample_data_path+'aqum_20140403_T+024_QM18.pp',
    ... sample_data_path+'aqum_20140404_T+000_QM18.pp']

    Now initialise class:

    >>> pp = PPData()

    Read pp data into iris cubes:

    >>> gcl = pp.readdata(ppfiles, start_datetime=start_dt,
    ... end_datetime=end_dt, short_name_list=['O3','NO2'],
    ... forecast_day=1)

    This produces ADAQData.gridded_cube_list cubes:

    >>> print(type(pp.gridded_cube_list))
    <class 'iris.cube.CubeList'>
    >>> print(pp.gridded_cube_list)
    0: mass_fraction_of_nitrogen_dioxide_in_air / (kg kg-1) \
(time: 30; grid_latitude: 182; grid_longitude: 146)
    1: mass_fraction_of_ozone_in_air / (kg kg-1) \
(time: 30; grid_latitude: 182; grid_longitude: 146)


    To just look at eg O3 cube:

    >>> o3_gridded_cube = pp.extract(short_name='O3',
    ... gridded=True,singlecube=True)
    >>> print("{:.3e}".format(o3_gridded_cube.data.max()))
    1.080e-07

    To extract site data, first read in site information:

    >>> import sites_info
    >>> sitesfilename = (config.SAMPLE_DATADIR +
    ...     'AURN_obs/aq_sites_GEMSRAQ_v4b_dev.txt')
    >>> obsdir = config.SAMPLE_DATADIR+'AURN_obs/'
    >>> sites = sites_info.SitesInfo()
    >>> sites_data = sites.read_from_file(sitesfilename,
    ... allsites=False,obsdir=obsdir)
    Number of sites:  5

    Extract all sites from this list:

    >>> sites_cube_list = pp.extract_sites(sites_data=sites_data)
    >>> print(sites_cube_list)
    0: mass_fraction_of_nitrogen_dioxide_in_air / (kg kg-1) \
(site_id: 5; time: 30)
    1: mass_fraction_of_ozone_in_air / (kg kg-1) \
(site_id: 5; time: 30)

    Extract a single cube for just Harwell ozone

    >>> har_cube = pp.extract(short_name='O3', abbrev='HAR', singlecube=True)
    >>> print("{:.3e}".format(har_cube.data.max()))
    6.767e-08

    3D data can also be read into gridded cubes:

    >>> start_dt = datetime.datetime(2014, 8, 31, 00)
    >>> end_dt = datetime.datetime(2014, 9, 4, 00)
    >>> sample_data_path = config.SAMPLE_DATADIR+'aqum_output/oper/'
    >>> ppfiles = [sample_data_path+'prodm_op_aqum_20140901_18.000.pp']
    >>> pp = PPData ()
    >>> pp.readdata( ppfiles, start_datetime=start_dt,
    ...     end_datetime=end_dt, short_name_list=['O3','p','T'] )
    [<iris 'Cube' of air_pressure / (Pa) \
(time: 2; model_level_number: 38; grid_latitude: 182; grid_longitude: 146)>,
    <iris 'Cube' of air_temperature / (K) \
(time: 2; model_level_number: 38; grid_latitude: 182; grid_longitude: 146)>,
    <iris 'Cube' of mass_fraction_of_ozone_in_air / (kg kg-1) \
(time: 2; model_level_number: 38; grid_latitude: 182; grid_longitude: 146)>]

    Fields which are derived from multiple stash codes, eg wind speed from
    u and v components of wind can also be read and immediately calculated:

    >>> start_dt = datetime.datetime(2017,6,3,0)
    >>> end_dt = datetime.datetime(2017,6,4,0)
    >>> ppfiles = [sample_data_path+'met/prods*.pp']
    >>> pp.readdata(ppfiles, start_datetime=start_dt,
    ... end_datetime=end_dt, short_name_list=['ws_10m'] ) # doctest: +ELLIPSIS
    [<iris 'Cube' of wind_speed / (m s-1) (time: 25; grid_latitude: 183; \
grid_longitude: 146)>]

    If trying to read in AQUM data for specific forecast ranges, eg Day 1,
    make use of the AQUMppFiles class to generate filenames and
    forecast_periods:

    >>> ppfiles = AQUMppFiles()
    >>> start_dt = datetime.datetime(2014, 4, 4, 0)
    >>> end_dt = datetime.datetime(2014, 4, 5, 6)
    >>> filenames = ppfiles.get_filenames(sample_data_path, start_dt, end_dt,
    ... forecast_day=1)
    >>> pp = PPData()
    >>> pp.readdata(filenames=filenames, label='AQUM',
    ... start_datetime=start_dt, end_datetime=end_dt,
    ... short_name_list=['O3'], forecast_day=1)
    [<iris 'Cube' of mass_fraction_of_ozone_in_air / (kg kg-1) \
(time: 31; grid_latitude: 182; grid_longitude: 146)>]

    If instead trying to read in AQUM data for a single forecast run, use
    forecast_day='forecast' and can also use AQUMppFiles class to generate
    filenames. For more information about this
    setting, see :func:`cube_time.extract_latest_forecast_days`.

    To read in files from forecast day 1 where available, otherwise from
    the final forecast can use get_filenames, or just use readdata directly,
    setting the keyword forecast_day='latest'. For more information about this
    setting, see :func:`cube_time.extract_latest_forecast_days`.

    To read in files with all available forecast days and then force each
    different forecast day into a different cube, use the keyword
    forecast_day='all'. For more information about this setting see
    :func:`cube_time.extract_all_forecast_days`.

    """

    DERIVED_SHORT_NAME_DICT = _DERIVED_SHORT_NAME_DICT

    def __init__(self, label='pp'):
        """
        Initiates class as a subset of adaq.ADAQData, plus other
        model-specific data
        """

        adaq_data.ADAQData.__init__(self)

        self.label = label #Label
        self.short_name_list = None #List of short names
        self.start_datetime = None #Starting datetime
        self.end_datetime = None #End datetime
        self.filenames = None #List of raw data filenames
        self.forecast_day = None #Forecast day
        self.sites_data = None #sites_data object
        self.model_level_number = None #model level number
        self.surface_only = False #only load surface data

    def readdata(self, filenames=None, short_name_list=None,
                 start_datetime=None, end_datetime=None, label=None,
                 forecast_day=None, model_level_number=None,
                 surface_only=False):

        """
        Create an Iris cube from the filename(s) specified.

        :param filenames: List of filenames to read from, can include wildcards.
        :param short_name_list: List of short names to read
        :param start_datetime: datetime formatted start date and time for cubes.
                               If not given, then minimum date taken from all
                               data read in.
        :param end_datetime: datetime formatted end date and time for cubes
                               If not given, then maximum date taken from all
                               data read in.
        :param label: String to be added to cube attributes as 'label'.
        :param forecast_day: This can be an integer, whose number corresponds
                             to the forecast day, eg 1 if using the first full
                             forecast day from each model run.
                             Alternatively, set to "forecast" to get all
                             leadtimes from a single model run.
                             Or set to "latest" to get the latest available
                             data, which would be day 1 forecasts where
                             available and then the forecast from a single run
                             at the end.
        :param model_level_number: Integer number of single model level to
                                   read in. If None, then all levels read.
        :param surface_only: Logical. If True, then only model_level_number=1
                             and all single level fields are read in.

        :returns: list of gridded cubes, one per short name.
        """

        #Change input values into self. values
        if filenames is not None:
            self.filenames = filenames
        if short_name_list is not None:
            self.short_name_list = short_name_list
        if forecast_day is not None:
            self.forecast_day = forecast_day
        if label is not None:
            self.label = label
        if start_datetime is not None:
            self.start_datetime = start_datetime
        if end_datetime is not None:
            self.end_datetime = end_datetime
        if model_level_number is not None:
            self.model_level_number = model_level_number
        if surface_only:
            self.surface_only = True

        #Raise errors if key information is missing
        if self.filenames is None:
            raise ValueError("ModelCube: no filename(s) specified")
        if self.short_name_list is None:
            raise ValueError("ModelCube: no short_name_list has been specified")

        #Set up time constraint using start and end datetimes
        if self.start_datetime is not None and self.end_datetime is not None:
            time_constraint = iris.Constraint(
                time=lambda c: c.point >= self.start_datetime
                and c.point <= self.end_datetime)
        elif self.start_datetime is not None:
            time_constraint = iris.Constraint(
                time=lambda c: c.point >= self.start_datetime)
        elif self.end_datetime is not None:
            time_constraint = iris.Constraint(
                time=lambda c: c.point <= self.end_datetime)
        else:
            time_constraint = iris.Constraint(time=lambda c: True)

        constraints = None

        #Set up forecast_day constraint
        if not (self.forecast_day is None or \
                self.forecast_day == 'forecast' or \
                self.forecast_day == 'latest' or \
                self.forecast_day == 'all'):
            fcst_day_constraint = iris.Constraint(forecast_day
                                                  =int(self.forecast_day))
            constraints = constraints & fcst_day_constraint

        #Set up short name constraint
        specifications = self.get_derived_sn_specs(short_name_list)
        dependencies = set(short_name_list)
        for spec in specifications:
            if 'components' in spec:
                dependencies.update(spec['components'])
        constraints &= iris.AttributeConstraint(
            short_name=lambda c: c in dependencies)

        #Set up model level constraint
        if self.model_level_number is not None:
            model_lev_constraint = iris.Constraint(model_level_number=
                                                   self.model_level_number)
            constraints = constraints & model_lev_constraint


        #Load cubes in - constrain using forecast period (leadtime) and stash.
        self.gridded_cube_list = iris.load(self.filenames, constraints,
                                           callback=self.__callback)

        if not isinstance(self.gridded_cube_list, iris.cube.CubeList):
            #If single cube read in, convert to a cubelist anyway
            self.gridded_cube_list = iris.cube.CubeList(self.gridded_cube_list)

        self.gridded_cube_list = self.gridded_cube_list.extract(time_constraint)

        #Extract forecast days required through the setting of forecast_day.
        if self.forecast_day == 'latest' or self.forecast_day == 'forecast':
            for i, cube in enumerate(self.gridded_cube_list):
                self.gridded_cube_list[i] = \
                    cube_time.extract_latest_forecast_days(
                        cube, self.forecast_day, start_dt=self.start_datetime)
        elif self.forecast_day == 'all':
            self.gridded_cube_list = cube_time.extract_all_forecast_days(
                self.gridded_cube_list)

        if not self.gridded_cube_list:
            raise ValueError('No cubes found')

        #Create any derived items requested
        self.derive(*specifications)

        #Remove any unwanted cubes used in derivation
        self.gridded_cube_list = self.gridded_cube_list.extract(
            iris.AttributeConstraint(
                short_name=lambda c: c in short_name_list))

        return self.gridded_cube_list

    def __callback(self, cube, field, filename):

        """
        Private method to remove unwanted coords and to check bounds
        Also:
            Adds units if unknown
            Adds standard name if needed
            Adds 'short_name' for easier identification
            Sets time points to be at the end of the meaning period
              (if the data is a mean)
        """

        # set some local variable from cube attributes
        stashcode = str(cube.attributes['STASH'])
        coord_list = cube.coords()
        cubename = cube.name()

        # Use our version of metarelate to set long names
        # and units if not set by Iris load
        if stashcode in _METAREL_DICT:
            # only set name if Iris failed to set a name
            # over time this should happen less often as the vocabulary
            # of metarelate expands
            stash_info = _METAREL_DICT[stashcode]

            if cubename == stashcode:
                # use the rename function
                cube.rename(stash_info['cube_name'])

            # set units if unknown
            if cube.units.is_unknown():
                cube.units = stash_info['units']

            if cubename == 'mass_fraction_of_lumped_chlorine_expressed_as_hydrogen_chloride':
                cube.rename('mass_fraction_of_xylene_in_air')
                cube.units = 'kg kg-1'


        if stashcode in _STASH_2_SHORT_NAME:
            cube.attributes['short_name'] = _STASH_2_SHORT_NAME[stashcode]
        else:
            cube.attributes['short_name'] = cubename

        #Ignore non-surface data if required
        if self.surface_only:
            if cube.coords('model_level_number'):
                if cube.coord('model_level_number').points[0] != 1:
                    print('Ignoring cube:', cube)
                    raise iris.exceptions.IgnoreCubeException
                #Note any single level data without a model_level_number coordinate
                #will be allowed to remain, so if at a higher single level, then
                #this will also remain in loaded cubes.

        # removal of forecast reference time needed at present
        # to prevent failures (otherwise produces 4d cube with
        #  forecast_period and forecast_reference_time on separate dims)
        # but first get runtime hour from the reference time for later use
        for coord in coord_list:
            if coord.standard_name == 'forecast_reference_time':
                runtime_hr = coord.units.num2date(coord.points[0]).hour
                cube.remove_coord('forecast_reference_time')

        # time bounds - ignore cube on load if
        # meaning period is greater than 1 hour
        #NB. should be able to do this with cell_methods, but
        #currently 24 hour means are also being set with
        #cell_method = "mean (1 hour)"
        tcoord = cube.coord('time')
        if tcoord.has_bounds() and \
           (tcoord.bounds[0][1] - tcoord.bounds[0][0]) > 1.00001:
            raise iris.exceptions.IgnoreCubeException

        #Ignore SO2 15mins
        if stashcode == 'm01s00i101':
            if not cube.cell_methods:
                raise iris.exceptions.IgnoreCubeException

        #If meaned data, then points on eventual time-dimension
        # should be at the end of the meaning period for pp data.
        # This corresponds to the end of the bound.
        # Currently only time and forecast_period are known to
        # be on time-dimension.
        # Note time-dimension not created at this stage as for single field.
        if field.lbproc == 128: #meaned data
            for coord_name in ['time', 'forecast_period']:
                coord = cube.coord(coord_name)
                if coord.has_bounds():
                    coord.points = coord.bounds[0][1]

        #Set up forecast_day coordinate
        if not cube.coords('forecast_day'):
            forecast_period = cube.coord('forecast_period').points[0]
            forecast_day = array_statistics.calc_forecast_day(forecast_period,
                                                              runtime_hr)
            cube.add_aux_coord(iris.coords.DimCoord(forecast_day,
                                                    long_name='forecast_day',
                                                    units='Days'))

        # Set label attribute of cube to value set in call to init
        cube.attributes['label'] = self.label

        if 'um_version' in cube.attributes:
            del cube.attributes['um_version']



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
                                                             'grid_longitude',
                                                             'grid_latitude')

        if self.forecast_day == 'all':
            #Sort out so that different forecast_days are put into different
            #cubes (otherwise day 1,3,5 are in one cube, 0,2,4 in another)
            self.sites_cube_list = cube_time.extract_all_forecast_days(
                self.sites_cube_list)


        return self.sites_cube_list

class AQUMppFiles(object):
    """
    Class for containing information about AQUM-specific pp-format filenames

    >>> import config
    >>> sampledir = config.SAMPLE_DATADIR+'aqum_output/oper/'
    >>> AQFiles = AQUMppFiles()
    >>> filenames = AQFiles.get_filenames(directory = sampledir,
    ...     start_datetime=datetime.datetime(2014,4,1),
    ...     end_datetime=datetime.datetime(2014,4,3),
    ...     forecast_day=1)
    >>> print(filenames) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
    ['/.../aqum_output/oper/aqum_20140330_T+024_QM18.pp',\
    '/.../aqum_output/oper/aqum_20140331_T+000_QM18.pp',\
    '/.../aqum_output/oper/aqum_20140331_T+012_QM18.pp',\
    '/.../aqum_output/oper/aqum_20140331_T+024_QM18.pp',\
    '/.../aqum_output/oper/aqum_20140401_T+000_QM18.pp',\
    '/.../aqum_output/oper/aqum_20140401_T+012_QM18.pp',\
    '/.../aqum_output/oper/aqum_20140401_T+024_QM18.pp',\
    '/.../aqum_output/oper/aqum_20140402_T+000_QM18.pp']
    """

    def __init__(self):

        self.directory = None
        self.start_datetime = None
        self.end_datetime = None
        self.runtimes = None
        self.forecast_day = None
        self.forecast_periods = None #leadtime, in hours (to match data in cube)
        self.fileinfo_list = []
        self.filenames = None
        self.file_pattern = None

    def get_filenames(self, directory=None, start_datetime=None,
                      end_datetime=None, forecast_day=None):
        """
        Gets the filenames from AQUM-style files.
        Uses the filename to determine the dates/hours that will be contained
        within the file, hence the returned filename list can be limited
        to only the required files

        :param directory: Directory containing files.
                          Could also be a directory with a wildcard at the end,
                          eg '/my/directory/aqum_2014033*.pp' to limit files
        :param start_datetime: Datetime formatted start date
        :param end_datetime: Datetime formatted end date
        :param forecast_day: This can be an integer, whose number corresponds
                             to the forecast day, eg 1 if using the first full
                             forecast day from each model run.
                             Alternatively, set to "forecast" to get all
                             leadtimes from a single model run.

        >>> import config
        >>> sampledir = config.SAMPLE_DATADIR+'aqum_output/oper/'
        >>> AQFiles = AQUMppFiles()
        >>> filenames = AQFiles.get_filenames(directory = sampledir,
        ...     start_datetime=datetime.datetime(2014,4,1),
        ...     end_datetime=datetime.datetime(2014,4,2))
        >>> print(filenames) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        ['/.../aqum_output/oper/aqum_20140330_T+024_QM18.pp', \
        '/.../aqum_output/oper/aqum_20140331_T+000_QM18.pp', \
        '/.../aqum_output/oper/aqum_20140331_T+012_QM18.pp', \
        '/.../aqum_output/oper/aqum_20140331_T+024_QM18.pp', \
        '/.../aqum_output/oper/aqum_20140401_T+000_QM18.pp']

        Could also pass in wildcard options to directory, eg:

        >>> AQFiles = AQUMppFiles()
        >>> filenames = AQFiles.get_filenames(
        ... directory = sampledir+'*20140401*.pp',
        ... start_datetime=datetime.datetime(2014,4,1),
        ... end_datetime=datetime.datetime(2014,4,3),
        ... forecast_day=1)
        >>> print(filenames) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        ['.../oper/aqum_20140401_T+000_QM18.pp',
        '.../oper/aqum_20140401_T+012_QM18.pp',
        '.../oper/aqum_20140401_T+024_QM18.pp']

        Filenames can be limited by the required forecast day, eg for
        the 1st full day of a forecast:

        >>> AQFiles = AQUMppFiles()
        >>> filenames = AQFiles.get_filenames(
        ... directory = sampledir,
        ... start_datetime=datetime.datetime(2014,4,1),
        ... end_datetime=datetime.datetime(2014,4,3),
        ... forecast_day=1)
        >>> print(filenames) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        ['.../oper/aqum_20140330_T+024_QM18.pp',
        '.../oper/aqum_20140331_T+000_QM18.pp',
        '.../oper/aqum_20140331_T+012_QM18.pp',
        '.../oper/aqum_20140331_T+024_QM18.pp',
        '.../oper/aqum_20140401_T+000_QM18.pp',
        '.../oper/aqum_20140401_T+012_QM18.pp',
        '.../oper/aqum_20140401_T+024_QM18.pp',
        '.../oper/aqum_20140402_T+000_QM18.pp']

        Or to get filenames corresponding to a single forecast run:

        >>> AQFiles = AQUMppFiles()
        >>> sampledir = config.SAMPLE_DATADIR+'aqum_output/oper_forecast/'
        >>> filenames = AQFiles.get_filenames(
        ... directory = sampledir,
        ... start_datetime=datetime.datetime(2014,3,28),
        ... forecast_day='forecast')
        >>> print(filenames) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        ['.../oper_forecast/prods_op_aqum_20140327_18.000.pp',
        '.../oper_forecast/prods_op_aqum_20140327_18.012.pp',
        '.../oper_forecast/prods_op_aqum_20140327_18.024.pp',
        '.../oper_forecast/prods_op_aqum_20140327_18.036.pp',
        ...,
        '.../oper_forecast/prods_op_aqum_20140327_18.096.pp',
        '.../oper_forecast/prods_op_aqum_20140327_18.108.pp',
        '.../oper_forecast/prods_op_aqum_20140327_18.120.pp']

        Or for filenames corresponding to the 'latest' available data:
        >>> AQFiles = AQUMppFiles()
        >>> sampledir = config.SAMPLE_DATADIR+'aqum_output/oper_forecast/'
        >>> filenames = AQFiles.get_filenames(
        ... directory = sampledir,
        ... start_datetime=datetime.datetime(2014,4,1),
        ... forecast_day='latest')
        >>> print(filenames) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        ['.../prods_op_aqum_20140330_18.024.pp',
        '.../prods_op_aqum_20140331_18.000.pp',
        '.../prods_op_aqum_20140331_18.012.pp',
        ...
        '.../prods_op_aqum_20140402_18.012.pp',
        '.../prods_op_aqum_20140402_18.024.pp',
        '.../prods_op_aqum_20140403_18.000.pp',
        '.../prods_op_aqum_20140403_18.012.pp',
        '.../prods_op_aqum_20140403_18.024.pp',
        '.../prods_op_aqum_20140403_18.036.pp',
        ...
        '.../prods_op_aqum_20140403_18.108.pp',
        '.../prods_op_aqum_20140403_18.120.pp']
        """

        #Set values into class
        if directory is not None:
            self.directory = directory
        if start_datetime is not None:
            self.start_datetime = start_datetime
        if end_datetime is not None:
            self.end_datetime = end_datetime
        if forecast_day is not None:
            self.forecast_day = forecast_day

        if not os.path.isdir(directory):
            #Path is not a directory - assume that this is a file pattern
            #to match instead - split up directory
            self.directory = os.path.dirname(directory)
            self.file_pattern = os.path.basename(directory)
        else:
            self.file_pattern = '*.pp'

        if not os.path.exists(self.directory):
            raise ValueError("Directory doesn't exist: "+ self.directory)

        if not self.directory.endswith('/'):
            self.directory += '/'

        listdir = glob.glob(self.directory + self.file_pattern)
        listdir.sort()

        #Set up variable to hold the date of the last day 1 in all the files
        #(note as the files are checked in reverse order, this
        # will be in the first file found)
        #This will be used for limiting files for forecast_day='latest'
        last_day1_date = None

        #Loop through files in reverse order,
        # to make finding latest files easier
        for filename in reversed(listdir):

            if not filename.endswith('.pp'):
                #Not a pp file - discard
                continue

            fileinfo = self.get_fileinfo(filename)

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

                #Discard file if file does not contain required forecast_day
                if fileinfo['forecast_days'] is not None:

                    if self.forecast_day == 'forecast':
                        if (fileinfo['day1_date'].date()
                                != self.start_datetime.date()):
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
                            last_day1_date = fileinfo['day1_date'].date()
                        if not (1 in fileinfo['forecast_days'] or \
                           last_day1_date == fileinfo['day1_date'].date()):
                            #If the forecast_day is not day 1, or if the day1
                            #date of the file is not the last_day1_date, then
                            #ignore this file.
                            continue
                    elif self.forecast_day == 'all':
                        #Keep all files
                        pass
                    else:
                        if (int(self.forecast_day) not in
                                fileinfo['forecast_days']):
                            continue

            if self.filenames is None:
                self.filenames = []
            self.filenames.append(filename)

            if self.runtimes is None:
                self.runtimes = []

            if fileinfo['runtime'] not in self.runtimes:
                self.runtimes.append(fileinfo['runtime'])

        if self.filenames is None:
            raise UserWarning("No files found with required data")

        #Now reverse back to get in sorted order again
        self.filenames.sort()

        return self.filenames



    def get_fileinfo(self, filename):
        """

        Returns dictionary of information gathered from the filename:
          * day1_date - datetime format which corresponds to the
                                 date of the first full day of data from
                                 this model forecast run.
          * data_start - datetime format of the start time
          * data_end - datetime format of the end time
          * filename - filename
          * forecast_days - list of integer day which data corresponds to.
            Note day 1 is first full day.
          * forecast_period - list of expected leadtimes in file (note leadtimes
                             correspond to the end of the period for means)
          * runtime - model runtime

        >>> filename = 'aqum_20130801_T+012_QM18.pp'
        >>> AQFiles = AQUMppFiles()
        >>> fileinfo = AQFiles.get_fileinfo(filename)
        >>> fileinfo == {'forecast_days': [1],
        ... 'runtime': 18,
        ... 'day1_date': datetime.datetime(2013, 8, 2, 0, 0),
        ... 'data_start': datetime.datetime(2013, 8, 2, 6, 0),
        ... 'data_end': datetime.datetime(2013, 8, 2, 18, 0),
        ... 'filename': 'aqum_20130801_T+012_QM18.pp',
        ... 'forecast_periods':
        ... [13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25] }
        True

        .. note::

            Assumes files in following format:

            prefix_yyyymmdd.pp (contains 24 hours of data)

            OR

            prefix_yyyymmdd_T+HHH_QMHH.pp (contains 12 hours of data)

            OR

            prods_op_aqum_20110701_00.000.pp (contains 12 hours of data)


        """

        fileinfo = {}
        splitfile = os.path.basename(filename).split('_')

        if len(splitfile) == 4 or len(splitfile) == 5:
            if len(splitfile) == 4:
                #File in format aqum_20130806_T+024_QM18.pp
                dtstr = splitfile[1][0:8]
                leadtime = int(splitfile[2][2:])
                runtime = int(splitfile[3][2:4])
            elif len(splitfile) == 5:
                #File in format prods_op_aqum_20110701_00.000.pp
                dtstr = splitfile[3]
                splitend = splitfile[4].split('.')
                leadtime = int(splitend[1])
                runtime = int(splitend[0])
            hours = 12
            hours_to_0z = (24 - runtime) % 24
            if hours_to_0z == 12:
                #For consistency with old casestudies,
                # use day1 as 13Z current day -12Z next day
                forecast_days = [leadtime//24]
            else:
                #Usually do 01Z-24Z from first full day
                if runtime > 0:
                    #First full day starts from 00Z on date after rundate
                    forecast_days = list(set([(runtime+leadtime)//24,
                                              (runtime+leadtime+hours-1)//24]))
                else:
                    #First full day starts from 00Z on rundate.
                    forecast_days = list(set([
                        (runtime+leadtime+24)//24,
                        (runtime+leadtime+24+hours-1)//24]))
        else:
            #File in format aqum_20060806.pp - assumes from 12Z run
            dtstr = splitfile[1][0:8]
            leadtime = 00
            runtime = 12
            hours = 24
            forecast_days = [1]
            hours_to_0z = 12

        file_date = datetime.datetime.strptime(dtstr+"00", "%Y%m%d%H")
        data_start_datetime = file_date + \
                              datetime.timedelta(hours=runtime+leadtime)
        data_end_datetime = data_start_datetime + \
                            datetime.timedelta(hours=hours)
        forecast_periods = list(np.arange(leadtime+1, leadtime+hours+2))

        day1_date = file_date + \
                    datetime.timedelta(hours=runtime) + \
                    datetime.timedelta(hours=hours_to_0z)

        fileinfo['day1_date'] = day1_date
        fileinfo['data_start'] = data_start_datetime
        fileinfo['data_end'] = data_end_datetime
        fileinfo['filename'] = filename
        fileinfo['forecast_days'] = forecast_days
        fileinfo['forecast_periods'] = forecast_periods
        fileinfo['runtime'] = runtime

        return fileinfo


def __mass_retrieve_global_str(short_name_list=None, forecast_day=None,
                               runtime=18, filenames="*prods*",
                               start_datetime=None):
    """
    Set up string containing global filter ready for moo select retrievals.

    :param short_name_list: list of short names that need to be retrieved
    :param forecast_day: This can be an integer, whose number corresponds
                         to the forecast day, eg 1 if using the first full
                         forecast day from each model run.
                         Alternatively, set to "forecast" to get all leadtimes
                         from a single model run.
    :param runtime: runtime of model. Usually 18 (default) or 0.
    :param filenames: String. File(s) to be retrieved, possibly by using
                      wildcards
    :returns: String which can be output to file for moo select filter.

    >>> global_str = __mass_retrieve_global_str(
    ... short_name_list=['O3','PM2p5','T_1p5','Ox'],
    ... forecast_day=1, runtime=18, filenames='*prods*')
    >>> print(global_str)
    begin_global
      stash=(03236,17221,34001,34004)
      lbft=(7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,\
30)
      pp_file="*prods*"
    end_global
    <BLANKLINE>
    >>> global_str = __mass_retrieve_global_str(
    ... short_name_list=['O3','PM2p5','T_1p5'],
    ... forecast_day='forecast', runtime=18, filenames='*prods*',
    ... start_datetime=datetime.datetime(2014,3,27))
    >>> print(global_str)
    begin_global
      stash=(03236,17221,34001)
      pp_file="*prods*20140326*"
    end_global
    <BLANKLINE>

    """

    global_str = "begin_global\n"

    #Add required stash codes
    if short_name_list is not None:
        global_str += '  stash='
        stashcodes = set()
        for short_name in short_name_list:
            stash = None
            if short_name in _SHORT_NAME_2_STASH:
                stash = _SHORT_NAME_2_STASH[short_name]
                #Convert from m12s34i567 to 34567
                stash = stash[4:6]+stash[7:10]
                stashcodes.add(stash)
            else:
                #Try and identify it as a derived short name instead
                spec = adaq_data.get_derived_sn_specs(
                    [short_name],
                    derived_short_name_dict=_DERIVED_SHORT_NAME_DICT)[0]
                for component in spec.get('components', []):
                    stash = _SHORT_NAME_2_STASH[component]
                    #Convert from m12s34i567 to 34567
                    stash = stash[4:6]+stash[7:10]
                    stashcodes.add(stash)
            if stash is None:
                #No stash codes were added
                raise ValueError("Stash code not found for "+short_name)
        if len(stashcodes) > 1:
            #Need brackets around stash list only
            #if more than one element in list.
            global_str += '('
        global_str += ",".join(sorted(stashcodes))
        if len(stashcodes) > 1:
            global_str += ')'
        global_str += '\n'

    #Limit by forecast_periods (leadtime)
    #NB only done for specific integer days.
    #For forecast_day='forecast', need to limit by filenames,
    #which are based on the required dates.
    if forecast_day is not None:
        if forecast_day == 'forecast':
            rundate = start_datetime.date()
            if int(runtime) > 12:
                #rundate is actually the previous day
                rundate = start_datetime.date() - datetime.timedelta(days=1)
            filenames += rundate.strftime("%Y%m%d")+'*'

        elif forecast_day == 'latest' or forecast_day == 'all':
            #For latest forecast day, or all forecast days,
            #just retrieve everything as difficult to
            #set up anything else for mass retrievals as don't know how long the
            #forecast runs are for.
            pass
        else:
            #Integer forecast day
            forecast_periods = array_statistics.calc_forecast_periods(
                int(forecast_day), int(runtime))
            if forecast_periods is not None:
                forecast_periods_str = [str(fp) for fp in forecast_periods]
                global_str += "  lbft=("+",".join(forecast_periods_str)+')\n'

    #Limit by filenames on mass
    if filenames is not None:
        global_str += '  pp_file="'+filenames+'"\n'

    global_str += 'end_global\n'

    return global_str

def __mass_retrieve_meaninst_str(start_dt=None, end_dt=None):
    """
    Set up two strings, containing moo select filters for meaned
    pp data and instantaneous data

    :param start_dt: datetime format of start time to retrieve
    :param end_dt: datetime format of end time to retrieve
    :returns: (mean_str, instant_str) - tuple of strings, first one
              containing moo select filter commands for meaned data,
              the second for instantaneous data.

    >>> mean_str, instant_str = __mass_retrieve_meaninst_str(
    ... start_dt=datetime.datetime(2014,11,5),
    ... end_dt=datetime.datetime(2014,11,9))
    >>> print(mean_str)
    begin
      lbproc=128
      T2>={2014/11/05 00:00:00}
      T2<={2014/11/09 00:00:00}
    end
    <BLANKLINE>

    >>> print(instant_str)
    begin
      lbproc=0
      T1>={2014/11/05 00:00:00}
      T1<={2014/11/09 00:00:00}
    end
    <BLANKLINE>

    """

    mean_str = "begin\n"
    instant_str = "begin\n"
    mean_str += '  lbproc=128\n' #Meaned data
    instant_str += '  lbproc=0\n' #Instantaneous data

    #Start and end dates...
    #For meaned fields, the actual time point is at the
    # end of the meaning period.
    # For this, use T2, which is the equivalent of
    # lbyrd, lbmond, lbdatd, lbmind in pp header.
    #For instantaneous data, the time is given by T1, which is
    # lbyr, lbmon, lbdat, lbhr, lbmin in pp header.

    if start_dt is not None:
        start_str = start_dt.strftime("%Y/%m/%d %H:%M:%S")
        #T1: validity date/time of field
        #    For meaned fields, this is the start of the meaning period
        mean_str += "  T2>={"+start_str+"}\n"
        instant_str += "  T1>={"+start_str+"}\n"

    if end_dt is not None:
        end_str = end_dt.strftime("%Y/%m/%d %H:%M:%S")
        #T2: For meaned fields, this is the end of the meaning period
        mean_str += "  T2<={"+end_str+"}\n"
        #T1: validity date/time of field
        instant_str += "  T1<={"+end_str+"}\n"

    mean_str += 'end\n'
    instant_str += 'end\n'

    return mean_str, instant_str


def pp_mass_retrieve(outputdir, runid=None, short_name_list=None,
                     start_datetime=None, end_datetime=None,
                     filenames="*prods*", operational=False, psuite=False,
                     forecast_day=None, runtime=18,
                     massdir=None, massretries=0,
                     massretrydelay=60, retrieve=True):
    """
    Function which will set up and perform mass retrieval for pp files.

    :param outputdir: String. Directory to retrieve files into.
                      Will be created if does not already exist.
    :param runid: String. Rose suite name.
    :param short_name_list: list of short names that need to be retrieved
    :param start_datetime: datetime format of start time from files
    :param end_datetime: datetime format of end time from files
    :param filenames: String. File(s) to be retrieved, possibly by using
                      wildcards
    :param operational: Logical or string. Set to True, or 'aqum' or 'aqeur'
                        to retrieve data from operational/parallel suite AQEUR
                        archive. Set to 'aqcops' to retrieve data from
                        operational/parallel suite aqcops archive. False will
                        retrieve from :/devfc instead.
    :param psuite: Set to False by default or set to string such as 'psuite37'
                   to get parallel suite data from the specified operational
                   archive. The psuite parameter will be ignored if operational
                   is set to False.
    :param forecast_day: This can be an integer, whose number corresponds
                         to the forecast day, eg 1 if using the first full
                         forecast day from each model run.
                         Alternatively, set to "forecast" to get all leadtimes
                         from a single model run.
                         Or set to "latest" to get the latest available
                         data, which would be day 1 forecasts where
                         available and then the forecast from a single run
                         at the end.
    :param runtime: runtime of model. Usually 18 (default) or 0.
    :param massdir: String. Directory on mass to retrieve from.
    :param massretries: number of times to retry mass retrieval
    :param massretrydelay: sleep time in seconds between each retry
    :param retrieve: Logical. If set to True retrieves files from mass.
                     If False, sets up all required moo select files and returns
                     the command to be run manually by user.
    :returns: **(mooselect_strs, commands)** where:

                * **mooselect_strs** - list (one element per retrieval)
                  containing the strings passed into the moo select file.
                * **commands** - list (one element per retrieval),
                  containing the unix commands to be run to perform
                  the retrieval.

    .. Note:: This code is only guaranteed to work with files archived by
              rose suite

    This code works by setting up a query file to be used by the moo
    select command.
    This enables a limited set of data to be retrieved from pp files on mass,
    hence saving on local disk space by only retrieving exactly the
    required data and also speeds up the time taken for the retrieval.
    Further information on this can be found on the
    `Moose User Guide pages on record level retrieval query syntax
    <http://www01/teams/mass_replacement/docs/user_guide.html#
    record-level-retrieval-query-syntax>`_

    Example of calling, although here retrieve is set to False.

    >>> mooselect_strs, commands = pp_mass_retrieve(
    ... '$DATADIR/mass_retrieve', runid='mi-ah183',
    ... filenames='*prods*', short_name_list=['O3','PM2p5','T_1p5'],
    ... start_datetime=datetime.datetime(2014,11,5),
    ... end_datetime=datetime.datetime(2014,11,9),
    ... forecast_day=1, runtime=18,
    ... operational=False, massretries=1, massretrydelay=5,
    ... retrieve=False) # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    Retrieving for period  2014-11-05 00:00:00 2014-11-09 00:00:00
    Moose Retrieval Command:
    moo select -f .../mass_retrieve/mooselect_20141105_20141109.txt \
    moose:/devfc/mi-ah183/field.pp .../mass_retrieve

    >>> print(mooselect_strs[0])
    begin_global
      stash=(03236,17221,34001)
      lbft=(7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,\
30)
      pp_file="*prods*"
    end_global
    begin
      lbproc=128
      T2>={2014/11/05 00:00:00}
      T2<={2014/11/09 00:00:00}
    end
    begin
      lbproc=0
      T1>={2014/11/05 00:00:00}
      T1<={2014/11/09 00:00:00}
    end
    <BLANKLINE>

    >>> print(commands) # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    ['moo select -f .../mass_retrieve/mooselect_20141105_20141109.txt
    moose:/devfc/mi-ah183/field.pp .../mass_retrieve']


    Example of retrieving parallel suite data

    >>> mooselect_strs, commands = pp_mass_retrieve(
    ... '$DATADIR/mass_retrieve', runid='para',
    ... filenames='*prods*', short_name_list=['O3','PM2p5','T_1p5'],
    ... start_datetime=datetime.datetime(2016,3,11),
    ... end_datetime=datetime.datetime(2016,3,12),
    ... forecast_day=1, runtime=18,
    ... operational='aqum', psuite='psuite37', massretries=1, massretrydelay=5,
    ... retrieve=False) # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    Retrieving for period  2016-03-11 00:00:00 2016-03-12 00:00:00
    Moose Retrieval Command:
    moo select -f .../mass_retrieve/mooselect_20160311_20160312.txt \
    moose:/opfc/atm/aqum/prods/psuite37.pp/ .../mass_retrieve

    """

    moo_commands = [] #Will hold list of shell commands
    mooselect_strs = [] #Will hold list of contents of query files

    #Make output directory if does not already exist
    outputdir = os.path.expandvars(outputdir)
    if not os.path.isdir(outputdir):
        os.makedirs(outputdir)

    #----
    #Set up moo select query file
    #Global attributes: stash, forecast_period, pp filename.
    global_str = __mass_retrieve_global_str(short_name_list, forecast_day,
                                            runtime, filenames, start_datetime)

    #---
    #From "moo si -l"
    #Multiple-get file-number limit: 1000
    #=> Can't retrieve more that 1000 files at once.
    #=> Assuming 6hrly files (eg for prodm) => 250 days
    #=> limit to 200 days per retrieval (>6months)
    #=> Need to loop over date-specific begin/end statements
    #   and do multiple retrievals
    retrievals_limit = 200 #days

    #Generate list of (start_dt,end_dt) tuples
    #that are required to retrieve separately.
    start_dt = start_datetime
    end_dt = start_dt
    start_end_list = []
    while end_dt <= end_datetime:
        end_dt += datetime.timedelta(days=retrievals_limit)
        if end_dt > end_datetime:
            #Don't go beyond required final date
            end_dt = end_datetime
        #Add this start/end date tuple to list
        start_end_list.append((start_dt, end_dt))
        if end_dt >= end_datetime:
            #Exit loop as end_datetime has already been included
            break
        start_dt = end_dt


    #Now loop over all retrieval periods:
    #Set up date-based mean and instantanteous strings
    #Write these and global strings to moo select query file
    #Set up directory(s) to retrieve from on mass
    #Perform mass retrieval
    for start_dt, end_dt in start_end_list:
        print("Retrieving for period ", start_dt, end_dt)

        #---
        #Set up moo select query strings for Meaned/Instantaneous data
        mean_str, instant_str = __mass_retrieve_meaninst_str(
            start_dt=start_dt, end_dt=end_dt)

        #---
        #Write query strings to moo select query file
        mooselect_str = global_str + mean_str + instant_str
        mooselect_strs.append(mooselect_str) #For returning to user

        selectfile = outputdir+'/mooselect.txt'
        if start_datetime is not None and end_datetime is not None:
            selectfile = outputdir+'/mooselect_' + \
                         start_dt.strftime("%Y%m%d")+'_' + \
                         end_dt.strftime("%Y%m%d")+'.txt'
        with open(selectfile, 'w') as fout:
            fout.write(mooselect_str)

        #---
        #Set up directory to retrieve data from on mass
        if massdir is None:
            op_mass_dir = None
            # If operational is a boolean value, then either
            # aqeur operational/parallel suite or user suite.
            if isinstance(operational, bool):
                if operational:
                    if filenames == '*prodm*':
                        op_mass_dir = 'moose:/opfc/atm/aqum/prodm/'
                    else:
                        op_mass_dir = 'moose:/opfc/atm/aqum/prods/'
            else:
                #string set to 'aqum', 'aqeur' or 'aqcops'
                if operational == 'aqum' or operational == 'aqeur':
                    if filenames == '*prodm*':
                        op_mass_dir = 'moose:/opfc/atm/aqum/prodm/'
                    else:
                        op_mass_dir = 'moose:/opfc/atm/aqum/prods/'
                elif operational == 'aqcops':
                    if filenames == '*prodm*':
                        op_mass_dir = 'moose:/opfc/atm/aqcops/prodm/'
                    else:
                        op_mass_dir = 'moose:/opfc/atm/aqcops/prods/'
                else:
                    raise ValueError('operational set incorrectly')
            if op_mass_dir is not None:
                #This has been set because operational or psuite
                if not isinstance(psuite, bool):
                    massdir = op_mass_dir + psuite + '.pp/'
                else:
                    #Start from 5 days before start_datetime to catch instances
                    #where need to retrieve forecast data from previous year's
                    #directory
                    years = range(
                        (start_datetime-datetime.timedelta(days=5)).year,
                        end_datetime.year+1)
                    massdir = ""
                    for year in years:
                        massdir_year = op_mass_dir + '%04d.pp ' % (year)
                        if len(years) > 1:
                            #For multiple years, check mass directory exists first
                            exists = shell_commands.check_mass_dir(massdir_year)
                            if not exists:
                                print('Mass dir does not exist:', massdir_year)
                            else:
                                #True or Unknown (possibly from system error)
                                massdir += massdir_year
                        else:
                            #Only one year
                            massdir = massdir_year
                    if not massdir:
                        raise ValueError('No existing mass directories found')
            else:
                #user suite
                massdir = 'moose:/devfc/'+runid+'/field.pp'
        #---
        #Do retrieval
        cmd = shell_commands.call_mass(selectfile, massdir, outputdir,
                                       massretries=massretries,
                                       massretrydelay=massretrydelay,
                                       retrieve=retrieve)
        moo_commands.append(cmd)

    return mooselect_strs, moo_commands


if __name__ == "__main__":

    import doctest
    doctest.testmod()
