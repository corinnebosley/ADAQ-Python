"""
Class for reading MACC ensemble netcdf files into an ADAQData class
(Note for reading grib-format MACC ensemble data, use ecgrib_data.py instead)
"""
from __future__ import division
from __future__ import absolute_import
from __future__ import print_function

from six.moves.builtins import range

import datetime
import os

import cf_units
import iris

import adaq_data
import array_statistics
import cube_time

MACC_2_STDNAME = {
    'o3_conc' : 'mass_concentration_of_ozone_in_air',
    'no2_conc' : 'mass_concentration_of_nitrogen_dioxide_in_air',
    'so2_conc' : 'mass_concentration_of_sulfur_dioxide_in_air',
    'co_conc' : 'mass_concentration_of_carbon_monoxide_in_air',
    'pm10_conc' : 'mass_concentration_of_pm10_ambient_aerosol_in_air',
    'pm2p5_conc' : 'mass_concentration_of_pm2p5_ambient_aerosol_in_air'
    }

MACC_2_SHORTNAME = {
    #First set for use with MACC_V2014 file formats:
    'o3_conc' : 'O3',
    'no2_conc' : 'NO2',
    'so2_conc' : 'SO2',
    'co_conc' : 'CO',
    'pm10_conc' : 'PM10',
    'pm2p5_conc' : 'PM2p5',
    #This set for use with CAMS file formats:
    'mass_concentration_of_ozone_in_air' : 'O3',
    'mass_concentration_of_nitrogen_dioxide_in_air' : 'NO2',
    'mass_concentration_of_sulfur_dioxide_in_air' : 'SO2',
    'mass_concentration_of_carbon_monoxide_in_air': 'CO',
    'mass_concentration_of_pm10_ambient_aerosol_in_air': 'PM10',
    'mass_concentration_of_pm2p5_ambient_aerosol_in_air': 'PM2p5',
    'Birch_Pollen_Grain' : 'BIRCH_POLLEN'
    }

class MaccEnsData(adaq_data.ADAQData):
    """
    Subclass of ADAQData, which contains extra functionality specific to
    MACC Ensemble Data (in netcdf format).
    Unique to this data format is the time is given as hours since a date
    as defined in the long_name for time. This class deals with this and
    converts it to a date-time format.

    **Example:**

    >>> import config
    >>> directory = config.SAMPLE_DATADIR+'macc_ens/CAMS/'

    Initialise class:

    >>> maccens = MaccEnsData()

    Determine which filenames should be read in according to required
    start and end dates:

    >>> filenames = maccens.get_filenames(directory,
    ...     start_datetime=datetime.datetime(2016, 6, 15),
    ...     end_datetime=datetime.datetime(2016, 6, 16),
    ...     forecast_day=1)
    >>> print(filenames[0]) #doctest: +ELLIPSIS
    /.../ENSEMBLE_FORECAST_SURFACE_BIRCH_POLLEN_0H_24H_201606140000000000.nc

    Now read the data into a gridded_cube_list. Note filenames,
    start_datetime, end_datetime and forecast_day have been saved
    in the maccens object so no need to give these again.
    This is 3D data, but here we only need to consider the lowest model
    level which is at 0m, so set level=0.
    Also limit which short_names, model and forecast_day to load:

    >>> gcl = maccens.readdata(level=0, short_name_list=['PM2p5','O3'],
    ... model='ENSEMBLE', forecast_day=1)

    Examine data:

    >>> print(maccens.gridded_cube_list)
    0: mass_concentration_of_ozone_in_air / (ug/m3) \
(time: 25; latitude: 400; longitude: 700)
    1: mass_concentration_of_pm2p5_ambient_aerosol_in_air / (ug/m3) \
(time: 25; latitude: 400; longitude: 700)

    Have a look at the ozone gridded cube:

    >>> o3 = maccens.extract(short_name='O3',gridded=True,singlecube=True)
    >>> print(o3)
    mass_concentration_of_ozone_in_air / (ug/m3) (time: 25; latitude: 400; longitude: 700)
         Dimension coordinates:
              time                                    x             -               -
              latitude                                -             x               -
              longitude                               -             -               x
         Auxiliary coordinates:
              forecast_day                            x             -               -
              forecast_period                         x             -               -
         Scalar coordinates:
              level: 0.0 m
         Attributes:
              Model: ENSEMBLE
              comment: Tools used: CDO V1.6.1rc6 and NCO V4.3.7
              institution: Data produced by Meteo France
              label: ENSEMBLE
              nco_openmp_thread_number: 1
              project: MACC-RAQ (http://macc-raq.gmes-atmosphere.eu)
              short_name: O3
              source: Data from ENSEMBLE model
              species: Ozone
              title: O3 Air Pollutant FORECAST at the Surface
              value: hourly values

    Can now extract data at specific sites. Firstly get some sites data:

    >>> import sites_info
    >>> sitesfilename = config.SAMPLE_DATADIR + \
    'AURN_obs/aq_sites_GEMSRAQ_v4b_dev.txt'
    >>> obsdir = config.SAMPLE_DATADIR+'AURN_obs/'
    >>> si = sites_info.SitesInfo()

    >>> sites_data = si.read_from_file(sitesfilename,
    ... allsites=False, obsdir=obsdir)
    Number of sites:  5

    Now extract these sites from MACC ensemble data, into a sites_cube_list:

    >>> scl = maccens.extract_sites(sites_data=sites_data)

    Have a look at the PM2p5 cube:

    >>> pm2p5_sites_cube = maccens.extract(short_name='PM2p5',singlecube=True)
    >>> print(pm2p5_sites_cube)
    mass_concentration_of_pm2p5_ambient_aerosol_in_air / (ug/m3) (site_id: 5; time: 25)
         Dimension coordinates:
              site_id                                                    x        -
              time                                                       -        x
         Auxiliary coordinates:
              abbrev                                                     x        -
              latitude                                                   x        -
              longitude                                                  x        -
              site_altitude                                              x        -
              site_name                                                  x        -
              site_type                                                  x        -
              forecast_day                                               -        x
              forecast_period                                            -        x
         Scalar coordinates:
              level: 0.0 m
         Attributes:
              Model: ENSEMBLE
              comment: Tools used: CDO V1.6.1rc6 and NCO V4.3.7
              institution: Data produced by Meteo France
              label: ENSEMBLE
              nco_openmp_thread_number: 1
              project: MACC-RAQ (http://macc-raq.gmes-atmosphere.eu)
              short_name: PM2p5
              source: Data from ENSEMBLE model
              species: PM2.5 Aerosol
              title: PM25 Air Pollutant FORECAST at the Surface
              value: hourly values

    Check sensible values:

    >>> print('{:.3f}'.format(o3.data.max()))
    142.360
    >>> print('{:.3f}'.format(pm2p5_sites_cube.data.mean()))
    5.021

    This MaccEns class can also be used to look at the older format data,
    (from MACC_V2014 directory on ftp site).
    Note this format has different filenames
    and slightly different attributes in the loaded cubes.

    >>> directory = config.SAMPLE_DATADIR+'macc_ens/'
    >>> maccens = MaccEnsData()

    >>> filenames = maccens.get_filenames(directory,
    ...     start_datetime=datetime.datetime(2014, 3, 27),
    ...     end_datetime=datetime.datetime(2014, 3, 28),
    ...     forecast_day=1)
    >>> print(filenames[0]) #doctest: +ELLIPSIS
    /.../macc_ens/HRES_ENS_2014032600+024.nc
    >>> gcl = maccens.readdata(level=0)
    >>> print(maccens.gridded_cube_list)
    0: mass_concentration_of_carbon_monoxide_in_air / (ug/m3) \
(time: 25; latitude: 400; longitude: 700)
    1: mass_concentration_of_nitrogen_dioxide_in_air / (ug/m3) \
(time: 25; latitude: 400; longitude: 700)
    2: mass_concentration_of_ozone_in_air / (ug/m3) \
(time: 25; latitude: 400; longitude: 700)
    3: mass_concentration_of_pm10_ambient_aerosol_in_air / (ug/m3) \
(time: 25; latitude: 400; longitude: 700)
    4: mass_concentration_of_pm2p5_ambient_aerosol_in_air / (ug/m3) \
(time: 25; latitude: 400; longitude: 700)
    5: mass_concentration_of_sulfur_dioxide_in_air / (ug/m3) \
(time: 25; latitude: 400; longitude: 700)
    >>> o3 = maccens.extract(short_name='O3',gridded=True,singlecube=True)
    >>> print(o3)
    mass_concentration_of_ozone_in_air / (ug/m3) (time: 25; latitude: 400; longitude: 700)
         Dimension coordinates:
              time                                    x             -               -
              latitude                                -             x               -
              longitude                               -             -               x
         Auxiliary coordinates:
              forecast_period                         x             -               -
         Scalar coordinates:
              forecast_day: 1 days
              level: 0.0 m
         Attributes:
              Content: Air Pollutant FORECAST at Levels:     0.     50.    250.    500.   1000....
              Model: ENS
              Service: MACC RAQ Ensemble
              Source: MACC RAQ
              Title: MACC Air Quality FORECAST
              label: ENS
              nco_openmp_thread_number: 1
              short_name: O3
              species: Ozone
              value: hourly values

    It can also be used to read in data from the new Data Server files in netcdf
    format:

    >>> directory = config.SAMPLE_DATADIR+'macc_ens/CAMS_DataServer'
    >>> maccens = MaccEnsData()
    >>> filenames = maccens.get_filenames(directory,
    ...     start_datetime=datetime.datetime(2016, 9, 14),
    ...     end_datetime=datetime.datetime(2016, 9, 16),
    ...     forecast_day=1, short_name_list=['O3'])
    >>> print(filenames[0]) #doctest: +ELLIPSIS
    /.../W_fr-meteofrance,MODEL,ENSEMBLE+FORECAST+SURFACE+O3+\
0H24H_C_LFPW_20160914000000.nc

    >>> gcl = maccens.readdata(level=0)
    >>> print(maccens.gridded_cube_list[0])
    mass_concentration_of_ozone_in_air / (ug/m3) (time: 48; latitude: 400; longitude: 700)
         Dimension coordinates:
              time                                    x             -               -
              latitude                                -             x               -
              longitude                               -             -               x
         Auxiliary coordinates:
              forecast_day                            x             -               -
              forecast_period                         x             -               -
         Scalar coordinates:
              level: 0.0 m
         Attributes:
              Model: ENSEMBLE
              comment: Tools used: CDO V1.6.1rc6 and NCO V4.3.7
              institution: Data produced by Meteo France
              label: ENSEMBLE
              nco_openmp_thread_number: 1
              project: MACC-RAQ (http://macc-raq.gmes-atmosphere.eu)
              short_name: O3
              source: Data from ENSEMBLE model
              species: Ozone
              title: O3 Air Pollutant FORECAST at the Surface
              value: hourly values

    To read data all from a  single forecast run, set the keyword
    forecast_day='forecast'. For more information about this setting, see
    :func:`cube_time.extract_latest_forecast_days`.
    Note can also pass the short_name_list and required model into get_filenames
    to limit the files that have to be considered on loading.

    Alternatively forecast_day can be set to 'latest'. This loads day 1
    forecasts where possible, expanding out to the full forecast range
    on the final run if needed. The routine get_filenames can be used to
    limit the required files down, or the readdata routine can be called
    directly. For more information about this
    setting, see :func:`cube_time.extract_latest_forecast_days`.

    To read in files with all available forecast days and then force each
    different forecast day into a different cube, use the keyword
    forecast_day='all'. For more information about this setting see
    :func:`cube_time.extract_all_forecast_days`.

    """

    def __init__(self, label=None):
        """
        Initiates a class from MACC Ensemble data as a subclass of ADAQData.
        """

        adaq_data.ADAQData.__init__(self)

        self.label = label #Label (else set from info in MACC file)
        self.short_name_list = None #List of short names
        self.start_datetime = None #Starting datetime
        self.end_datetime = None #End datetime
        self.level = None #Vertical level
        self.sites_data = None #sites_data object
        self.filenames = None #List of filenames
        self.forecast_day = None #List of forecast days
        self.model = None #Which model to read in

        #Set up itime - this will be a counter which will be used in callbacks
        #to ensure that the temporary itime coordinate is strictly monotonic,
        #but will be incremented everytime the callback is called.
        self.itime = 0

    def readdata(self, filenames=None, short_name_list=None, level=None,
                 start_datetime=None, end_datetime=None, label=None,
                 model=None, forecast_day=None):
        """
        Create a gridded_cube_list from filenames specified.

        :param filenames: list of files to read, can include wildcards etc
        :param short_name_list: list of short names (None gets everything)
        :param level: model level height - eg 0 for the surface at 0m.
        :param start_datetime: datetime object start date
        :param end_datetime: datetime object end date
        :param label: label for the cubes
        :param model: model name to constrain by
        :param forecast_day: This can be an integer, whose number corresponds
                             to the forecast day, eg 1 if using the first full
                             forecast day from each model run.
                             Alternatively, set to "forecast" to get all
                             leadtimes from a single model run.
                             Or set to "latest" to get the latest available
                             data, which would be day 1 forecasts where
                             available and then the forecast from a single run
                             at the end.
        """

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
        if model is not None:
            self.model = model

        if self.filenames is None:
            raise ValueError('No filenames given')

        constraints = None
        time_constraints = None

        if self.short_name_list is not None:
            sname_constraint = iris.AttributeConstraint(
                short_name=lambda c: c in self.short_name_list)

            constraints = constraints & sname_constraint

        if self.model is not None:
            model_constraint = iris.AttributeConstraint(Model=self.model)
            constraints = constraints & model_constraint

        if self.level is not None:
            level_constraint = iris.Constraint(level=self.level)
            constraints = constraints & level_constraint

        if not (self.forecast_day is None or \
                self.forecast_day == 'forecast' or \
                self.forecast_day == 'latest' or \
                self.forecast_day == 'all'):

            day_constraint = iris.Constraint(
                forecast_day=int(self.forecast_day))
            constraints = constraints & day_constraint

        if self.start_datetime is not None:
            time_constraint = iris.Constraint(
                time=lambda c: c.point >= self.start_datetime)

            time_constraints = time_constraints & time_constraint

        if self.end_datetime is not None:
            time_constraint = iris.Constraint(
                time=lambda c: c.point <= self.end_datetime)

            time_constraints = time_constraints & time_constraint


        #Load data into gridded_cube_list
        self.gridded_cube_list = iris.load(self.filenames,
                                           constraints,
                                           callback=self.__callback)


        #Join time-dimensions to get single cubes for each short name
        self.gridded_cube_list = self.gridded_cube_list.concatenate()

        #Extract times
        self.gridded_cube_list = self.gridded_cube_list.extract(
            time_constraints)

        #Extract forecast days required through the setting of forecast_day.
        if self.forecast_day == 'latest' or self.forecast_day == 'forecast':
            for i, cube in enumerate(self.gridded_cube_list):
                cube.remove_coord('itime')
                self.gridded_cube_list[i] = \
                    cube_time.extract_latest_forecast_days(
                        cube, self.forecast_day, start_dt=self.start_datetime)
        elif self.forecast_day == 'all':
            self.gridded_cube_list = cube_time.extract_all_forecast_days(
                self.gridded_cube_list)

        #Add coord sys to lat and lon coords
        lat_lon_coord_sys = iris.coord_systems.GeogCS(
            semi_major_axis=iris.fileformats.pp.EARTH_RADIUS)

        for cube in self.gridded_cube_list:
            if cube.coord('longitude'):
                if not cube.coord('longitude').coord_system:
                    cube.coord('longitude').coord_system = lat_lon_coord_sys
            if cube.coord('latitude'):
                if not cube.coord('latitude').coord_system:
                    cube.coord('latitude').coord_system = lat_lon_coord_sys
            #Add bounds to time coordinate
            if cube.coord('time'):
                if len(cube.coord('time').points) > 1:
                    if not cube.coord('time').has_bounds():
                        #Assume points are given at the end of the
                        #meaning period
                        cube.coord('time').guess_bounds(bound_position=1.0)

        return self.gridded_cube_list

    def __callback(self, cube, field, filename):
        """
        Private method - an iris callback function.
        This removes the history attribute (otherwise cubes may not merge).
        Standard name and short name and label attributes are added.
        Time coordinates are corrected to be in date-time format rather
        than hours since forecast time.
        """

        #Remove any attributes that may cause problems when
        #trying to merge cubes from multiple days
        if 'history' in cube.attributes:
            del cube.attributes['history']
        if 'FORECAST' in cube.attributes:
            del cube.attributes['FORECAST']
        if 'summary' in cube.attributes:
            del cube.attributes['summary']

        cubename = cube.name()
        if cubename == 'Not Defined' and 'species' in cube.attributes:
            cubename = cube.attributes['species'].replace(' ', '_')
            cube.rename(cubename)
        if cubename in MACC_2_SHORTNAME:
            cube.attributes['short_name'] = MACC_2_SHORTNAME[cubename]
        if cubename in MACC_2_STDNAME:
            cube.rename(MACC_2_STDNAME[cubename])

        if 'source' in cube.attributes:
            cube.attributes['Model'] = cube.attributes['source'].split(' ')[2]
        if self.label is None:
            cube.attributes['label'] = cube.attributes['Model']
        else:
            cube.attributes['label'] = self.label


        cube = self.correct_timecoord(cube)

        #Fix units as u from ug/m3 get confused on loading
        #NB, still issues warning
        if field.units == u'\xb5g/m3':
            cube.units = 'ug/m3'
            if 'invalid_units' in cube.attributes:
                del cube.attributes['invalid_units']

        #Sort out longitudes to be from -180-+180 to ensure montonically
        #increasing so can be a dimension coordinate rather than auxillary
        if isinstance(cube.coord('longitude'), iris.coords.AuxCoord):
            mask = cube.coord('longitude').points > 180.
            cube.coord('longitude').points[mask] -= 360.
            iris.util.promote_aux_coord_to_dim_coord(cube, 'longitude')



    def correct_timecoord(self, cube):
        """
        Convert time coordinate to gregorian calendar, instead of hours, with
        reference time given in long_name.
        Also gives forecast_period coordinate and ensures time coord has bounds.
        Also adds forecast_day coordinate.
        """

        time_unit = cf_units.Unit('hours since epoch', calendar='gregorian')
        tcoord = cube.coord('time')

        #Get reference time, taken from long_name which is
        #eg u'FORECAST time from 2015021900'
        assert tcoord.long_name[:18] == 'FORECAST time from'
        timestr = tcoord.long_name.split()[-1]
        if len(timestr) == 8:
            #CAMS format
            ref_time = datetime.datetime.strptime(timestr, "%Y%m%d")
        else: #len=10
            #Old MACC_V2014 format
            ref_time = datetime.datetime.strptime(timestr, "%Y%m%d%H")

        tpoints = [time_unit.date2num(
            ref_time+datetime.timedelta(hours=int(time)))
                   for time in tcoord.points]
        tcoord_new = iris.coords.DimCoord(tpoints, standard_name='time',
                                          units=time_unit)
        #Remove old time coordinate
        cube.remove_coord('time')
        #Replace with new time coordinate
        if self.forecast_day == 'forecast' or self.forecast_day == 'latest':

            #Need to cope with cubes having multiple times in
            #(unlike pp reading where each cube only has one field, ie one time)
            #Add extra coordinate 'itime' with unique points to ensure cubes can
            #be concatenated - This will be removed soon after returning from
            #callback.
            #The itime points will be a strictly monontic array of integers,
            #starting from zero.
            #itime will become the dim coord and
            #time left as an aux coord for now

            cube.add_aux_coord(tcoord_new, 0)

            #Generate an list of integer points, starting from the last value of
            #self.itime+1 and going up by 1 for each point.
            itime_pts = [self.itime + i for i in range(0, len(tpoints))]
            #Now increment self.itime so next time this routine is called,
            #the same selection of integers are not chosen.
            self.itime += len(tpoints)
            cube.add_dim_coord(iris.coords.DimCoord(itime_pts,
                                                    long_name='itime'),
                               0)

        else:
            #Add time on as the dim coord
            cube.add_dim_coord(tcoord_new, 0)

        #Also add new coordinate, forecast_period, which contains same
        #points as original time coordinate
        cube.add_aux_coord(iris.coords.DimCoord(tcoord.points,
                                                standard_name='forecast_period',
                                                units=cf_units.Unit('hours')),
                           0)

        #Also add forecast_day coordinate
        forecast_days = array_statistics.calc_forecast_day(tcoord.points,
                                                           ref_time.hour)
        cube.add_aux_coord(iris.coords.AuxCoord(forecast_days,
                                                long_name='forecast_day',
                                                units=cf_units.Unit('days')),
                           0)

        return cube

    def extract_sites(self, sites_data=None):
        """
        Extract site-specific data from gridded_cube_list into
        sites_cube_list, given

        :param sites_data: site information from :class:`sites_info.SitesInfo`
        """

        if self.sites_data is None:
            self.sites_data = sites_data
        if sites_data is None:
            if self.sites_data is None:
                raise ValueError("No sites requested")
            else:
                sites_data = self.sites_data

        self.sites_cube_list = self.extract_scl_from_gridded(sites_data,
                                                             'longitude',
                                                             'latitude')

        if self.forecast_day == 'all':
            #Sort out so that different forecast_days are put into different
            #cubes (otherwise day 1,3,5 are in one cube, 0,2,4 in another)
            self.sites_cube_list = cube_time.extract_all_forecast_days(
                self.sites_cube_list)

        return self.sites_cube_list

    def get_fileinfo(self, filename):
        """
        Get information from filename such as start and endtime and
        forecast_day, returns information in a dictionary.

        There are four formats of filename currently supported:

          * <model>_yyyymmddhh+hhh.nc, eg HRES_ENS_2014032600+012.nc
            These are distinguished by the '+' in the filename
            Each of these files are assumed to contain a single time.
            Each file contains multiple species.
          * <model>_<species>_hhH_hhH_yyyymmddHHMM000000.nc
            These are distinguished by not having '+' in the filename
            Each of these files contain 24 (or 25) hours, with the forecast
            leadtimes covered being given by hhH_hhH.
            Each file contains a single species.
          * <model(part1)>_<model(part2)>_<height>_yyyymmdd.nc
            These are similar to the last format, but all species are in the
            same file. These are as used for hindcast verification.
          * W_fr-meteofrance,MODEL,<modelname>+<datatype>+<levels>+<species>+
            <forecastperiod>_C_LFPW_yyyymmddHHMMSS.nc'.
            Each file contains a single species for a 24 or 25 hour period.
            Where <modelname> is the model name - ENSEMBLE or other CAMS member,
            eg CHIMERE; <datatype> is usually FORECAST (or possibly ANALYSIS
            - not tested); <levels> is SURFACE or ALLLEVELS, <species> is
            O3, CO, NO2, SO2, PM25, PM10, PANS, NMVOC, NO, NH3, or BIRCHPOLLEN;
            <forecastperiod> is 0H24H ,25H48H ,49H72H or 73H96.
            For full information about this, see
            http://www.regional.atmosphere.copernicus.eu/doc/\
Guide_Numerical_Data_CAMS_new.pdf



        Assumes forecast day 1 runs from 01Z-24Z.

        >>> maccens = MaccEnsData()
        >>> filename = '/path/to/data/HRES_ENS_2014032600+008.nc'
        >>> fileinfo = maccens.get_fileinfo(filename)
        >>> for key in sorted(fileinfo):
        ...     print(key, fileinfo[key])
        data_end 2014-03-26 08:00:00
        data_start 2014-03-26 08:00:00
        day1_date 2014-03-26
        forecast_days [1]
        forecast_periods [8]
        model HRES_ENS

        >>> filename = ('/path/to/data/ENSEMBLE_FORECAST_SURFACE_O3'
        ... '_25H_48H_201606140000000000.nc')
        >>> fileinfo = maccens.get_fileinfo(filename)
        >>> for key in sorted(fileinfo):
        ...     print(key, fileinfo[key])
        data_end 2016-06-16 00:00:00
        data_start 2016-06-15 01:00:00
        day1_date 2016-06-14
        forecast_days [2]
        forecast_periods [25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, \
37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48]
        model ENSEMBLE
        short_name O3

        >>> filename = ('/path/to/data/HRES_ENS_SFC_20150617.nc')
        >>> fileinfo = maccens.get_fileinfo(filename)
        >>> for key in sorted(fileinfo):
        ...     print(key, fileinfo[key])
        data_end 2015-06-18 00:00:00
        data_start 2015-06-17 00:00:00
        day1_date 2015-06-17
        forecast_days [0, 1]
        forecast_periods [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, \
13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]
        model HRES_ENS

        >>> filename = ('/path/to/data/W_fr-meteofrance,MODEL,ENSEMBLE+FORECAST'
        ... '+SURFACE+O3+0H24H_C_LFPW_20160914000000.nc')
        >>> fileinfo = maccens.get_fileinfo(filename)
        >>> for key in sorted(fileinfo):
        ...     print(key, fileinfo[key])
        data_end 2016-09-15 00:00:00
        data_start 2016-09-14 00:00:00
        day1_date 2016-09-14
        forecast_days [0, 1]
        forecast_periods [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, \
13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]
        model ENSEMBLE
        short_name O3


        """

        fileinfo = {}
        #Remove path if attached
        filename = os.path.basename(filename)
        split_file = filename.split('_')

        dt_string = split_file[-1][:-3]

        #In CAMS formatted files as opposed to older MACC data?
        cams_format = False

        if split_file[1].split(',')[0] == 'fr-meteofrance':
            #W_fr-meteofrance,MODEL,ENSEMBLE+FORECAST+SURFACE+O3+0H24H
            #_C_LFPW_20160914000000.nc
            #For more info, see
            #http://www.regional.atmosphere.copernicus.eu/doc/
            #Guide_Numerical_Data_CAMS_new.pdf
            cams_format = True
            ref_time_dt = datetime.datetime.strptime(dt_string, "%Y%m%d%H%M%S")
            comma_split = split_file[1].split(',')
            composite = comma_split[2]
            composite_split = composite.split('+')
            leadtime_start = int(composite_split[4].split('H')[0])
            leadtime_end = int(composite_split[4].split('H')[1])

            fileinfo['model'] = composite_split[0]
            fileinfo['short_name'] = composite_split[3]
            if fileinfo['short_name'] == 'PM25':
                fileinfo['short_name'] = 'PM2p5'

        elif '+' in dt_string:
            #Old format MACC_V2014
            #HRES_ENS_2014032600+008.nc'
            ref_time_dt = datetime.datetime.strptime(dt_string.split('+')[0],
                                                     "%Y%m%d%H")
            leadtime = int(dt_string.split('+')[1])
            leadtime_delta = datetime.timedelta(hours=leadtime)

            fileinfo['data_start'] = ref_time_dt + leadtime_delta
            #As only one hour in file, data_end is the same as data_start:
            fileinfo['data_end'] = fileinfo['data_start']
            fileinfo['forecast_periods'] = [leadtime]
            forecast_day = 1+(fileinfo['data_start']
                              -ref_time_dt.replace(hour=1)).days
            fileinfo['forecast_days'] = [forecast_day]
            fileinfo['model'] = '_'.join(split_file[:-1])


        elif len(dt_string) == 8:
            #CAMS hindcast format <model>_<levels>_yyyymmdd.nc
            cams_format = True
            ref_time_dt = datetime.datetime.strptime(dt_string, "%Y%m%d")
            #Assumption - all files are T+0 - T+24 (25 hours)
            leadtime_start = 0
            leadtime_end = 24
            fileinfo['model'] = '_'.join(split_file[:2])
        else:
            #CAMS routinely running format
            # <model>_<species>_hhH_hhH_yyyymmddHHMM000000.nc
            cams_format = True
            ref_time_dt = datetime.datetime.strptime(dt_string[:10],
                                                     "%Y%m%d%H")

            leadtime_start = int(split_file[-3][:-1])
            leadtime_end = int(split_file[-2][:-1])

            fileinfo['short_name'] = split_file[-4]
            if fileinfo['short_name'] == 'PM25':
                fileinfo['short_name'] = 'PM2p5'
            fileinfo['model'] = split_file[0]

        if cams_format:

            fileinfo['data_start'] = ref_time_dt + \
                                     datetime.timedelta(hours=leadtime_start)
            fileinfo['data_end'] = ref_time_dt + \
                                   datetime.timedelta(hours=leadtime_end)

            fileinfo['forecast_periods'] = list(range(leadtime_start,
                                                      leadtime_end+1))
            if leadtime_start == 0:
                #Most files only have a single forecast day in.
                #But the first file has T+0-T+24.
                #T+0 corresponds to day 0, while T+1-T+24 is day 1.
                #So also add day 0 into forecast_days list.
                fileinfo['forecast_days'] = [0, leadtime_end//24]
            else:
                fileinfo['forecast_days'] = [leadtime_end//24]

        #Get date of first full day of forecast for this model run.
        if ref_time_dt.hour == 0:
            #First full day is for today
            fileinfo['day1_date'] = ref_time_dt.date()
        else:
            fileinfo['day1_date'] = ref_time_dt.date() \
                                    + datetime.timedelta(days=1)

        return fileinfo

    def get_filenames(self, directory, start_datetime=None, end_datetime=None,
                      forecast_day=None, short_name_list=None, model=None):

        """
        Get file names from a directory, possibly limited by start and end dates
        (in datetime format) or forecast_day

        :param directory: string, directory containing files
        :param start_datetime: datetime object start date
        :param end_datetime: datetime object end date
        :param forecast_day: This can be an integer, whose number corresponds
                             to the forecast day, eg 1 if using the first full
                             forecast day from each model run.
                             Alternatively, set to "forecast" to get all
                             leadtimes from a single model run.
        :param short_name_list: List of short names to try and limit files by
        :param model: Which model to limit files by

        >>> import config
        >>> directory = config.SAMPLE_DATADIR+'macc_ens/'
        >>> maccens = MaccEnsData()
        >>> filenames = maccens.get_filenames(directory,
        ...     start_datetime=datetime.datetime(2014,3,26,21),
        ...     end_datetime=datetime.datetime(2014,3,27,3),
        ...     forecast_day=1)
        >>> print(filenames) # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
        ['.../macc_ens/HRES_ENS_2014032600+021.nc', \
        '.../macc_ens/HRES_ENS_2014032600+022.nc', \
        '.../macc_ens/HRES_ENS_2014032600+023.nc', \
        '.../macc_ens/HRES_ENS_2014032600+024.nc', \
        '.../macc_ens/HRES_ENS_2014032700+001.nc', \
        '.../macc_ens/HRES_ENS_2014032700+002.nc', \
        '.../macc_ens/HRES_ENS_2014032700+003.nc']

        An example of getting all files from a single forecast run:

        >>> maccens = MaccEnsData()
        >>> filenames = maccens.get_filenames(directory,
        ... start_datetime=datetime.datetime(2014, 3, 26),
        ... forecast_day='forecast')
        >>> print(filenames) # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
        ['.../macc_ens/HRES_ENS_2014032600+000.nc',
        '.../macc_ens/HRES_ENS_2014032600+001.nc',
        '.../macc_ens/HRES_ENS_2014032600+002.nc',
        ...
        '.../macc_ens/HRES_ENS_2014032600+094.nc',
        '.../macc_ens/HRES_ENS_2014032600+095.nc',
        '.../macc_ens/HRES_ENS_2014032600+096.nc']

        """

        if forecast_day is not None:
            self.forecast_day = forecast_day
        if start_datetime is not None:
            self.start_datetime = start_datetime
        if end_datetime is not None:
            self.end_datetime = end_datetime
        if short_name_list is not None:
            self.short_name_list = short_name_list
        if model is not None:
            self.model = model

        self.filenames = []
        if not directory.endswith('/'):
            directory += '/'
        listdir = os.listdir(directory)
        listdir.sort()

        #Set up variable to hold the date of the last day 1 in all the files
        #(note as the files are checked in reverse order, this
        # will be in the first file found)
        #This will be used for limiting files for forecast_day='latest'
        last_day1_date = None

        #Loop through files in reverse order,
        # to make finding latest files easier
        for filename in reversed(listdir):

            if not filename.endswith('.nc'):
                #Not a netcdf file - discard
                continue

            fileinfo = self.get_fileinfo(filename)

            #Discard file if doesn't contain required short-name
            if self.short_name_list is not None:
                if 'short_name' in fileinfo:
                    if fileinfo['short_name'] not in self.short_name_list:
                        continue

            #Discard file if not for required model
            if self.model is not None:
                if fileinfo['model'] != self.model:
                    continue

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

            self.filenames.append(directory+filename)
            #Now reverse back to get in sorted order again
            self.filenames.sort()

        return self.filenames

if __name__ == '__main__':

    import doctest
    doctest.testmod()
