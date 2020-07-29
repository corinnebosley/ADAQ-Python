"""
Class for reading from nimrod data files into an ADAQData class.
"""
from __future__ import print_function

from six.moves.builtins import str
import os
import datetime
import glob
import warnings

import numpy as np
import iris
import cf_units

import adaq_data
import cube_time
import shell_commands

#: Conversion between nimrod names and CF standard names
NIMROD_2_STDNAME = {
    'screen sulphur dioxide' : 'mass_concentration_of_sulphur_dioxide_in_air',
    'screen nitrogen dioxide' : 'mass_concentration_of_nitrogen_dioxide_in_air',
    'screen ozone' : 'mass_concentration_of_ozone_in_air',
    'screen PM10' : 'mass_concentration_of_pm10_ambient_aerosol_in_air',
    'screen PM2.5' : 'mass_concentration_of_pm2p5_ambient_aerosol_in_air',
    'screen carbon monoxide' : 'mass_concentration_of_carbon_monoxide_in_air',
    'screen nitric oxide' : 'mass_concentration_of_nitrogen_monoxide_in_air'
    }

#: Conversion between nimrod names and shortnames
NIMROD_2_SHORTNAME = {'screen sulphur dioxide'  : 'SO2',
                      'screen nitrogen dioxide' : 'NO2',
                      'screen ozone'            : 'O3',
                      'screen PM10'             : 'PM10',
                      'screen PM2.5'            : 'PM2p5',
                      'screen carbon monoxide'  : 'CO',
                      'screen nitric oxide'     : 'NO',
                      'daily air quality'       : 'DAQI'}



class NimrodData(adaq_data.ADAQData):
    """
    Subclass of ADAQData, which contains extra functionality specific to
    nimrod data.

    **Example:**

    >>> import config
    >>> from collections import OrderedDict
    >>> directory = config.SAMPLE_DATADIR+'aqum_output/nimrod/'

    Initialise class:

    >>> ND = NimrodData()

    Determine which filenames should be read in according to required
    start and end dates:

    >>> filenames = ND.get_filenames(directory,
    ...     start_datetime=datetime.datetime(2014, 4, 1, 12),
    ...     end_datetime=datetime.datetime(2014, 4, 2, 12),
    ...     forecast_day=1)
    >>> print(filenames[0]) # doctest: +ELLIPSIS
    /.../aqum_output/nimrod/201403311800_u1096_ng_aqum_airquality_2km

    Now read the data from nimrod files into a gridded_cube_list,
    by giving required filenames. Also limit to just
    raw AQUM data, day 1 forecast only and just a few species:

    >>> gcl = ND.readdata(filenames,data_type='raw',forecast_day=1,
    ... short_name_list = ['SO2','O3','PM2p5'])

    .. Note:: If filenames set using ND.get_filenames, then don't need to \
    pass filenames in as an argument to readdata.

    Check that the expected cubes were generated:

    >>> cube_name_list = ['mass_concentration_of_sulphur_dioxide_in_air',
    ... 'mass_concentration_of_ozone_in_air',
    ... 'mass_concentration_of_pm2p5_ambient_aerosol_in_air']
    >>> for cube in ND.gridded_cube_list:
    ...    assert(cube.name() in cube_name_list)

    Note: In Python 3, cubes may be generated in a different order to your
    input list, so comparison with an expected list could be ineffective.

    Have a look at the first cube only:

    >>> print(ND.gridded_cube_list[0])  # doctest: +NORMALIZE_WHITESPACE
    mass_concentration_of_sulphur_dioxide_in_air / (ug/m3) \
(time: 25; projection_y_coordinate: 704; projection_x_coordinate: 548)
         Dimension coordinates:
              time                                              x \
                           -                             -
              projection_y_coordinate                           - \
                           x                             -
              projection_x_coordinate                           - \
                           -                             x
         Auxiliary coordinates:
              forecast_reference_time                           x \
                           -                             -
         Scalar coordinates:
              experiment_number: 0
              forecast_day: 1
              height: 1.65 m
         Attributes:
              data_type: raw
              field_code: 771
              label: Nimrod
              nimrod_version: 2
              num_model_levels: 1
              short_name: SO2
              source: aqum

    Note this has added forecast_day as a coordinate,
    plus data_type, label and short_name as attributes.

    Now extract data at specific sites. Firstly, get some sites data:

    >>> import sites_info
    >>> sitesfilename = (config.SAMPLE_DATADIR +
    ...     'AURN_obs/aq_sites_GEMSRAQ_v4b_dev.txt')
    >>> obsdir = config.SAMPLE_DATADIR+'AURN_obs/'
    >>> sites = sites_info.SitesInfo()
    >>> sites_data = sites.read_from_file(sitesfilename,
    ... allsites=False, obsdir=obsdir)
    Number of sites:  5

    Now extract these sites from nimrod data, into a sites_cube_list:

    >>> scl = ND.extract_sites(sites_data=sites_data)

    Have a look at the SO2:

    >>> so2_sites_cube = ND.extract(short_name='SO2',singlecube=True)
    >>> print(so2_sites_cube)   # doctest: +NORMALIZE_WHITESPACE
    mass_concentration_of_sulphur_dioxide_in_air / (ug/m3) \
(site_id: 5; time: 25)
         Dimension coordinates:
              site_id                                              x        -
              time                                                 -        x
         Auxiliary coordinates:
              abbrev                                               x        -
              latitude                                             x        -
              longitude                                            x        -
              projection_x_coordinate                              x        -
              projection_y_coordinate                              x        -
              site_altitude                                        x        -
              site_name                                            x        -
              site_type                                            x        -
              forecast_reference_time                              -        x
         Scalar coordinates:
              experiment_number: 0
              forecast_day: 1
              height: 1.65 m
         Attributes:
              data_type: raw
              field_code: 771
              label: Nimrod
              nimrod_version: 2
              num_model_levels: 1
              short_name: SO2
              source: aqum

    Check sensible values:

    >>> gridded_cube = ND.extract(short_name='SO2',gridded=True,singlecube=True)
    >>> print("{:.2f}".format(gridded_cube.data.max()))
    149.90
    >>> so2_sites_cube.data.max()
    16.702557

    Alternatively, if want to read in data from a single forecast run, use the
    keyword setting forecast_day='forecast' when calling get_filenames and
    readdata. For more information about this setting, see
    :func:`cube_time.extract_latest_forecast_days`.

    Or to read in all the latest available data, use forecast_day='latest'.
    Again for information about this, see
    :func:`cube_time.extract_latest_forecast_days`.

    Finally, to read in files with all available forecast days and then force
    each different forecast day into a different cube, use the keyword
    forecast_day='all'. For more information about this setting see
    :func:`cube_time.extract_all_forecast_days`.

    """


    def __init__(self, label='Nimrod'):
        """
        Initiates a class from (new format) nimrod data
        as a subclass of ADAQData.
        """

        adaq_data.ADAQData.__init__(self)

        self.label = label #Label
        self.short_name_list = None #List of short names
        self.start_datetime = None #Starting datetime
        self.end_datetime = None #End datetime
        self.filenames = None #List of raw data filenames
        self.forecast_day = None #Required forecast days, either an integer
                                 #(typically 1), or string, eg 'forecast'
        self.sites_data = None #sites_data object
        self.data_type = None #'raw'/'sppo' or None (ie all available)

    def readdata(self, filenames=None, data_type=None, short_name_list=None,
                 forecast_day=None, start_datetime=None, end_datetime=None,
                 label=None):
        """
        Create gridded_cube_list from filenames specified.

        :param filenames: list of files to read, can include wildcards etc
        :param data_type: 'raw' or 'sppo' or None (gets everything available)
        :param short_name_list: list of short names (None gets everything)
        :param forecast_day: This can be an integer, whose number corresponds
                             to the forecast day, eg 1 if using the first full
                             forecast day from each model run.
                             Alternatively, set to "forecast" to get all
                             leadtimes from a single model run.
                             Or set to "latest" to get the latest available
                             data, which would be day 1 forecasts where
                             available and then the forecast from a single run
                             at the end.
        :param start_datetime: datetime object start date
        :param end_datetime: datetime object end date
        :param label: label for the cubes

        Note nimrod data in iris are given in `iris.fileformats.nimrod` and
        a description of the format can be found at
        http://fcm1.metoffice.com/projects/PostProc/wiki/PostProcDocFiles

        """

        if filenames is not None:
            self.filenames = filenames
        if data_type is not None:
            self.data_type = data_type
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

        if self.filenames is None:
            raise ValueError('No filenames given')

        constraints = None

        if not(self.forecast_day is None or self.forecast_day == 'forecast' or \
           self.forecast_day == 'latest' or self.forecast_day == 'all'):
            fcst_day_constraint = iris.Constraint(
                forecast_day=int(self.forecast_day))
            constraints = constraints & fcst_day_constraint

        if self.short_name_list is not None:
            sname_constraint = iris.AttributeConstraint(
               short_name=lambda c: c in self.short_name_list)
            constraints = constraints & sname_constraint

        if self.data_type is not None:
            datatype_constraint = iris.AttributeConstraint(data_type
                                                           =self.data_type)
            constraints = constraints & datatype_constraint

        if self.start_datetime is not None:
            time_constraint = iris.Constraint(
                time=lambda c: c.point >= self.start_datetime)
            constraints = constraints & time_constraint

        if self.end_datetime is not None:
            time_constraint = iris.Constraint(time=lambda c:
                                              c.point <= self.end_datetime)
            constraints = constraints & time_constraint


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

        return self.gridded_cube_list


    def __callback(self, cube, field, filename):
        """
        Private method to determine if data is from Raw or SPPO data,
        plus to add leadtime coordinate and forecast_day,
        correct units and scale data if required.
        Uses field data from `iris.fileformats.nimrod`.
        Nimrod format is described at
        http://fcm1.metoffice.com/projects/PostProc/wiki/PostProcDocFiles
        """

        #------------------------
        #Determine if Raw or SPPO Data
        if field.spare2 == 160:
            #raw data
            cube.attributes['data_type'] = 'raw'
        elif field.spare2 == 128:
            cube.attributes['data_type'] = 'sppo'

        #------------------------
        #Calculate forecast day.
        #A forecast day is considered to be 01Z-24Z.
        #Day 1 forecast is therefore the first full forecast day
        #So runtime of 18Z gives day 0 for 18-24Z,
        # day 1 for 01-24Z following day.

        meaningperiod = np.float(field.period_minutes)

        #For DAQI, 24hr meaningperiod is mis-leading for
        # calculating forecast day
        #meanperiod is normally 1 hour for other species
        if cube.name() == 'daily air quality' and meaningperiod == 60.*24:
            meaningperiod = 0.

        #Note there is an assumption that there are no 12Z runtimes
        #for nimrod data (we don't have any nimrod output from historical
        #12Z AQUM runs).
        if field.dt_hour == 12:
            raise ValueError('12Z runs not setup for calculating forecast_day'
                             'correctly')

        dt = datetime.datetime(field.dt_year, field.dt_month, field.dt_day,
                               field.dt_hour, field.dt_minute)
        vt = datetime.datetime(field.vt_year, field.vt_month, field.vt_day,
                               field.vt_hour, field.vt_minute)
        vd = (vt-datetime.timedelta(minutes=meaningperiod)).date() #valid-date
        dd = dt.date() #data-date
        forecast_day = (vd-dd).days

        if not cube.coords('forecast_day'):
            cube.add_aux_coord(iris.coords.DimCoord(forecast_day,
                                                    long_name='forecast_day'))

        #------------------------
        #Correct data scaling
        #(see nimrod file format description for more detail)
        cube.data = cube.data * field.MKS_data_scaling + field.data_offset

        #------------------------
        #Correct units
        if str(field.units.strip()) == 'ug/m3E1':
            if field.MKS_data_scaling == np.float32(0.1):
                cube.units = 'ug/m3'
                #Remove invalid units
                del cube.attributes['invalid_units']

        #------------------------
        #Set standard and short names if known
        cubename = cube.name()
        if cubename in NIMROD_2_STDNAME:
            cube.rename(NIMROD_2_STDNAME[cubename])
        if cubename in NIMROD_2_SHORTNAME:
            cube.attributes['short_name'] = NIMROD_2_SHORTNAME[cubename]
        else:
            cube.attributes['short_name'] = cubename

        #Give label to cube
        cube.attributes['label'] = self.label

        # ------------------------------------------------------
        # Change standard calendar to a gregorian calendar to
        # avoid problems with newer version of iris.
        # Get the time coordinate with a standard calendar
        tcoord_standard = cube.coord('time')

        # Create a new cf_unit object with the gregorian calendar
        greg_unit = cf_units.Unit(tcoord_standard.units.name,
                                  calendar='gregorian')

        # Make a new time coordinate with the same points and bounds
        # but using a gregorian calendar
        greg_tcoord = iris.coords.DimCoord(
            tcoord_standard.points,
            bounds=tcoord_standard.bounds,
            standard_name=tcoord_standard.standard_name,
            units=greg_unit)

        # delete old time coordinate and add the new one
        cube.remove_coord('time')
        cube.add_aux_coord(greg_tcoord)


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

        self.sites_cube_list = self.extract_scl_from_gridded(
            sites_data,
            'projection_x_coordinate',
            'projection_y_coordinate')

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

        >>> ND = NimrodData()
        >>> filename = '/path_to_data/201403271800_u1096_ng_aqum_airquality_2km'
        >>> fileinfo = ND.get_fileinfo(filename)
        >>> fileinfo == {'forecast_times':
        ... {0: (datetime.datetime(2014, 3, 27, 19, 0),
        ...      datetime.datetime(2014, 3, 28, 0, 0)),
        ...  1: (datetime.datetime(2014, 3, 28, 1, 0),
        ...      datetime.datetime(2014, 3, 29, 0, 0)),
        ...  2: (datetime.datetime(2014, 3, 29, 1, 0),
        ...      datetime.datetime(2014, 3, 30, 0, 0)),
        ...  3: (datetime.datetime(2014, 3, 30, 1, 0),
        ...      datetime.datetime(2014, 3, 31, 0, 0)),
        ...  4: (datetime.datetime(2014, 3, 31, 1, 0),
        ...      datetime.datetime(2014, 4, 1, 0, 0)),
        ...  5: (datetime.datetime(2014, 4, 1, 1, 0),
        ...      datetime.datetime(2014, 4, 2, 0, 0))},
        ... 'forecast_days': [0, 1, 2, 3, 4, 5],
        ... 'day1_date': datetime.date(2014, 3, 28),
        ... 'data_end': datetime.datetime(2014, 4, 2, 0, 0),
        ... 'data_start': datetime.datetime(2014, 3, 27, 18, 0)}
        True

        An example from a potential 00Z nimrod file:

        >>> filename = '/path_to_data/201403280000_u1096_ng_aqum_airquality_2km'
        >>> fileinfo = ND.get_fileinfo(filename)
        >>> fileinfo == {'forecast_times':
        ... {1: (datetime.datetime(2014, 3, 28, 1, 0),
        ...      datetime.datetime(2014, 3, 29, 0, 0)),
        ...  2: (datetime.datetime(2014, 3, 29, 1, 0),
        ...      datetime.datetime(2014, 3, 30, 0, 0)),
        ...  3: (datetime.datetime(2014, 3, 30, 1, 0),
        ...      datetime.datetime(2014, 3, 31, 0, 0)),
        ...  4: (datetime.datetime(2014, 3, 31, 1, 0),
        ...      datetime.datetime(2014, 4, 1, 0, 0)),
        ...  5: (datetime.datetime(2014, 4, 1, 1, 0),
        ...      datetime.datetime(2014, 4, 2, 0, 0))},
        ... 'forecast_days': [1, 2, 3, 4, 5],
        ... 'day1_date': datetime.date(2014, 3, 28),
        ... 'data_end': datetime.datetime(2014, 4, 2, 0, 0),
        ... 'data_start': datetime.datetime(2014, 3, 28, 0, 0)}
        True

        .. note:: If the filename contains 'airquality' then there is an
                  assumption that the file contains 5 days of forecasts.

        """

        fileinfo = {}
        #Remove path if attached
        filename = os.path.basename(filename)
        split_file = filename.split('_')
        #Get start time - first section of filename
        try:
            fileinfo['data_start'] = datetime.datetime.strptime(split_file[0],
                                                                "%Y%m%d%H%M")
        except:
            #Failed to convert date
            #- must be invalid file as not named correctly.
            return fileinfo

        dt = fileinfo['data_start'] #data time
        if dt.hour == 0:
            day1_date = dt.date()
        else:
            day1_date = dt.date()+datetime.timedelta(days=1)

        fileinfo['day1_date'] = day1_date
        fileinfo['data_end'] = None
        fileinfo['forecast_days'] = None
        fileinfo['forecast_times'] = {}

        if 'airquality' in split_file:
            #Can now make various assumptions:
            #In particular that it is a 5 day forecast
            n_fc_days = 5

            #Get a list of days which this file corresponds to.
            #So 1 would be the first full day of the forecast

            fileinfo['forecast_days'] = list(np.arange(n_fc_days)+1)
            if dt.hour > 0:
                #If dt.hour > 0, then there is also a day0 forecast component
                fileinfo['forecast_days'].insert(0, 0)

            #Also calculate the end datetime
            if dt.hour <= 12:
                forecast_hrs = n_fc_days*24 + dt.hour
                fileinfo['data_end'] = fileinfo['data_start'] \
                                        +datetime.timedelta(hours=forecast_hrs)
            else:
                forecast_hrs = n_fc_days*24 + (24-dt.hour)
                fileinfo['data_end'] = fileinfo['data_start'] \
                                        +datetime.timedelta(hours=forecast_hrs)

        if fileinfo['forecast_days'] is not None:
            for day in fileinfo['forecast_days']:
                #Date/Times for this forecast day run from 1Z-24Z
                #Set up a tuple giving the first and last valid datetime for
                #data on this day
                fc_date = day1_date + datetime.timedelta(days=int(day)-1)
                day_start = datetime.datetime(fc_date.year, fc_date.month,
                                              fc_date.day, 1)
                day_end = day_start + datetime.timedelta(hours=23) #0Z
                if day == 0:
                    #First valid data output one hour after model run starts
                    day_start = fileinfo['data_start'] + \
                                datetime.timedelta(hours=1)
                fileinfo['forecast_times'][day] = (day_start, day_end)

        return fileinfo


    def get_filenames(self, directory=None, start_datetime=None,
                      end_datetime=None, forecast_day=None, recursive=True,
                      directory_list=None):
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
        :param recursive: Recursively look through all subdirectories for files.
                          Note this is often required for nimrod as different
                          days are located in different directories.
        :param directory_list: list of files from directory


        Usage:

        >>> import config
        >>> directory = config.SAMPLE_DATADIR+'aqum_output/nimrod'
        >>> ND = NimrodData()
        >>> filenames = ND.get_filenames(directory,
        ...     start_datetime=datetime.datetime(2014, 4, 1),
        ...     end_datetime=datetime.datetime(2014, 4, 3),
        ...     forecast_day=1)
        >>> for filename in filenames: print(os.path.basename(filename))
        201403301800_u1096_ng_aqum_airquality_2km
        201403311800_u1096_ng_aqum_airquality_2km
        201404011800_u1096_ng_aqum_airquality_2km

        For getting filenames for a single forecast run:

        >>> filenames = ND.get_filenames(directory,
        ...     start_datetime=datetime.datetime(2014, 4, 1),
        ...     end_datetime=datetime.datetime(2014, 4, 3),
        ...     forecast_day='forecast')
        >>> print(filenames) # doctest: +ELLIPSIS
        ['.../201403311800_u1096_ng_aqum_airquality_2km']

        Or for the latest available data:

        >>> filenames = ND.get_filenames(directory,
        ...     start_datetime=datetime.datetime(2014, 4, 3),
        ...     end_datetime=datetime.datetime(2014, 4, 8),
        ...     forecast_day='latest')
        >>> print(filenames) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        ['.../201404021800_u1096_ng_aqum_airquality_2km',
        '.../201404031800_u1096_ng_aqum_airquality_2km',
        '.../201404041800_u1096_ng_aqum_airquality_2km',
        '.../201404051800_u1096_ng_aqum_airquality_2km']


        """

        if forecast_day is not None:
            self.forecast_day = forecast_day
        if start_datetime is not None:
            self.start_datetime = start_datetime
        if end_datetime is not None:
            self.end_datetime = end_datetime

        self.filenames = []
        if directory_list is None:
            if directory is None:
                raise IOError('Need either directory_list or directory '
                              'parameter set')
            if not directory.endswith('/'):
                directory += '/'
            if recursive:
                #Recursively search through all directories for files
                listdir = [os.path.join(dp, f)
                           for dp, dn, fn in os.walk(directory) for f in fn]
            else:
                #Only consider this directory
                listdir = [os.path.join(directory, filename)
                           for filename in os.listdir(directory)]
        else:
            listdir = directory_list
        listdir.sort()

        for filename in listdir:

            fileinfo = self.get_fileinfo(filename)

            if not fileinfo:
                #Not a valid file
                continue

            #Discard file if all data is before required start time
            if start_datetime is not None:
                if fileinfo['data_end'] is not None:
                    if fileinfo['data_end'] < start_datetime:
                        continue

            #Discard file if all data is after required end time
            if end_datetime is not None:
                if fileinfo['data_start'] is not None:
                    if fileinfo['data_start'] > end_datetime:
                        continue

            if self.forecast_day is not None:

                #Discard file if file does not contain required forecast_day
                if fileinfo['forecast_days'] is not None:
                    if self.forecast_day == 'forecast':
                        if fileinfo['day1_date'] != self.start_datetime.date():
                            continue
                    elif self.forecast_day == 'latest':
                        #Day 1-5 are in all files.
                        #Only need to get files where day 1 >= startdate
                        #as will at least need one file with day 1 in.
                        #Can not constrain further eg on end date
                        if fileinfo['day1_date'] < self.start_datetime.date():
                            continue
                    elif self.forecast_day == 'all':
                        #Keep all files
                        pass
                    else:
                        if (int(self.forecast_day) not in
                                fileinfo['forecast_days']):
                            continue

                if self.forecast_day != 'forecast' and \
                   self.forecast_day != 'latest' and \
                   self.forecast_day != 'all':
                    #Determine if required forecast day is within required
                    #start/end dates - discard if not
                    #Note this may give an extra file than required if the
                    #required end time is at 00Z.
                    if fileinfo['forecast_times']:
                        datematch = False
                        fd = int(self.forecast_day)

                        if fd in fileinfo['forecast_times']:
                            #Check that start time is before the end of the file
                            #and end time is after the start of the file
                            if start_datetime <= fileinfo['forecast_times'][fd][1] \
                               and end_datetime >= fileinfo['forecast_times'][fd][0]:
                                datematch = True
                        if not datematch:
                            continue

            self.filenames.append(filename)

        return self.filenames

def __mass_retrieve_filenames(forecast_day, start_datetime, end_datetime,
                              runtime):
    """
    Generate list of filenames that need to be retrieved from mass

    :param start_datetime: datetime format of start time from files
    :param end_datetime: datetime format of end time from files
    :param forecast_day: This can be an integer, whose number corresponds
                         to the forecast day, eg 1 if using the first full
                         forecast day from each model run.
                         Alternatively, set to "forecast" to get all leadtimes
                         from a single model run.
                         Or set to "latest" to get the latest available
                         data, which would be day 1 forecasts where
                         available and then the full forecast from a single
                         model forecast run at the end.
    :param runtime: runtime of model. Usually 18 (default) or 0.

    :returns: Dictionary, whose keys are years that the filename matches and
              whose values are a list of filenames for this year.

    Example:

    >>> filenames_dict = __mass_retrieve_filenames(
    ... forecast_day=1,
    ... start_datetime=datetime.datetime(2015,12,31),
    ... end_datetime=datetime.datetime(2016,1,3),
    ... runtime=18)

    >>> for key in sorted(filenames_dict.keys()):
    ...     print(key, filenames_dict[key])
    2015 ['201512291800_u1096_ng_aqum_airquality_2km.Z', \
'201512301800_u1096_ng_aqum_airquality_2km.Z', \
'201512311800_u1096_ng_aqum_airquality_2km.Z']
    2016 ['201601011800_u1096_ng_aqum_airquality_2km.Z']

    """

    filename_suffix = '00_u1096_ng_aqum_airquality_2km'

    #----
    # Get list of required filenames.
    # Assumption: required dates are within startdate-5days - enddate+5days

    assert start_datetime.date() <= end_datetime.date()

    #Generate a list of possible filenames, by manually creating fileanames
    #within startdate-5days - enddate+5days
    possible_filenames = []
    date = start_datetime.date() - datetime.timedelta(days=5)
    while date <= end_datetime.date() + datetime.timedelta(days=5):
        #zfill - ensure 2 characters in runtime
        possible_filenames.append(date.strftime('%Y%m%d')
                                  + str(runtime).zfill(2)
                                  + filename_suffix)
        date += datetime.timedelta(days=1)


    #Now setup nimrod data object so can call get_filenames and therefore
    #extract only the actually required filenames
    nd = NimrodData()
    nd.get_filenames(start_datetime=start_datetime,
                     end_datetime=end_datetime,
                     forecast_day=forecast_day,
                     directory_list=possible_filenames)

    #Split up this list of filenames to put in a dictionary
    #based on keys of year
    #Also add suffix of '.Z' onto the end of each filename
    #as the files are gzipped within mass
    filenames_dict = {}
    for filename in nd.filenames:
        year = filename[:4]
        if year in filenames_dict:
            filenames_dict[year].append(filename + '.Z')
        else:
            filenames_dict[year] = [filename + '.Z']

    return filenames_dict

def nimrod_mass_retrieve(outputdir,
                         start_datetime=None, end_datetime=None,
                         forecast_day=None, runtime=18,
                         massretries=0, massretrydelay=60,
                         retrieve=True):
    """
    Function which will set up and perform mass retrieval for operational
    nimrod files.

    :param outputdir: String. Directory to retrieve files into.
                      Will be created if does not already exist.
    :param start_datetime: datetime format of start time from files
    :param end_datetime: datetime format of end time from files
    :param forecast_day: This can be an integer, whose number corresponds
                         to the forecast day, eg 1 if using the first full
                         forecast day from each model run.
                         Alternatively, set to "forecast" to get all leadtimes
                         from a single model run.
                         Or set to "latest" to get the latest available
                         data, which would be day 1 forecasts where
                         available and then the full forecast from a single
                         model forecast run at the end.
    :param runtime: runtime of model. Usually 18 (default) or 0.
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

    This code works by determining the required filenames and then setting up
    a mass query file with these filenames in which is used by the moo select
    command. Following this the files are untarred and then unzipped.
    A separate mass retrieval is done for each required year.

    Example of calling, although here retrieve is set to False.

    >>> mooselect_strs, mass_cmds = nimrod_mass_retrieve(
    ... '$DATADIR/mass_retrieve/nimrod',
    ... start_datetime=datetime.datetime(2015,12,31),
    ... end_datetime=datetime.datetime(2016,1,3),
    ... forecast_day=1,
    ... runtime=18,
    ... massretries=0,
    ... massretrydelay=60,
    ... retrieve=False) # doctest: +ELLIPSIS
    Moose Retrieval Command:
    moo select -f .../mooselect_2015.txt moose:/opfc/atm/postpro/prods/2015.tar\
 moose:/opfc/atm/postpro/prods/2016.tar .../mass_retrieve/nimrod
    Moose Retrieval Command:
    moo select -f .../mooselect_2016.txt moose:/opfc/atm/postpro/prods/2016.tar\
 moose:/opfc/atm/postpro/prods/2017.tar .../mass_retrieve/nimrod

    Contents of the query files:

    >>> for mooselect_str in mooselect_strs:
    ...     print(mooselect_str)
    begin
      filename=("201512291800_u1096_ng_aqum_airquality_2km.Z",\
"201512301800_u1096_ng_aqum_airquality_2km.Z",\
"201512311800_u1096_ng_aqum_airquality_2km.Z")
    end
    begin
      filename="201601011800_u1096_ng_aqum_airquality_2km.Z"
    end

    Mass retrieval commands:

    >>> for cmd in mass_cmds:
    ...     print(cmd) # doctest: +ELLIPSIS
    moo select -f .../mooselect_2015.txt moose:/opfc/atm/postpro/prods/2015.tar\
 moose:/opfc/atm/postpro/prods/2016.tar .../nimrod
    moo select -f .../mooselect_2016.txt moose:/opfc/atm/postpro/prods/2016.tar\
 moose:/opfc/atm/postpro/prods/2017.tar .../nimrod

    """

    #Make output directory if does not already exist
    outputdir = os.path.expandvars(outputdir)
    if not os.path.isdir(outputdir):
        os.makedirs(outputdir)

    #Delete any tar or .Z files left from previous attempts
    returncode = shell_commands.call_shell('rm -f ' + outputdir + '*.tar')
    if returncode != 0:
        warnings.warn('Unable to delete any .tar files')
    returncode = shell_commands.call_shell('rm -f ' + outputdir + '*.Z')
    if returncode != 0:
        warnings.warn('Unable to delete any .Z files')

    filenames_dict = __mass_retrieve_filenames(
        forecast_day=forecast_day,
        start_datetime=start_datetime,
        end_datetime=end_datetime,
        runtime=runtime)

    #Store mass commands and strings used
    mass_cmds = []
    mooselect_strs = []

    #Loop over years to retrieve different years separately to
    #limit number of files retrieved at once.
    for year, filenames in filenames_dict.items():

        # Set up moo select filter file
        if len(filenames) == 1:
            mooselect_str = 'begin\n  filename="'
            mooselect_str += filenames[0]
            mooselect_str += '"\nend'
        else:
            mooselect_str = 'begin\n  filename=("'
            mooselect_str += '","'.join(filenames)
            mooselect_str += '")\nend'

        mooselect_strs.append(mooselect_str)

        selectfile = outputdir + '/mooselect_' + year + '.txt'

        with open(selectfile, 'w') as fout:
            fout.write(mooselect_str)

        #Setup archive to retrieve from on mass.
        #Retrieve files from the archive corresponding to the files' year and
        #also from the following year's archive if it is available.
        #This is because some files are saved in a UKV archive file with date
        #a few days later, which for files close to the year end may be the
        #following year.
        #To avoid mass failures, check if next year's archive is available by
        #comparing to current year. The archive will not be available until
        #the new year starts but it is assumed that it will be available by
        #the time this code is executed at the start of a new year.
        nextyear = str(int(year)+1)
        if int(nextyear) > datetime.datetime.now().year:
            massdir = 'moose:/opfc/atm/postpro/prods/' + year +'.tar '
        else:
            #Should exist
            massdir = ('moose:/opfc/atm/postpro/prods/' + year +'.tar '
                       'moose:/opfc/atm/postpro/prods/' + nextyear +'.tar')

        # Do retrieval
        mass_cmd = shell_commands.call_mass(selectfile, massdir, outputdir,
                                            massretries=massretries,
                                            massretrydelay=massretrydelay,
                                            retrieve=retrieve)
        mass_cmds.append(mass_cmd)

        if retrieve:

            #Un-tar the files and then remove .tar file
            print('Untarring files...')
            for filename in glob.glob(outputdir + '/*.tar'):
                print(filename)
                returncode = shell_commands.call_shell(
                    'cd '+outputdir+'; tar -xvf ' + filename)
                if returncode != 0:
                    warnings.warn('Unable to un-tar file: ' + filename)
                else:
                    print('Deleting ' + filename)
                    returncode = shell_commands.call_shell('rm ' + filename)
                    if returncode != 0:
                        warnings.warn('Unable to delete file: ' + filename)

            #Unzip the files and then remove .Z file
            print('Gunzipping files...')
            for filename in glob.glob(outputdir + '/*.Z'):
                print(filename)
                #First remove any previously unzipped versions of this file
                if os.path.exists(filename[:-2]):
                    print('Deleting pre-existing ' + filename[:-2])
                    returncode = shell_commands.call_shell(
                        'rm ' + filename[:-2])
                    if returncode != 0:
                        warnings.warn('Unable to delete file: ' + filename[:-2])
                returncode = shell_commands.call_shell('gunzip ' + filename)
                if returncode != 0:
                    warnings.warn('Unable to gunzip file: ' + filename)

    return mooselect_strs, mass_cmds



if __name__ == '__main__':

    import doctest
    doctest.testmod()
