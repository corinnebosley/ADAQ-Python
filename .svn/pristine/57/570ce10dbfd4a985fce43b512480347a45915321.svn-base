"""
Class for reading from CAMS air quality observation files into an
ADAQData class.
"""
from __future__ import print_function

from six.moves.builtins import str
import datetime
import os
import glob
import warnings

import cf_units
import iris
import iris.coord_categorisation as coord_cat
import numpy as np

import config
import adaq_functions
import sites_info
import adaq_data

#: Conversion between CAMS parameters and CF standard names
CAMSAQOBS_2_STDNAME = {
    'o3'   : 'mass_concentration_of_ozone_in_air',
    'no2'  : 'mass_concentration_of_nitrogen_dioxide_in_air',
    'co'   : 'mass_concentration_of_carbon_monoxide_in_air',
    'pm10' : 'mass_concentration_of_pm10_ambient_aerosol_in_air',
    'pm2p5': 'mass_concentration_of_pm2p5_ambient_aerosol_in_air',
    'so2'  : 'mass_concentration_of_sulfur_dioxide_in_air'}

#: Conversion between CAMS parameters and short_name
CAMSAQOBS_2_SHORTNAME = {'o3'    : 'O3',
                         'no2'   : 'NO2',
                         'co'    : 'CO',
                         'pm10'  : 'PM10',
                         'pm2p5' : 'PM2p5',
                         'so2'   : 'SO2'}


class CAMSAQObsData(adaq_data.ADAQData):
    """
    Subclass of ADAQData, which contains extra functionality specific to
    Air Quality Observation data from CAMS.

    This dataset is given in ascii text files, with a single file
    per day, which includes all species. The expected headers are:

    .. code-block:: none

        STATION;LAT;LON;ALT(m);PARAMETER;YEAR;MONTH;DAY;HOUR;AVERAGING_PERIOD(h);CONCENTRATION(kg/m3)

    The station code is converted to both an abbrev and the station name.
    Units are all converted to ug/m3.
    Currently assumes all sites are of type 'RURAL'.

    **Example:**

    Firstly set up required filenames (this could be a directory, or a list of
    files, or even use wildcards), plus start and end date-times:

    >>> filenames = config.SAMPLE_DATADIR+'CAMSAQ_obs'
    >>> start_dt = datetime.datetime(2015,4,4,0)
    >>> end_dt = datetime.datetime(2015,4,5,0)

    Set up observation data class and read data in:

    >>> od = CAMSAQObsData()
    >>> scl = od.readdata(filenames, start_datetime=start_dt,
    ... end_datetime=end_dt) # doctest: +ELLIPSIS
    Reading .../obsmacc4ana_20150403.csv
    Reading .../obsmacc4ana_20150404.csv
    Converting to cubes

    Sort the data (od.sites_cube_list) alphabetically before printing
    (only required for doctest):

    >>> scl = iris.cube.CubeList(sorted(
    ... od.sites_cube_list, key=lambda cube: cube.name()))
    >>> print(scl)
    0: mass_concentration_of_carbon_monoxide_in_air / (ug/m3) \
(site_id: 803; time: 25)
    1: mass_concentration_of_nitrogen_dioxide_in_air / (ug/m3) \
(site_id: 803; time: 25)
    2: mass_concentration_of_ozone_in_air / (ug/m3) \
(site_id: 803; time: 25)
    3: mass_concentration_of_pm10_ambient_aerosol_in_air / (ug/m3) \
(site_id: 803; time: 25)
    4: mass_concentration_of_pm2p5_ambient_aerosol_in_air / (ug/m3) \
(site_id: 803; time: 25)
    5: mass_concentration_of_sulfur_dioxide_in_air / (ug/m3) \
(site_id: 803; time: 25)

    The required species could also be used to limit the list returned,
    by using the short_name_list keyword.
    Similarly, the sites can be limited, by using first few letters of
    station names, for example to get Great Britain (GB) and France (FR),
    use abbrev_prefix_list = ['GB','FR']):

    >>> scl = od.readdata(filenames, start_datetime=start_dt,
    ... end_datetime=end_dt, short_name_list=['O3'],
    ... abbrev_prefix_list = ['GB','FR']) # doctest: +ELLIPSIS
    Reading .../obsmacc4ana_20150403.csv
    Reading .../obsmacc4ana_20150404.csv
    Converting to cubes

    >>> print(od.sites_cube_list[0])
    mass_concentration_of_ozone_in_air / (ug/m3) (site_id: 123; time: 25)
         Dimension coordinates:
              site_id                                    x          -
              time                                       -          x
         Auxiliary coordinates:
              abbrev                                     x          -
              latitude                                   x          -
              longitude                                  x          -
              site_altitude                              x          -
              site_classification_CO                     x          -
              site_classification_NO                     x          -
              site_classification_NO2                    x          -
              site_classification_NOx                    x          -
              site_classification_O3                     x          -
              site_classification_PM10                   x          -
              site_classification_PM2p5                  x          -
              site_classification_SO2                    x          -
              site_name                                  x          -
              site_type                                  x          -
         Attributes:
              label: Obs
              short_name: O3
              source: CAMSAQObs
         Cell methods:
              mean: time (1 hour)

    >>> print('{:.2f}'.format(np.nanmax(od.sites_cube_list[0].data)))
    102.00
    >>> print('{:.2f}'.format(np.nanmin(od.sites_cube_list[0].data)))
    1.00
    >>> print(od.sites_cube_list[0].coord('site_name')[:10])
    ... # doctest: +NORMALIZE_WHITESPACE  +ELLIPSIS
    AuxCoord(array(['FR05074', 'FR21031', 'FR12020', 'GB0038R', 'FR23182',
           'FR05040', 'GB0045R', 'FR05010', 'GB0617A', 'GB0728A'],
          dtype='...7'), standard_name=None, units=Unit('1'),
          long_name='site_name')

    To extract SitesInfo data from this, for example for use when extracting
    model data at these sites:

    >>> si = sites_info.SitesInfo()
    >>> sites_data = si.read_from_sites_cube(od.sites_cube_list[0])
    >>> fmtstr = sites_info.format_string(sites_data[:5])
    >>> for site in sites_data[:5]:
    ...    print(fmtstr.format(*site))
       10070.13948960,FR05074,  49.490, 0.100703,    5.000, 0, 3, 3, 3, 5, \
4, 0, 4,FR05074,   URBAN_BACKGROUND
       11083.13844560,FR21031,  48.446, 0.110833,  140.000, 0, 5, 5, 5, 4, \
8, 0, 0,FR21031,SUBURBAN_BACKGROUND
       17972.13363030,FR12020,  43.630, 0.179722,  240.000, 0, 1, 1, 1, 2, \
3, 0, 0,FR12020,   RURAL_BACKGROUND
       18125.14079370,GB0038R,  50.794, 0.181250,  125.000, 0, 1, 1, 1, 3, \
0, 0, 2,GB0038R,   RURAL_BACKGROUND
       22142.13797110,FR23182,  47.971, 0.221422,   85.000, 0, 4, 4, 4, 5, \
6, 0, 0,FR23182,   URBAN_BACKGROUND

    """

    def __init__(self, label='Obs'):
        """
        Initiates class as a subset of adaq_data.ADAQData,
        plus other observation-specific data.
        """

        adaq_data.ADAQData.__init__(self)

        self.label = label #Label
        self.short_name_list = None #List of short names
        self.start_datetime = None #Starting datetime
        self.end_datetime = None #End datetime
        self.filenames = None #List of raw data filenames
        self.abbrev_prefix_list = None
        self.sites_data = None
        self.site_types = None
        self.classes_file = None # File containing classification data
        self.classes_data = None # Data from classifications file
        self.abbrev2latlon = {} #Dictionary containing conversion from
                                #abbrev to (lat,lon)
        self.site_types = None

    def readdata(self, filenames=None, short_name_list=None,
                 sites_data=None,
                 start_datetime=None, end_datetime=None,
                 abbrev_prefix_list=None,
                 classes_file=None, site_types=None):
        """
        Create sites_cube_list from filenames specified.

        :param filenames: list of files to read, can include wildcards etc.
                          If this is a directory, then reads all files in
                          directory that end with .csv.
                          Note files are expected to be of the format
                          ``*yyyymmdd.csv``
        :param short_name_list: list of short names
                          (Setting to None [default] returns all available
                          parameters).
        :param start_datetime: Start time of data in datetime format.
        :param end_datetime: End time of data in datetime format.
        :param abbrev_prefix_list: List of strings. Only return sites whose
                                   abbrev starts with an element of this list.
        :param classes_file: Filename of file containing information about
                             different site types and classes. Uses file
                             from sample data directory if not given.
        :param site_types: List of strings. Only return sites whose
                           site type matches an element of this list.

        :returns: sites_cube_list

        """

        #Set class attributes on basis of keywords if set
        if filenames is not None:
            self.filenames = filenames
        #Raise an error if no filenames are set
        #(may have been set not as part of readdata keywords)
        if self.filenames is None:
            raise ValueError('No filenames given')

        if short_name_list is not None:
            self.short_name_list = short_name_list
        if start_datetime is not None:
            self.start_datetime = start_datetime
        if end_datetime is not None:
            self.end_datetime = end_datetime
        if abbrev_prefix_list is not None:
            self.abbrev_prefix_list = abbrev_prefix_list
        if sites_data is not None:
            self.sites_data = sites_data
        if classes_file is not None:
            self.classes_file = classes_file
        if site_types is not None:
            self.site_types = site_types

        #Read in classifications file
        self.__read_sites_classes()

        #Convert to filenames rather than directory if required
        if os.path.isdir(self.filenames):
            self.filenames += '/*'
        self.filenames = sorted(glob.glob(self.filenames))

        #Read all data in
        data = None
        #Set up fixed dtype to ensure that all files have same
        #types of data - particularly for PARAMETER which can get
        #reduced to S3 if PM10/PM2p5 not present in file.
        file_dtype = [('STATION', str, 7), ('LAT', '<f8'), ('LON', '<f8'),
                      ('ALTm', '<i8'), ('PARAMETER', str, 5), ('YEAR', '<i8'),
                      ('MONTH', '<i8'), ('DAY', '<i8'), ('HOUR', '<i8'),
                      ('AVERAGING_PERIODh', '<i8'),
                      ('CONCENTRATIONkgm3', '<f8')]
        for filename in self.filenames:

            if filename.endswith('.csv'):

                #Also check datestamp on file
                #(makes reading much quicker if don't need to read all files!
                file_dt_stamp = datetime.datetime.strptime(
                    filename.split('_')[-1], "%Y%m%d.csv")
                #File runs from 01Z on this date, to 24Z on this date
                file_dt_start = file_dt_stamp + datetime.timedelta(hours=1)
                file_dt_end = file_dt_stamp + datetime.timedelta(hours=24)
                #Compare to required start/end datetimes
                #- ignore this file if not in required range
                if self.start_datetime is not None:
                    if file_dt_end < self.start_datetime:
                        continue
                if self.end_datetime is not None:
                    if file_dt_start > self.end_datetime:
                        continue

                print("Reading", filename)

                #Read file in

                filedata = np.genfromtxt(filename, delimiter=';',
                                         dtype=file_dtype, names=True)

                if filedata.size:
                    #Only add to data if data is available
                    #note - file might be available, but only
                    #contains headers and no data.
                    if data is None:
                        data = filedata
                    else:
                        data = np.append(data, filedata)

        #Convert to cubelist
        if data is None:
            warnings.warn("No data read in")
            return self.sites_cube_list

        print("Converting to cubes")
        cubelist = self.__data_to_cubelist(data)

        #Apply iris constraints to limit returned data
        constraints = iris.Constraint()
        if self.start_datetime is not None:
            time_constraint = iris.Constraint(
                time=lambda c: c.point >= self.start_datetime)
            constraints = constraints & time_constraint

        if self.end_datetime is not None:
            time_constraint = iris.Constraint(
                time=lambda c: c.point <= self.end_datetime)
            constraints = constraints & time_constraint

        if self.short_name_list is None:
            sname_constraint = iris.AttributeConstraint()
        else:
            sname_constraint = iris.AttributeConstraint(
                short_name=lambda c: c in self.short_name_list)

        if self.abbrev_prefix_list is not None:
            abbrev_constraint = iris.Constraint(
                abbrev=lambda c: str(c).startswith(
                    tuple(self.abbrev_prefix_list)))
            constraints = constraints & abbrev_constraint

        if self.sites_data is not None:
            site_abbrevs = self.sites_data['abbrev']
            sites_constraint = iris.Constraint(
                abbrev=lambda c: c.point in site_abbrevs)
            constraints = constraints & sites_constraint

        if self.site_types is not None:
            types_constraint = iris.Constraint(
                site_type=lambda c: c.point in self.site_types)
            constraints = constraints & types_constraint

        self.sites_cube_list = cubelist.extract(constraints & sname_constraint)

        return self.sites_cube_list

    #----------------------------------------------
    #PRIVATE methods
    # - not expected to be used outside this module
    #----------------------------------------------

    def __read_sites_classes(self):
        """
        Read site classifications file into self.classes_data.
        Also set up self.abbrev2latlon as a dictionary to convert between
        site codes (abbrevs) and (lat,lon) - this will be used to avoid having
        two sites with same site code, but slightly different locations.

        >>> od = CAMSAQObsData()
        >>> od._CAMSAQObsData__read_sites_classes()
        >>> print(len(list(od.abbrev2latlon.keys())))
        4298
        >>> print(od.abbrev2latlon['AL0201A'])
        ('41.3303', '19.8218')
        """
        if self.classes_file is None:
            self.classes_file = config.SAMPLE_DATADIR + \
                                'CAMSAQ_obs/classes_MACC2013.txt'

        #Read data in initially to determine data types
        data = np.genfromtxt(self.classes_file,
                             dtype=None, names=True)

        #Convert string data types to str
        dtypes = []
        for name in data.dtype.names:
            dtype = data.dtype[name]
            if str(dtype)[:2] == '|S':
                dtypes.append((name, str, int(str(dtype)[2:])))
            else:
                dtypes.append((name, dtype))

        #Read in a second time with corrected dtypes
        self.classes_data = np.genfromtxt(self.classes_file,
                                          dtype=dtypes, names=True)

        for site in self.classes_data:
            self.abbrev2latlon[site['code']] = (site['lat'], site['lon'])


    def __data_to_cubelist(self, data):
        """
        Private method.
        Convert data, a numpy ndarray (as read from np.genfromtxt)
        into cubelist - one cube per species
        """
        # Method
        # ------
        # 1. Loop through all data to get a dictionary of site information
        #    (keys are site ids), plus all the datetimes in the data.
        #    These can be used to get the number of unique sites and times.
        # 2. Read in site classifications file to set up site_type in sites_dict
        # 3. Set up an array of size nsites x ntimes for each species
        #    Loop through all data again to populate this array.
        #    Any missing data is left as np.nan
        #    (Not all sites have data for all species)
        # 4. For each species, convert this array to a cube
        #    Add time and site id coordinates, as required by sites_cube_list
        #    Add bounds to the time coordinate to indicate that the point is
        #    at the end of the meaning period
        #    Also add the extra coordinate as setup in the sites_dict
        #    Finally add required attributes and convert to ug/m3.

        #1.
        dt_array, site_id_array, sites_dict = self.__get_sites_dt(data)

        #2.
        sites_dict = self.__set_sites_classes(sites_dict)

        #3.
        species_dict = self.__get_species_dict(data, dt_array, site_id_array)

        #4.
        cubelist = self.__create_cubes(dt_array, site_id_array,
                                       sites_dict, species_dict)

        return cubelist


    def __get_sites_dt(self, data):
        """
        Step 1: Get time and site information ready to set up arrays from.

        Loop through all data to get a dictionary of site information
        (keys are site ids), plus all the datetimes in the data.
        These can be used to get the number of unique sites and times.

        :returns: dt_array, site_id_array, sites_dict

        * dt_array is a numpy array of unique date-times
        * site_id_array is a numpy array of unique site-ids
        * sites_dict is a dictionary of site information with keys
          which are the site id

        """
        dt_list = []
        sites_dict = {}
        for line in data:

            if self.short_name_list is not None:
                short_name = CAMSAQOBS_2_SHORTNAME[line['PARAMETER']]
                if short_name not in self.short_name_list:
                    #Don't need this line
                    continue

            if self.abbrev_prefix_list is not None:
                if not str(line['STATION']).startswith(\
                    tuple(self.abbrev_prefix_list)):
                    #Don't need this site/line
                    continue

            #Fix for particular site which has two different lat/lon locations:
            #Choose location for this station which is present in most
            #observation files and overwrite lat/lons to match this.
            #Overwrite station locations with lat/lon read from classifications
            #file - this should avoid duplication of site locations.
            if line['STATION'] in self.abbrev2latlon:
                line['LAT'] = self.abbrev2latlon[line['STATION']][0]
                line['LON'] = self.abbrev2latlon[line['STATION']][1]

            else:
                self.abbrev2latlon[line['STATION']] = (line['LAT'], line['LON'])

            site_id = adaq_data.generate_siteids([line['LON']],
                                                 [line['LAT']])[0]
            if site_id not in sites_dict:
                sites_dict[site_id] = {'abbrev' : line['STATION'],
                                       'lat' : line['LAT'],
                                       'lon' : line['LON'],
                                       'altitude' : float(line['ALTm'])
                                      }

            dt = datetime.datetime(line['YEAR'], line['MONTH'],
                                   line['DAY'], line['HOUR'])
            dt_list.append(dt)

        #Also ensure that if sites_data is set, then these sites are
        #also included
        if self.sites_data is not None:
            for site in self.sites_data:
                if 'site_id' in site.dtype.names:
                    site_id = site['site_id']
                else:
                    site_id = adaq_data.generate_siteids([site['longitude']],
                                                         [site['latitude']])[0]
                if site_id not in sites_dict:
                    #Not already included, so add this data
                    sites_dict[site_id] = {}
                    for info in site.dtype.names:
                        if info == 'latitude':
                            sites_dict[site_id]['lat'] = site[info]
                        elif info == 'longitude':
                            sites_dict[site_id]['lon'] = site[info]
                        elif info == 'site_altitude':
                            sites_dict[site_id]['altitude'] = site[info]
                        else:
                            sites_dict[site_id][info] = site[info]

        #Get the unique, sorted lists,
        #as numpy arrays of date-times and site-ids
        dt_array = np.array(np.unique(dt_list))
        site_id_array = np.array(sorted(sites_dict.keys()))

        return dt_array, site_id_array, sites_dict


    def __set_sites_classes(self, sites_dict):
        """
        Step 2: Read site classifications file and add site_type (string)
        and site_classes (dictionary, key of species, values are integers)
        to sites_dict for each site.

        :returns: sites_dict - dictionary of site information with keys
                  which are the site id
        """

        if self.classes_data is None:
            self.__read_sites_classes()

        area_conversion_dict = {'urb': 'URBAN',
                                'sub': 'SUBURBAN',
                                'rur': 'RURAL'}
        site_conversion_dict = {'tra': 'TRAFFIC',
                                'bac': 'BACKGROUND',
                                'ind': 'INDUSTRIAL'}

        site_classes_species2shortname = {'o3': 'O3', 'no2': 'NO2',
                                          'pm10': 'PM10', 'nox': 'NOx',
                                          'no': 'NO', 'so2': 'SO2',
                                          'co': 'CO', 'pm25': 'PM2p5'}

        for site_dict in sites_dict.values():
            abbrev = site_dict['abbrev']
            site_indices = np.where(self.classes_data['code'] == abbrev)
            #Check found a match
            if site_indices[0].size:
                #Should only be one matching location
                site_indice = site_indices[0][0]
                site_dict['site_indice'] = site_indice

                #Add site types
                area = self.classes_data['area'][site_indice]
                if area in area_conversion_dict:
                    area = area_conversion_dict[area]
                else:
                    area = 'UNKNOWN'
                site = self.classes_data['site'][site_indice]
                if site in site_conversion_dict:
                    site = site_conversion_dict[site]
                else:
                    site = 'UNKNOWN'
                site_dict['site_type'] = area + '_' + site

                for k, v in site_classes_species2shortname.items():
                    site_dict['site_class_'+v] = int(
                        self.classes_data[k][site_indice])

            else:
                #No match found - set to unknown
                site_dict['site_type'] = 'UNKNOWN'
                for k, v in site_classes_species2shortname.items():
                    site_dict['site_class_'+v] = 0


        return sites_dict


    def __get_species_dict(self, data, dt_array, site_id_array):
        """
        Step 3:
        Set up an array of size nsites x ntimes for each species
        Loop through all data again to populate this array.
        Any missing data is left as np.nan
        (Not all sites have data for all species)

        :returns: species_dict which is a dictionary whose keys are
                  species.
                  Each value is an array of size nsites x ntimes,
                  which is filled with data values.
        """
        species_dict = {}
        for line in data:

            if self.short_name_list is not None:
                short_name = CAMSAQOBS_2_SHORTNAME[line['PARAMETER']]
                if short_name not in self.short_name_list:
                    #Don't need this line
                    continue

            if self.abbrev_prefix_list is not None:
                if not str(line['STATION']).startswith( \
                    tuple(self.abbrev_prefix_list)):
                    #Don't need this site/line
                    continue

            if self.sites_data is not None:
                #Limit by abbreviation/station name
                if line['STATION'] not in self.sites_data['abbrev']:
                    #Don't need this site/line
                    continue

            if line['PARAMETER'] not in species_dict:
                species_dict[line['PARAMETER']] = np.zeros((len(site_id_array),
                                                            len(dt_array)))
                species_dict[line['PARAMETER']][:] = np.nan

            site_id = adaq_data.generate_siteids([line['LON']],
                                                 [line['LAT']])[0]
            i_site_id = np.where(site_id_array == site_id)[0]

            dt = datetime.datetime(line['YEAR'], line['MONTH'],
                                   line['DAY'], line['HOUR'])
            i_dt = np.where(dt_array == dt)[0]

            species_dict[line['PARAMETER']][i_site_id, i_dt] = \
                                                       line['CONCENTRATIONkgm3']

        return species_dict

    def __create_cubes(self, dt_array, site_id_array,
                       sites_dict, species_dict):
        """
        Step 4: Convert to cubes:
        For each species, convert its array to a cube
        Add time and site id coordinates, as required by sites_cube_list
        Add bounds to the time coordinate to indicate that the point is
        at the end of the meaning period
        Also add the extra coordinate as setup in the sites_dict
        Finally add required attributes and convert to ug/m3.

        :returns: cubelist of all species found in species_dict.

        """
        #---
        #Convert to cubes
        time_unit = cf_units.Unit('hours since epoch', calendar='gregorian')
        lat_lon_coord_system = iris.coord_systems.GeogCS(
            semi_major_axis=iris.fileformats.pp.EARTH_RADIUS)

        cubelist = iris.cube.CubeList()
        for param in species_dict.keys():

            short_name = CAMSAQOBS_2_SHORTNAME[param]
            if self.short_name_list is not None:
                if short_name not in self.short_name_list:
                    #Don't need to make a cube for this short name
                    continue

            #Convert site_id and time into coordinates
            site_id_coord = iris.coords.DimCoord(site_id_array,
                                                 long_name='site_id',
                                                 var_name='site_id')
            dt_num_array = time_unit.date2num(dt_array)
            #Add bounds to time coord
            #- the time point is at the end of the meaning period
            #- assume meaning period =1
            bounds = np.array((dt_num_array-1, dt_num_array)).transpose()
            time_coord = iris.coords.DimCoord(dt_num_array,
                                              standard_name='time',
                                              units=time_unit,
                                              bounds=bounds)

            #Create cube
            std_name = CAMSAQOBS_2_STDNAME[param]
            units = 'kg/m3'
            try:
                cube = iris.cube.Cube(species_dict[param],
                                      standard_name=std_name,
                                      units=units,
                                      dim_coords_and_dims=[(site_id_coord, 0),
                                                           (time_coord, 1)])
            # This is a fix for some standard names being removed from CF1.7.
            # If the standard name is no longer in the CF list, use it as a
            # long name instead.  It will still be the name of the cube.
            except ValueError:
                cube = iris.cube.Cube(species_dict[param],
                                      long_name=std_name,
                                      units=units,
                                      dim_coords_and_dims=[(site_id_coord, 0),
                                                           (time_coord, 1)])

            #Add other attributes as coordinates

            #Longitude
            lons = [sites_dict[k]['lon'] for k in site_id_array]
            lon_coord = iris.coords.AuxCoord(lons,
                                             standard_name='longitude',
                                             units='degrees',
                                             coord_system=lat_lon_coord_system)
            cube.add_aux_coord(lon_coord, 0)

            #Latitude
            lats = [sites_dict[k]['lat'] for k in site_id_array]
            lat_coord = iris.coords.AuxCoord(lats,
                                             standard_name='latitude',
                                             units='degrees',
                                             coord_system=lat_lon_coord_system)
            cube.add_aux_coord(lat_coord, 0)

            #Other site coordinates
            abbrevs = [sites_dict[k]['abbrev'] for k in site_id_array]
            cube.add_aux_coord(iris.coords.AuxCoord(
                np.array(abbrevs, dtype=abbrevs[0].dtype),
                long_name='abbrev'), 0)
            cube.add_aux_coord(iris.coords.AuxCoord(
                np.array(abbrevs, dtype=abbrevs[0].dtype),
                long_name='site_name'), 0)
            altitudes = [sites_dict[k]['altitude'] for k in site_id_array]
            cube.add_aux_coord(iris.coords.AuxCoord(
                altitudes, long_name='site_altitude',
                units='m'), 0)

            #Site types
            site_types = [sites_dict[k]['site_type'] for k in site_id_array]
            cube.add_aux_coord(iris.coords.AuxCoord(
                np.array(site_types, dtype=(str, 19)),
                long_name='site_type'), 0)
            #Site classifications
            site_class_species = ['O3', 'NO2', 'PM10', 'NOx',
                                  'NO', 'SO2', 'CO', 'PM2p5']
            for species in site_class_species:
                site_class = [sites_dict[k]['site_class_'+species]
                              for k in site_id_array]
                cube.add_aux_coord(iris.coords.AuxCoord(
                    site_class, long_name='site_classification_'+species),
                                   0)

            #Add attributes
            cube.attributes['short_name'] = short_name
            cube.attributes['label'] = self.label
            cube.attributes['source'] = 'CAMSAQObs'

            #Add cell method
            cell_method = iris.coords.CellMethod('mean', coords='time',
                                                 intervals='1 hour')
            cube.add_cell_method(cell_method)

            #Finally, convert to ug/m3
            cube.convert_units('ug/m3')

            #Add to cubelist
            cubelist.append(cube)

        #Check site names are unique - raise error otherwise as this will
        #cause problems later in adaq python.
        site_names = cubelist[0].coord('site_name').points
        if len(site_names) != len(np.unique(site_names)):
            for site_name in np.unique(site_names):
                if len(np.where(site_names == site_name)[0]) > 1:
                    print('Site name duplicated:', site_name)
            raise IOError('Site names are not unique.')

        return cubelist


def country_from_abbrev(coord, point):
    """
    Return the country string from the abbreviation
    Given by letters upto the first numerical value in abbrev.
    eg if point = 'ES1445A'
    country_str = 'ES'
    """

    new_point = np.array(point[:2], dtype=(str, 4))

    return new_point


def limit_sites(ini_dict, n_sites_per_country=3, limit_species=None,
                output_sites_file=None):
    """
    Limit the number of sites used.
    This is based on the sites with the maximum data availablity for required
    limit_species.
    There is a limit of sites per country, given by n_sites_per_country.
    The generated list of sites is output to file,
    along with a plot of locations.

    :param ini_dict: Dictionary of a :class:`inifile` object
                     Should contain:

                      * 'obs_fmt' - currently only 'camsaqobs' allowed
                      * 'obs_dir' - directory contain obs files
                      * 'short_name_list'
                      * 'start_datetime'
                      * 'end_datetime'
                      * 'plot_dir' - for output of text file and map
                      * 'site_types' - list of required site types

    :param n_sites_per_country: Maximum number of sites per country
    :param limit_species: List of species to use for determining which sites
                          to use. Only sites with data for these species will
                          be considered. If None, then uses short_name_list
                          from ini_dict.
    :param output_sites_file: Filename to write output sites list to. If
                              set to None, then defaults to
                              ini_dict['plot_dir']+'/limited_site_locations.txt'

    :returns: sites_data - numpy ndarray containing site information data
              from a :class:`sites_info.SitesInfo` object.

    >>> ini_dict = {}
    >>> ini_dict['obs_fmt'] = 'camsaqobs'
    >>> ini_dict['obs_dir'] =  config.SAMPLE_DATADIR+'CAMSAQ_obs'
    >>> ini_dict['short_name_list'] = ['O3']
    >>> ini_dict['site_types_list'] = ['RURAL_BACKGROUND', 'URBAN_BACKGROUND']
    >>> ini_dict['start_datetime'] = datetime.datetime(2015,4,4,00)
    >>> ini_dict['end_datetime'] = datetime.datetime(2015,4,5,00)
    >>> ini_dict['plot_dir'] = config.CODE_DIR+'adaqdocs/figures/'
    >>> sites_data = limit_sites(ini_dict, n_sites_per_country=3)
    ... # doctest: +ELLIPSIS
    Creating observation data at ...
    Reading .../obsmacc4ana_20150403.csv
    Reading .../obsmacc4ana_20150404.csv
    Converting to cubes
    Total sites:  51
    Total countries:  19
    Total RURAL_BACKGROUND sites:  35
    Total URBAN_BACKGROUND sites:  16
    Saved to file  .../limited_site_locations.png
    Written to file .../limited_site_locations.txt
    >>> fmtstr = sites_info.format_string(sites_data)
    >>> for site in sites_data[:3]:
    ...    print(fmtstr.format(*site))
       28889.13064500,ES1754A,ES,  40.645,   0.289,     428.0, 0, 0, 0, 0, 1, \
0, 0, 0,ES1754A,RURAL_BACKGROUND
       29092.14229850,GB0045R,GB,  52.298,   0.291,       5.0, 0, 1, 1, 1, 4, \
0, 0, 3,GB0045R,RURAL_BACKGROUND
       44083.13105920,ES1379A,ES,  41.059,   0.441,     363.0, 0, 0, 0, 0, 2, \
0, 0, 0,ES1379A,RURAL_BACKGROUND

    .. image:: ../adaqdocs/figures/limited_site_locations.png
      :scale: 75%

    """

    if limit_species is not None:
        ini_dict['short_name_list'] = limit_species

    od = adaq_functions.get_obs(ini_dict, None)

    #Get example site cube from observation site cube
    sites_cube = od.sites_cube_list[0]


    #Get array of total number of valid times for each site across all species
    n_validtimes = np.zeros((len(sites_cube.coord('site_id').points)))
    for cube in od.sites_cube_list:
        n_validtimes += np.sum(np.isfinite(cube.data), axis=1)

    #Add country coordinate
    coord_cat.add_categorised_coord(sites_cube, 'country_abbrev',
                                    sites_cube.coord('abbrev'),
                                    country_from_abbrev,
                                    units='')

    #Find maximum data coverage in each country
    countries = sites_cube.coord('country_abbrev').points
    unique_countries = np.unique(countries)
    required_abbrevs = []
    for country in unique_countries:
        #Find indices for this country
        country_indices = np.where(countries == country)
        #Get the indices in order which would pick out this country
        #in order of n_validtimes, from min to max
        sort_indices = np.argsort(n_validtimes[country_indices])
        #Indices to actually use for this country
        country_indices = country_indices[0][sort_indices[
            0-n_sites_per_country:]]
        #Now get the abbrevs from a od sites_cube_list
        abbrevs = list(sites_cube[country_indices].coord('abbrev').points)
        required_abbrevs += abbrevs

    new_sites_cube = sites_cube.extract(iris.Constraint(
        abbrev=required_abbrevs))

    #Print some information about chosen sites.
    print('Total sites: ', len(new_sites_cube.coord('abbrev').points))
    print('Total countries: ', len(unique_countries))
    site_types = new_sites_cube.coord('site_type').points
    for site_type in np.unique(site_types):
        print('Total ' + site_type + ' sites: ', end=' ')
        print(len(np.where(site_types == site_type)[0]))

    #Convert into a SitesInfo object so can get sites_data returned,
    #along with plotting locations and writing to a site list
    si = sites_info.SitesInfo()
    sites_data = si.read_from_sites_cube(new_sites_cube)
    if not os.path.isdir(ini_dict['plot_dir']):
        os.makedirs(ini_dict['plot_dir'])
    si.plot_location_map(ini_dict['plot_dir']+'/limited_site_locations.png',
                         label=True, labelsize='x-small')
    if output_sites_file is None:
        output_sites_file = ini_dict['plot_dir']+'/limited_site_locations.txt'
    si.write_to_file(output_sites_file)

    return sites_data


if __name__ == "__main__":

    import doctest
    doctest.testmod()
