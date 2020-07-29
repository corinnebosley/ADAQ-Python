"""
Class for reading from pollen observation files into an ADAQData class.
"""
from __future__ import print_function

import os
import glob
import math
import re
import datetime
import warnings
import pandas as pd
import numpy as np
import iris
import cf_units

import adaq_data

POLLEN_SHORTNAME_2_LATIN = {'grass_pollen': 'Poaceae',
                            'birch_pollen': 'Betula',
                            'oak_pollen'  : 'Quercus'}
POLLEN_SHORTNAME_2_LONGNAME = {k: 'grain_concentration_of_'+v.lower()+'_pollen_in_air'
                               for k, v in POLLEN_SHORTNAME_2_LATIN.items()}

POLLEN_DIAMETER = {'grass_pollen': 35e-6} #meters
POLLEN_DENSITY = {'grass_pollen': 1000.} #kg/m3


def convert_cubelist_units(cubelist, pollen_units='grains m-3'):
    """
    Convert the units of all cubes in a cubelist to and from mass/m3 and
    grains/m3. This can only convert those species which are given in
    POLLEN_KG_PER_GRAIN dictionary.

    :param cubelist: Iris CubeList containing cubes to convert.
    :param pollen_units: Required units to convert to. Currently allowed values
                         are 'grains m-3', 'kg m-3', 'g m-3', 'ug m-3' (or also
                         given as '/m3').
    :returns: Iris cubelist, with units converted.

    Note accessing the input cube list after calling this function will have the
    new units and data.

    """

    #Note, don't use convert_units to grains/m3 as cf_units 'grain'=6.47989e-5kg

    #Ensure consistent naming convention
    pollen_units = pollen_units.replace('/m3', ' m-3')

    for cube in cubelist:

        if cube.units != pollen_units:

            short_name = cube.attributes['short_name']

            if short_name in POLLEN_DIAMETER and short_name in POLLEN_DENSITY:
                # mass (kg) = density (kg/m3) * volume (m3)
                #           = density * (4/3)*pi*r3
                kg_per_grain = POLLEN_DENSITY[short_name] * (4./3.) * math.pi * \
                               (POLLEN_DIAMETER[short_name]/2.)**3

                print('Converting '+short_name+' to '+pollen_units)

                if pollen_units == 'grains m-3':
                    if cube.units in [cf_units.Unit('g m-3'),
                                      cf_units.Unit('ug m-3')]:
                        #Convert to kg/m3
                        cube.convert_units('kg m-3')
                    if cube.units == cf_units.Unit('kg m-3'):
                        #Convert to grains/m3 using mass of a single grain,
                        #dependant on pollen species
                        cube.data = cube.data / kg_per_grain
                        cube.units = 'grains m-3'

                elif pollen_units in ['kg m-3', 'g m-3', 'ug m-3']:
                    cube.data = cube.data * kg_per_grain
                    cube.units = cf_units.Unit('kg m-3')
                    cube.convert_units(pollen_units)

                else:
                    warnings.warn('Unable to convert to '+pollen_units)

    return cubelist


class PollenObsData(adaq_data.ADAQData):
    """
    Subclass of ADAQData, which contains extra functionality specific to pollen
    observations.

    This dataset is given in csv file format, with a single file per site per
    year, which includes all species in grains/m3.
    An example of the expected format is:

    .. code-block:: none

       start_time,end_time,Corylus,Alnus,Salix,Betula,Fraxinus,Ulmus,Quercus,Platanus,Poaceae,\
Urtica,Artemisia,Ambrosia,Notes
       2017-03-14 09:00:00,2017-03-14 11:00:00,13.333,26.667,0.000,0.000,13.333,13.333,0.000,\
0.000,0.000,0.000,0.000,0.000,42 pine pollen
       2017-03-14 11:00:00,2017-03-14 13:00:00,0.000,0.000,0.000,0.000,0.000,13.333,0.000,\
0.000,0.000,0.000,0.000,0.000,42 pine pollen

    **Example:**

    >>> import config
    >>> pollendir = config.SAMPLE_DATADIR + 'pollen_obs/'

    Firstly read in site information data:

    >>> import sites_info
    >>> sitesfile = pollendir + 'PollenSites_short.txt'
    >>> sites = sites_info.SitesInfo()
    >>> sites_data = sites.read_from_file(sitesfile, allsites=True)
    Number of sites:  3

    Set up filenames - this could include wildcards:

    >>> file_pattern = pollendir + '*_bihourly_conc_20*.csv'

    Set up observations data object and read in data for required species and
    time period (if no time period is set, then reads all times in input files):

    >>> od = PollenObsData()
    >>> scl = od.readdata(file_pattern,
    ... short_name_list=['grass_pollen', 'birch_pollen'], sites_data=sites_data,
    ... start_datetime=datetime.datetime(2017,6,1),
    ... end_datetime=datetime.datetime(2017,6,3))

    Examine output (the cube list will be in the order of the requested
    short_name_list):

    >>> print(od.sites_cube_list)
    0: grain_concentration_of_poaceae_pollen_in_air / (grains m-3) (site_id: 3; time: 48)
    1: grain_concentration_of_betula_pollen_in_air / (grains m-3) (site_id: 3; time: 48)

    >>> grass_cube = od.sites_cube_list[0]
    >>> print(grass_cube)
    grain_concentration_of_poaceae_pollen_in_air / (grains m-3) (site_id: 3; time: 48)
         Dimension coordinates:
              site_id                                                   x        -
              time                                                      -        x
         Auxiliary coordinates:
              abbrev                                                    x        -
              latitude                                                  x        -
              longitude                                                 x        -
              site_altitude                                             x        -
              site_name                                                 x        -
              site_type                                                 x        -
         Attributes:
              label: Obs
              short_name: grass_pollen
              source: PollenObs
         Cell methods:
              mean: time (1 hour)

    Note this contains information for all sites in the sites file, although not
    all of these observation files may exist (eg Bath):

    >>> print(grass_cube.coord('site_name').points)
    ['Exeter' 'Eskdalemuir' 'Bath']

    Note as the raw data is measured bihourly (in local time, although the
    data in the csv files are in UTC), the assumption is that the concentrations
    are consistent for the two hours. For example 26.667 is valid for the period
    10Z - 12Z in the following example. Note also the time points are given at
    the end of each bounded period to be consistent with other datasets.

    >>> print(grass_cube[0].data[10:20])
    [26.667 26.667 13.333 13.333 66.667 66.667 40.    40.     0.     0.   ]

    >>> print(grass_cube.coord('time')[10:20])
    DimCoord([2017-06-01 11:00:00, 2017-06-01 12:00:00, 2017-06-01 13:00:00,
           2017-06-01 14:00:00, 2017-06-01 15:00:00, 2017-06-01 16:00:00,
           2017-06-01 17:00:00, 2017-06-01 18:00:00, 2017-06-01 19:00:00,
           2017-06-01 20:00:00], bounds=[[2017-06-01 10:00:00, 2017-06-01 11:00:00],
           [2017-06-01 11:00:00, 2017-06-01 12:00:00],
           [2017-06-01 12:00:00, 2017-06-01 13:00:00],
           [2017-06-01 13:00:00, 2017-06-01 14:00:00],
           [2017-06-01 14:00:00, 2017-06-01 15:00:00],
           [2017-06-01 15:00:00, 2017-06-01 16:00:00],
           [2017-06-01 16:00:00, 2017-06-01 17:00:00],
           [2017-06-01 17:00:00, 2017-06-01 18:00:00],
           [2017-06-01 18:00:00, 2017-06-01 19:00:00],
           [2017-06-01 19:00:00, 2017-06-01 20:00:00]], standard_name='time', calendar='gregorian')

    Can also change the units to match NAME output of g/m3:

    >>> print(np.nanmax(od.sites_cube_list[0].data),od.sites_cube_list[0].units)
    80.0 grains m-3
    >>> cubelist = convert_cubelist_units(od.sites_cube_list, pollen_units='g/m3')
    Converting grass_pollen to g m-3
    >>> print('{:.2e}'.format(np.nanmax(od.sites_cube_list[0].data)),od.sites_cube_list[0].units)
    1.80e-06 g m-3

    And also test the conversion back to grains/m3:

    >>> cubelist = convert_cubelist_units(od.sites_cube_list, pollen_units='grains/m3')
    Converting grass_pollen to grains m-3
    >>> print('{:.1f}'.format(np.nanmax(od.sites_cube_list[0].data)),od.sites_cube_list[0].units)
    80.0 grains m-3

    """

    def __init__(self, label='Obs'):
        """
        Initiates class as a subset of adaq_data.ADAQData,
        plus other observation-specific information.
        """
        adaq_data.ADAQData.__init__(self)
        self.label = label #Label
        self.short_name_list = None #List of short names
        self.start_datetime = None #Starting datetime
        self.end_datetime = None #End datetime
        self.filenames = None #List of raw data filenames
        self.df_dict = {} #Dictionary of raw data frames, with keys of site names
        self.df_start_dt = None #Earliest start time available from raw data frames
        self.df_end_dt = None #Latest end time available from raw data frames
        self.tcoord = None #Time coordinate used in all cubes
        self.sites_data = None #Numpy array containing site information data
        self.sort_site_indices = None #Numpy array giving the index of the
                                      #sites_data array within the cube

    def readdata(self, file_pattern=None, short_name_list=None, sites_data=None,
                 start_datetime=None, end_datetime=None):
        """
        Create sites_cube_list from files specified.

        :param file_pattern: String, files to read in, may contain wildcards.
        :param short_name_list: List of short  names.
        :param start_datetime:  Start time of data in datetime format.
        :param end_datetime: End time of data in datetime format.

        :returns: sites_cube_list
        """

        if file_pattern is None:
            raise ValueError('No file_pattern given')
        if short_name_list is not None:
            self.short_name_list = short_name_list
        if self.short_name_list is None:
            raise ValueError('No short_name_list given')
        if start_datetime is not None:
            self.start_datetime = start_datetime
        if end_datetime is not None:
            self.end_datetime = end_datetime
        if sites_data is not None:
            self.sites_data = sites_data
        if self.sites_data is None:
            raise ValueError('No sites_data given')


        #Read raw dataframes in from files
        self.__read_files(file_pattern)

        if self.df_dict:

            #Create blank (nan data) cubes
            for short_name in self.short_name_list:
                if short_name in POLLEN_SHORTNAME_2_LATIN:
                    sites_cube = self.__create_blank_cube(short_name)
                    self.sites_cube_list.append(sites_cube)

            #Fill cubes in with data from input files
            self.__fill_cubes()

        return self.sites_cube_list


    def __read_files(self, file_pattern):
        """
        Read files into a dictionary of pandas dataframes, whose keys
        are site names. If multiple files (eg with different dates), for
        the same site, then these are appended onto the same dataframe.
        Note, this expects filenames (after the directory), to start
        with the site name.

        :param file_pattern: String, files to read in, may contain wildcards.

        >>> import config
        >>> pollendir = config.SAMPLE_DATADIR + 'pollen_obs/'
        >>> import sites_info
        >>> sitesfile = pollendir + 'PollenSites_short.txt'
        >>> sites = sites_info.SitesInfo()
        >>> sites_data = sites.read_from_file(sitesfile, allsites=True)
        Number of sites:  3

        >>> file_pattern = pollendir + 'Exeter_bihourly_conc_20*.csv'
        >>> od = PollenObsData()
        >>> od.sites_data = sites_data
        >>> od._PollenObsData__read_files(file_pattern)
        >>> print(od.df_dict['Exeter']) # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
                       start_time             end_time  ...  Ambrosia            Notes
        0     2017-03-14 09:00:00  2017-03-14 11:00:00  ...       0.0   42 pine pollen
        1     2017-03-14 11:00:00  2017-03-14 13:00:00  ...       0.0   42 pine pollen
        2     2017-03-14 13:00:00  2017-03-14 15:00:00  ...       0.0   42 pine pollen
        ...
        2049  2017-09-04 02:00:00  2017-09-04 04:00:00  ...       0.0              NaN
        2050  2017-09-04 04:00:00  2017-09-04 06:00:00  ...       0.0              NaN
        2051  2017-09-04 06:00:00  2017-09-04 08:00:00  ...       0.0              NaN
        <BLANKLINE>
        [2052 rows x 15 columns]

        This routine also sets the earliest start date and latest end date
        found within the files into the variables self.df_start_dt, and
        self.df_end_dt respectively.

        >>> print(od.df_start_dt, od.df_end_dt)
        2017-03-14 09:00:00 2017-09-04 08:00:00
        """

        self.filenames = []
        sitenames = []
        for filename in glob.glob(file_pattern):
            for site_name in self.sites_data['site_name']:
                if re.search('^'+site_name, os.path.basename(filename)):
                    self.filenames.append(filename)
                    sitenames.append(str(site_name))

        for site_name, filename in zip(sitenames, self.filenames):
            df = pd.read_csv(filename)
            #Get start and end time of the file data
            df_start_dt = datetime.datetime.strptime(df['start_time'].min(), '%Y-%m-%d %H:%M:%S')
            df_end_dt = datetime.datetime.strptime(df['end_time'].max(), '%Y-%m-%d %H:%M:%S')
            #Don't keep this data if none of it is in required time range
            if self.start_datetime is not None:
                if df_end_dt < self.start_datetime:
                    continue
            if self.end_datetime is not None:
                if df_start_dt > self.end_datetime:
                    continue
            #Store earliest start and latest end time from all files
            if self.df_start_dt is None or df_start_dt < self.df_start_dt:
                self.df_start_dt = df_start_dt
            if self.df_end_dt is None or df_end_dt > self.df_end_dt:
                self.df_end_dt = df_end_dt
            #Store in dictionary
            if site_name not in self.df_dict:
                self.df_dict[site_name] = df
            else:
                self.df_dict[site_name] = self.df_dict[site_name].append(df)

        if not self.df_dict:
            warnings.warn('No pollen observation files with required dates in')



    def __create_blank_cube(self, short_name):
        """
        Build up a blank cube full of nan data for required short name.

        :param short_name: String, short name to be given to the cube. This is
                           also used to translate to long name which includes
                           the latin name.

        :returns: Iris cube, which is a valid ADAQData site cube, all data is
                  set to np.nan

        >>> import config
        >>> pollendir = config.SAMPLE_DATADIR + 'pollen_obs/'
        >>> import sites_info
        >>> sitesfile = pollendir + 'PollenSites_short.txt'
        >>> sites = sites_info.SitesInfo()
        >>> sites_data = sites.read_from_file(sitesfile, allsites=True)
        Number of sites:  3
        >>> od = PollenObsData()
        >>> od.df_start_dt = datetime.datetime(2018,5,1)
        >>> od.df_end_dt = datetime.datetime(2018,5,3)
        >>> od.sites_data = sites_data
        >>> cube = od._PollenObsData__create_blank_cube('grass_pollen')
        >>> print(cube)
        grain_concentration_of_poaceae_pollen_in_air / (grains m-3) (site_id: 3; time: 48)
             Dimension coordinates:
                  site_id                                                   x        -
                  time                                                      -        x
             Auxiliary coordinates:
                  abbrev                                                    x        -
                  latitude                                                  x        -
                  longitude                                                 x        -
                  site_altitude                                             x        -
                  site_name                                                 x        -
                  site_type                                                 x        -
             Attributes:
                  label: Obs
                  short_name: grass_pollen
                  source: PollenObs
             Cell methods:
                  mean: time (1 hour)

        """

        #Set up time coord
        if self.tcoord is None:
            #Only need to create this once

            start_dt = (self.df_start_dt if self.start_datetime is None
                        else max([self.df_start_dt, self.start_datetime]))
            end_dt = (self.df_end_dt if self.end_datetime is None
                      else min([self.df_end_dt, self.end_datetime]))
            #Points are at the end of each bounded period
            dt_array = pd.date_range(start_dt+datetime.timedelta(hours=1),
                                     end_dt, freq='1H').to_pydatetime()
            time_unit = cf_units.Unit('hours since epoch', calendar='gregorian')
            dt_num_array = time_unit.date2num(dt_array)
            bounds = np.array((dt_num_array-1, dt_num_array)).transpose()
            self.tcoord = iris.coords.DimCoord(dt_num_array,
                                               standard_name='time',
                                               units=time_unit,
                                               bounds=bounds)

        #Set up sites id coord
        site_id_array = np.array(adaq_data.generate_siteids(
            self.sites_data['longitude'], self.sites_data['latitude']))
        #Get indices such that the site_id_array will be monotonic and sorted
        self.sort_site_indices = np.argsort(site_id_array)
        site_id_coord = iris.coords.DimCoord(
            site_id_array[self.sort_site_indices],
            long_name='site_id',
            var_name='site_id')

        #Set up nan data
        data = np.full((len(site_id_array), len(self.tcoord.points)), np.nan)

        #Set up cube
        cube = iris.cube.Cube(data,
                              long_name=POLLEN_SHORTNAME_2_LONGNAME[short_name],
                              dim_coords_and_dims=[(site_id_coord, 0),
                                                   (self.tcoord, 1)],
                              attributes={'short_name': short_name,
                                          'label': self.label,
                                          'source': 'PollenObs'},
                              units='grains m-3')

        #Add other site attributes as coordinates
        lat_lon_coord_system = iris.coord_systems.GeogCS(
            semi_major_axis=iris.fileformats.pp.EARTH_RADIUS)
        lat_coord = iris.coords.AuxCoord(
            self.sites_data['latitude'][self.sort_site_indices],
            standard_name='latitude',
            units='degrees',
            coord_system=lat_lon_coord_system)
        cube.add_aux_coord(lat_coord, 0)

        lon_coord = iris.coords.AuxCoord(
            self.sites_data['longitude'][self.sort_site_indices],
            standard_name='longitude',
            units='degrees',
            coord_system=lat_lon_coord_system)
        cube.add_aux_coord(lon_coord, 0)

        cube.add_aux_coord(iris.coords.AuxCoord(
            self.sites_data['abbrev'][self.sort_site_indices],
            long_name='abbrev'), 0)
        cube.add_aux_coord(iris.coords.AuxCoord(
            self.sites_data['site_name'][self.sort_site_indices],
            long_name='site_name'), 0)
        cube.add_aux_coord(iris.coords.AuxCoord(
            self.sites_data['site_altitude'][self.sort_site_indices],
            long_name='site_altitude', units='m'), 0)
        cube.add_aux_coord(iris.coords.AuxCoord(
            self.sites_data['site_type'][self.sort_site_indices],
            long_name='site_type'), 0)

        #Add cell method
        cube.add_cell_method(
            iris.coords.CellMethod('mean', coords='time', intervals='1 hour'))

        return cube

    def __fill_cubes(self):
        """
        Use the raw data to fill in the blank cubes with data, for all sites and
        required species.

        :returns: Iris cube list of sites cubes.

        Tested as part of the PollenObsData() class doctests.
        """

        #Loop through sites in the same order as they appear in the cubes
        for isite, site in enumerate(self.sites_data[self.sort_site_indices]):

            if str(site['site_name']) not in self.df_dict:
                continue
            df = self.df_dict[str(site['site_name'])]


            tbounds_dt = self.tcoord.units.num2date(self.tcoord.bounds)

            #Fill in cubes with data
            for index, row in df.iterrows():
                start_dt = datetime.datetime.strptime(row['start_time'],
                                                      '%Y-%m-%d %H:%M:%S')
                end_dt = datetime.datetime.strptime(row['end_time'],
                                                    '%Y-%m-%d %H:%M:%S')
                #Find matching times in cube
                cube_dt_indices_start = np.where((tbounds_dt[:, 0] >= start_dt))[0]
                cube_dt_indices_end = np.where((tbounds_dt[:, 1] <= end_dt))[0]
                if cube_dt_indices_start.size and cube_dt_indices_end.size:
                    #Match found - set the value for all data points
                    #Assume a constant concentration for all of these hours

                    #Fill cubes for all required short_names
                    for sp_cube in self.sites_cube_list:
                        short_name = sp_cube.attributes['short_name']
                        latin_name = POLLEN_SHORTNAME_2_LATIN[short_name]

                        sp_cube.data[isite,
                                     cube_dt_indices_start[0]:cube_dt_indices_end[-1]+1] = \
                                     row[latin_name]

        return self.sites_cube_list


if __name__ == '__main__':

    import doctest
    doctest.testmod()
