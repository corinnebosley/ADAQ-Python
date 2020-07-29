"""
Contains AURNData class for reading in AURN and MARGA observation data both
fixed format and .csv
"""
from __future__ import division
from __future__ import print_function

from six.moves.builtins import str
from six.moves.builtins import zip

import datetime
import os
import glob
import csv
import warnings
import numpy as np

import iris
import iris.iterate
import cf_units

import adaq_data

#: AURN_2_STDNAME is used to get correct standard/long name from the
#: short name. This is specific to AURN and MARGA data  and if for example
#: phenomena in AURN files changed, this would need to be revisited
AURN_2_STDNAME = {
    'O3' : 'mass_concentration_of_ozone_in_air',
    'NO'          : 'mass_concentration_of_nitrogen_monoxide_in_air',
    'NO2'         : 'mass_concentration_of_nitrogen_dioxide_in_air',
    'CO'          : 'mass_concentration_of_carbon_monoxide_in_air',
    'PM10'        : 'mass_concentration_of_pm10_ambient_aerosol_in_air',
    'PM2p5'       : 'mass_concentration_of_pm2p5_ambient_aerosol_in_air',
    'SO2'         : 'mass_concentration_of_sulfur_dioxide_in_air',
    'jNO2'        : 'photolysis_rate_of_nitrogen_dioxide',
    'C5H8'        : 'mass_concentration_of_isoprene_in_air',
    'NOx_as_NO2'  : ('mass_concentration_of_nox_'
                     'expressed_as_nitrogen_dioxide_in_air'),
    'PM10_Ca2'    : ('mass_concentration_of_calcium_'
                     'in_pm10_ambient_aerosol_in_air'),
    'PM10_Cl'     : ('mass_concentration_of_chloride_'
                     'in_pm10_ambient_aerosol_in_air'),
    'PM10_K'      : ('mass_concentration_of_potassium_'
                     'in_pm10_ambient_aerosol_in_air'),
    'PM10_Mg2'    : ('mass_concentration_of_magnesium_'
                     'in_pm10_ambient_aerosol_in_air'),
    'PM10_NA'     : ('mass_concentration_of_sodium_in_'
                     'pm10_ambient_aerosol_in_air'),
    'PM10_NH4'    : ('mass_concentration_of_ammonium_'
                     'in_pm10_ambient_aerosol_in_air'),
    'PM10_NO3'    : ('mass_concentration_of_nitrate_'
                     'in_pm10_ambient_aerosol_in_air'),
    'PM10_SO42'   : ('mass_concentration_of_sulphate_'
                     'in_pm10_ambient_aerosol_in_air'),
    'PM2p5_Ca2'   : ('mass_concentration_of_calcium_'
                     'in_pm2p5_ambient_aerosol_in_air'),
    'PM2p5_Cl'    : ('mass_concentration_of_chloride_'
                     'in_pm2p5_ambient_aerosol_in_air'),
    'PM2p5_K'     : ('mass_concentration_of_potassium_'
                     'in_pm2p5_ambient_aerosol_in_air'),
    'PM2p5_Mg2'   : ('mass_concentration_of_magnesium_'
                     'in_pm2p5_ambient_aerosol_in_air'),
    'PM2p5_NA'    : ('mass_concentration_of_sodium_'
                     'in_pm2p5_ambient_aerosol_in_air'),
    'PM2p5_NH4'   : ('mass_concentration_of_ammonium_'
                     'in_pm2p5_ambient_aerosol_in_air'),
    'PM2p5_NO3'   : ('mass_concentration_of_nitrate_'
                     'in_pm2p5_ambient_aerosol_in_air'),
    'PM2p5_SO4_2' : ('mass_concentration_of_sulphate_in_'
                     'pm2p5_ambient_aerosol_in_air'),
    'HCl'         : 'mass_concentration_of_hydrochloric_acid_in_air',
    'HNO2'        : 'mass_concentration_of_nitrous_acid_in_air',
    'HNO3'        : 'mass_concentration_of_nitric_acid_in_air',
    'NH3'         : 'mass_concentration_of_ammonia_in_air'
    }

# Conversion between AURN and MARGA parameters and short_name
# This list is not exhaustive and there may be additional species available from
# data e.g. VOCs such as ethane and additional species such as black carbon.
# Note + and - are not included on the ends of any short_names as these are
# lost on reading in using genfromtxt
# These would also need to be added to AURN_2_STDNAME and
# AURN_LONGNAME_2_SHORTNAME to be utilised.
AURN_LONGNAME_2_SHORTNAME = {
    'Ozone'                                     : 'O3',
    'Nitric oxide'                              : 'NO',
    'Nitrogen dioxide'                          : 'NO2',
    'Carbon monoxide'                           : 'CO',
    'PM10 particulate matter (Hourly measured)' : 'PM10',
    'PM10 Particulate matter (Hourly)'          : 'PM10',
    'PM2.5 particulate matter (Hourly measured)': 'PM2p5',
    'PM2.5 Particulate matter (Hourly)'         : 'PM2p5',
    'Sulphur dioxide'                           : 'SO2',
    'isoprene'                                  : 'C5H8',
    'Nitrogen oxides as nitrogen dioxide'       : 'NOx_as_NO2',
    'Nitrogen oxides'                           : 'NOx_as_NO2',
    'calcium in PM10'                           : 'PM10_Ca2',
    'chloride in PM10'                          : 'PM10_Cl',
    'potassium in PM10'                         : 'PM10_K',
    'magnesium in PM10"'                        : 'PM10_Mg2',
    'sodium in PM10'                            : 'PM10_NA',
    'ammonium in PM10'                          : 'PM10_NH4',
    'nitrate in PM10'                           : 'PM10_NO3',
    'sulphate in PM10'                          : 'PM10_SO42',
    'calcium in PM2.5'                          : 'PM2p5_Ca2',
    'chloride in PM2.5'                         : 'PM2p5_Cl',
    'potassium in PM2.5'                        : 'PM2p5_K',
    'magnesium in PM2.5'                        : 'PM2p5_Mg2',
    'sodium in PM2.5'                           : 'PM2p5_NA',
    'ammonium in PM2.5'                         : 'PM2p5_NH4',
    'nitrate in PM2.5'                          : 'PM2p5_NO3',
    'sulphate in PM2.5'                         : 'PM2p5_SO4_2',
    'gaseous hydrochloric acid'                 : 'HCl',
    'gaseous nitrous acid'                      : 'HNO2',
    'gaseous nitric acid'                       : 'HNO3',
    'gaseous ammonia'                           : 'NH3',
    'gaseous sulphur dioxide'                   : 'SO2'
    }

def ffdate2dt(ffdate):
    """ Converts from fixed format (%Y%m%d) to datetime format """
    return datetime.datetime.strptime(ffdate.decode(), "%Y%m%d")

def fftime2tdelta(time):
    """ Converts from fixed format (HH) to time-delta """
    return datetime.timedelta(hours=int(time.decode()))

def ffdate2dtcsv(ffdate):
    """ Converts from fixed format (%Y-%m-%d) or (%d-%m-%Y) to
    datetime format """
    ffdate = ffdate.decode() #Convert to unicode not bytes
    yyyy = ffdate.split('-')[0]
    if len(yyyy) == 4:
        dt = datetime.datetime.strptime(ffdate, "%Y-%m-%d")
    elif len(yyyy) == 2:
        dt = datetime.datetime.strptime(ffdate, "%d-%m-%Y")

    return dt

def fftime2tdeltacsv(time):
    """ Converts from fixed format (HH:MM:SS) or (HH:MM) to
    time-delta for marga data taking (HH)"""
    return datetime.timedelta(hours=int(time.decode().split(':')[0]))

def genfromtxt_csv(filename, headers, dtypes=None):
    """
    Reads data in - this is a numpy ndarray, which will
    allow columns to be accessed by their names, eg data['NH3']

    :param filename: Input filename to read
    :param headers: List of header names to overwrite exisiting headers in file
    :param dtypes: numpy.dtype to describe data type for each header
    :returns: data, numpy ndarray
    """
    try:
        data = np.genfromtxt(filename,
                             dtype=dtypes,
                             skip_header=5,
                             delimiter=',',
                             names=headers,
                             missing_values='No data',
                             converters={'date_time':ffdate2dtcsv,
                                         'time':fftime2tdeltacsv})

    # If a value error is returned, the last line of the file i.e.
    # in the case of MARGA data where the last line includes
    # text will be excluded from the numpy ndarray.
    except ValueError:
        data = np.genfromtxt(filename,
                             dtype=dtypes,
                             skip_header=5,
                             skip_footer=1,
                             delimiter=',',
                             names=headers,
                             missing_values='No data',
                             converters={'date_time':ffdate2dtcsv,
                                         'time':fftime2tdeltacsv})

    return data

class AURNData(adaq_data.ADAQData):

    """
    Subclass of ADAQData, which contains extra functionality specific to AURN
    observation data. This will read both fixed format and .csv files for AURN
    and MARGA data.
    Converts all data into a sites_cube_list.

    Usage:

    Get sites information:

    >>> import sites_info
    >>> import config
    >>> obsdir = config.SAMPLE_DATADIR+'AURN_obs/'
    >>> sitefile = obsdir+'aq_sites_GEMSRAQ_v4b_dev.txt'
    >>> sites = sites_info.SitesInfo()
    >>> sites_data = sites.read_from_file(sitefile,
    ... allsites=False, obsdir=obsdir)
    Number of sites:  5

    Setup dates:

    >>> start_dt = datetime.datetime(2014, 8, 10, 00)
    >>> end_dt = datetime.datetime(2014, 8, 20, 00)

    Setup AURN Data object:

    >>> od = AURNData()

    To read data from raw fixed-format AURN observation files and immediately
    convert to a sites_cube_list:

    >>> scl = od.readdata(obsdir, short_name_list=['O3','NO2'],
    ... sites=sites_data, start_datetime=start_dt,
    ... end_datetime=end_dt) # doctest: +ELLIPSIS
    Reading obs data files
    Found obs for  .../ABD_....txt
    Found obs for  .../ACTH_....txt
    Found obs for  .../AH_....txt
    Found obs for  .../HAR_....txt
    Found obs for  .../YW_....txt
    Creating obs cubes

    Examine observation data object:

    >>> print(od) # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    <class '....AURNData'> - Subclass of ADAQData - Contains:
    sites_cube_list:
    0: mass_concentration_of_nitrogen_dioxide_in_air / (ug/m3) \
(site_id: 5; time: 241)
    1: mass_concentration_of_ozone_in_air / (ug/m3) \
(site_id: 5; time: 241)
    gridded_cube_list:
    < No cubes >
    trajectory_cube_list:
    < No cubes >

    To then extract an AURNData object containing only data with abbrev='HAR'
    and short_name='O3':

    >>> har_od = od.extract(short_name='O3',abbrev='HAR')
    >>> print(har_od) # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    <class '....AURNData'> - Subclass of ADAQData - Contains:
    sites_cube_list:
    0: mass_concentration_of_ozone_in_air / (ug/m3) (time: 241)
    gridded_cube_list:
    < No cubes >
    trajectory_cube_list:
    < No cubes >

    This actually contains just one cube, so could be retrieved as a
    single cube:

    >>> har_cube = od.extract(short_name='O3',abbrev='HAR',singlecube=True)
    >>> print(har_cube) #doctest: +ELLIPSIS
    mass_concentration_of_ozone_in_air / (ug/m3) (time: 241)
         Dimension coordinates:
              time                                    x
         Scalar coordinates:
              abbrev: HAR
              latitude: 51.57110977 degrees
              longitude: -1.326666594 degrees
              site_altitude: 137 m
              site_id: 35867333.14...
              site_name: Harwell
              site_type: RURAL
         Attributes:
              label: Obs
              short_name: O3
              source: AURN
         Cell methods:
              mean: time (1 hour)

    To output this data now to fixed format files:

    >>> outputdir = config.TEST_DIR + '/aurn_fixed_format'
    >>> od.write_ff_files(outputdir) # doctest: +ELLIPSIS
    Written to .../aurn_fixed_format/2014_fixed_format/ABD_2014.txt
    Written to .../aurn_fixed_format/2014_fixed_format/ACTH_2014.txt
    Written to .../aurn_fixed_format/2014_fixed_format/AH_2014.txt
    Written to .../aurn_fixed_format/2014_fixed_format/HAR_2014.txt
    Written to .../aurn_fixed_format/2014_fixed_format/YW_2014.txt

    Check the output of one of these files:

    >>> outputfile = outputdir + '/2014_fixed_format/HAR_2014.txt'
    >>> with open(outputfile,"r") as fin:
    ...     for iline, line in enumerate(fin):
    ...         print(line.strip())
    ...         if iline > 3: break
    locationID date     time       O3_ug/m3   NO2_ug/m3
    HAR        20140810 00         40.00000     7.10000
    HAR        20140810 01         42.00000     4.40000
    HAR        20140810 02              nan         nan
    HAR        20140810 03         31.00000     3.30000


    The code can also read in from multiple directories and merge all the
    data from different files to a single observations cube. For example:

    >>> obsdir_list = [config.SAMPLE_DATADIR+'AURN_obs/raw_csv',
    ... config.SAMPLE_DATADIR+'AURN_obs/MARGA_obs/fixed_csv']

    >>> od = AURNData()
    >>> scl = od.readdata(obsdir_list=obsdir_list,
    ... short_name_list=['O3','NO2','PM10_NO3'],
    ... sites=sites_data,
    ... start_datetime=datetime.datetime(2015,1,1,0),
    ... end_datetime=datetime.datetime(2015,1,3,0)) # doctest: +ELLIPSIS
    Reading obs data files
    Found obs for  .../AURN_obs/raw_csv/ABD_2015.csv
    Found obs for  .../AURN_obs/raw_csv/ACTH_2015.csv
    Found obs for  .../AURN_obs/MARGA_obs/fixed_csv/ACTH_MARGA_2015.csv
    Found obs for  .../AURN_obs/raw_csv/AH_2015.csv
    Found obs for  .../AURN_obs/raw_csv/HAR_2015.csv
    Found obs for  .../AURN_obs/MARGA_obs/fixed_csv/HAR_MARGA_2015.csv
    Found obs for  .../AURN_obs/raw_csv/YW_2015.csv
    Creating obs cubes
    Creating missing data cubes

    >>> print(od) # doctest: +ELLIPSIS
    <class '....AURNData'> - Subclass of ADAQData - Contains:
    sites_cube_list:
    0: mass_concentration_of_nitrate_in_pm10_ambient_aerosol_in_air / \
(ug/m3) (site_id: 5; time: 49)
    1: mass_concentration_of_nitrogen_dioxide_in_air / (ug/m3) \
(site_id: 5; time: 49)
    2: mass_concentration_of_ozone_in_air / (ug/m3) (site_id: 5; time: 49)
    gridded_cube_list:
    < No cubes >
    trajectory_cube_list:
    < No cubes >

    """

    def __init__(self, label='Obs'):

        """
        Initiates class as a subset of adaq_data.ADAQData, plus other
        observation-specific data
        """

        adaq_data.ADAQData.__init__(self)
        self.obsdir_list = None #List of directories contain obs,
                           #first directory is highest priority to get data from
        self.obsdir_fixed_csv = None  #Directory containing fixed csv files
        self.raw_data = dict()  #Raw data read in from observation files
        self.raw_units = dict()
        self.utc_format = '%Y-%m-%d %H:%M:%S'
        self.start_datetime = None
        self.end_datetime = None
        self.short_name_list = None
        self.sites = None
        self.label = label

    def readdata(self, obsdir=None, obsdir_list=None,
                 short_name_list=None, sites=None,
                 start_datetime=None, end_datetime=None, obsdir_fixed_csv=None):

        """
        Method to reads in raw data from fixed-format files and converts
        individual sites to sitecubes before returns a cubelist containing
        one cube per species, of which each cube is in sitecube format.

        :param obsdir: Directory containing observations.
        :param obsdir_list: List of directories containing observations.
                            First directory in list is highest priority for
                            getting observations from, last directory is lowest
                            priority so will only use data from here if not
                            available in any other directories.
        :param short_name_list: List of short names.
        :param sites: numpy ndarray containing site information data
                       from a :class:`sites_info.SitesInfo` object.
        :param start_datetime: Datetime format for start time.
        :param end_datetime: Datetime format for end time.
        """
        self.obsdir_list = obsdir_list #Directories containing observations
        self.obsdir_fixed_csv = obsdir_fixed_csv
        self.start_datetime = start_datetime
        self.end_datetime = end_datetime
        self.short_name_list = short_name_list
        self.sites = sites

        if self.short_name_list is None:
            raise ValueError("AURNData: no short_name_list has been set")
        if self.sites is None:
            raise ValueError("AURNData: sites has not been set")
        if self.obsdir_list is None:
            self.obsdir_list = obsdir
        if isinstance(self.obsdir_list, str):
            #Convert to a list if only a single directory passed in
            self.obsdir_list = [self.obsdir_list]
        if not self.obsdir_list:
            raise ValueError("AURNData: No obsdir_list has been set")
        if self.start_datetime > self.end_datetime:
            raise IOError('Start time is after end time')

        print('Reading obs data files')
        self.raw_data, self.raw_units = self.__read_alldata()

        print('Creating obs cubes')
        allcubes = iris.cube.CubeList()

        #Set up dictionary to hold site information
        #about sites which do not have any data - this could be because
        #this species is not observed at this site, or because there
        #is no observations file for this site.
        #Dictionary keys are species, values are list of sites (each
        #site is a dictionary containing site information).
        missingsites = {}

        for site in self.sites:

            for species in self.short_name_list:
                data_for_species = False

                if self.raw_data[site['abbrev']]: #has data

                    datacubes = iris.cube.CubeList()
                    for data in self.raw_data[site['abbrev']]:

                        if species in data.dtype.names:

                            data_for_species = True
                            datacube = self.__create_cube(species, site, data)
                            datacubes.append(datacube)

                    if data_for_species:
                        #Merge cubes together.
                        #First cube has highest priority,
                        #so set mergedcube to this initially
                        mergedcube = datacubes[0]
                        mergedcubedata = mergedcube.data[0, :]
                        for cube in datacubes[1:]:
                            #Only one site per cube, so safe to extract
                            #just first site and all times
                            cubedata = cube.data[0, :]
                            for i, (mergedvalue, cubevalue) in \
                                enumerate(zip(mergedcubedata, cubedata)):
                                #Note all cubes have same number of times,
                                #so same length of data
                                if not np.isfinite(mergedvalue):
                                    #The merged cube has nan value, so overwrite
                                    #with value from current cube
                                    mergedcubedata[i] = cubevalue

                        allcubes.append(mergedcube)

                if not data_for_species:

                    #No observations for this site/species
                    if species in missingsites:
                        missingsites[species].append(site)
                    else:
                        missingsites[species] = [site]

        if missingsites:
            #Need to add some nan cubes in
            print('Creating missing data cubes')
            for species, sites in missingsites.items():
                if species in AURN_2_STDNAME:
                    cube_name = AURN_2_STDNAME[species]
                    #Get a single example cube for this species
                    # to base new cubes on
                    species_cube_list = allcubes.extract(
                        iris.Constraint(cube_name))
                    if species_cube_list:
                        #This species is available for at least one site
                        species_cube = species_cube_list[0]
                        #Create a new cube for each missing site
                        for site in sites:
                            cube = self.__create_blank_cube(site, species_cube)
                            allcubes.append(cube)


        self.sites_cube_list = allcubes.concatenate()


        return self.sites_cube_list


    def write_ff_files(self, outputdir):
        """
        Output AURNdata to fixed-format (space-delimited) file

        :param outputdir: Directory to output files to.
                          Subdirectories will be created below this level
                          to output different years to different directories.

        Output filenames will be of the form:
        outputdir/year_fixed_format/siteabbrev_year.txt
        """


        if not os.path.isdir(outputdir):
            os.makedirs(outputdir)

        #Set up timearray
        #Cubes should all have same times in them
        tcoord = self.sites_cube_list[0].coord('time')
        tpts_dt = np.array([tcoord.units.num2date(dt) for dt in tcoord.points])

        #Set up outputdata dictionary
        #keys will be site abbrev, values will be 2d arrays of (time, shortname)
        outputdata = {}
        for abbrev in self.sites['abbrev']:
            outputdata[abbrev] = np.full((len(tpts_dt),
                                          len(self.short_name_list)),
                                         np.nan)

        #Set up headers array ready for outputting
        headers = np.zeros((len(self.short_name_list)+3), dtype=(str, 12))
        headers[:3] = ['locationID', 'date', 'time']


        for speciescube in self.sites_cube_list:
            short_name = speciescube.attributes['short_name']
            i_short_name = np.where(np.array(self.short_name_list)
                                    == short_name)[0][0]
            #Fill in remaining headers
            units = str(speciescube.units)
            headers[i_short_name+3] = short_name + '_' + units
            for sitecube in speciescube.slices_over('site_id'):
                abbrev = sitecube.coord('abbrev').points[0]
                #Check time coord is as expected
                assert sitecube.coord('time') == tcoord
                outputdata[abbrev][:, i_short_name] = sitecube.data

        #Now start outputting data

        years = np.array([dt.year for dt in tpts_dt])
        for year in np.unique(years):
            #Set up output subdirectory
            outputdir_year = outputdir+'/'+str(year)+'_fixed_format/'
            if not os.path.isdir(outputdir_year):
                os.makedirs(outputdir_year)

            for abbrev, dataarray in sorted(outputdata.items()):
                #Check if some data exists for this year
                indices = np.where(years == year)
                if np.sum(np.isfinite(dataarray[indices[0], :])) == 0:
                    #If all nan, then don't output file
                    print(abbrev + ' - all nan data, no file output')
                    continue
                filename = outputdir_year + abbrev + '_' + str(year) + '.txt'
                with open(filename, 'w') as fout:
                    header_fmt = ("%-10s %-8s %-7s"
                                  + len(self.short_name_list)*" %11s" + '\n')
                    fout.write(header_fmt % tuple(headers))

                    for idt, dt in enumerate(tpts_dt):

                        if dt.year == year:
                            #LocationID
                            output_str = "%-10s " % abbrev

                            #Date Time columns
                            dt_str = dt.strftime("%Y%m%d %H")
                            output_str += dt_str + '     '

                            #Loop over individual short_names
                            #and output their data
                            fmt = len(self.short_name_list)*" %11.5f"
                            output_str += fmt % tuple(dataarray[idt, :])

                            fout.write(output_str+ '\n')

                print('Written to '+filename)


    #----------------------------------------------
    #PRIVATE methods
    # - not expected to be used outside this module
    #----------------------------------------------

    def __add_site_coords(self, cube, site):
        """
        Add site-specific coordinates to a cube
        """

        #Add dim coord of site_id
        site_id = adaq_data.generate_siteids([site['longitude']],
                                             [site['latitude']])[0]
        cube.add_aux_coord(iris.coords.DimCoord(
            site_id, long_name='site_id', var_name='site_id'))
        cube = iris.util.new_axis(cube, 'site_id')

        #Set up lat/lon coords and add as aux coords
        lat_lon_coord_system = iris.coord_systems.GeogCS(
            semi_major_axis=iris.fileformats.pp.EARTH_RADIUS)

        lon_coord = iris.coords.DimCoord([site['longitude']],
                                         standard_name='longitude',
                                         units='degrees',
                                         coord_system=lat_lon_coord_system)

        lat_coord = iris.coords.DimCoord([site['latitude']],
                                         standard_name='latitude',
                                         units='degrees',
                                         coord_system=lat_lon_coord_system)

        cube.add_aux_coord(lat_coord, 0)
        cube.add_aux_coord(lon_coord, 0)

        #Add other site coordinates to cube as aux coords
        cube.add_aux_coord(iris.coords.AuxCoord(
            np.array(site['site_name'], dtype=self.sites['site_name'].dtype),
            long_name='site_name'), 0)
        if 'abbrev' in site.dtype.names:
            cube.add_aux_coord(iris.coords.AuxCoord(
                np.array(site['abbrev'], dtype=self.sites['abbrev'].dtype),
                long_name='abbrev'), 0)
        if 'site_type' in site.dtype.names:
            cube.add_aux_coord(iris.coords.AuxCoord(
                np.array(site['site_type'],
                         dtype=self.sites['site_type'].dtype),
                long_name='site_type'), 0)
        if 'site_altitude' in site.dtype.names:
            cube.add_aux_coord(iris.coords.AuxCoord(
                site['site_altitude'], long_name='site_altitude',
                units='m'), 0)

        return cube


    def __create_cube(self, species, site, data):
        """
        Put all data into a sites_cube for this particular site and species
        """
        abbrev = site['abbrev']

        #Deal with time
        #- get date/times from data and convert these to an iris coord.
        time_array = [datetime.datetime.strptime(str(dt), self.utc_format)
                      for dt in data['date_time']]
        time_unit = cf_units.Unit('hours since epoch', calendar='gregorian')
        time_array = time_unit.date2num(time_array)
        time_coord = iris.coords.DimCoord(time_array, standard_name='time',
                                          units=time_unit)
        #Add bounds - points is at the end of meaning time
        #So 1-2Z meaning period is recorded at point of 2Z.
        if len(time_array) > 1:
            time_coord.guess_bounds(bound_position=1.0)

        #Set up cube
        data_array = data[species]
        cube = iris.cube.Cube(data_array)

        # set name for this species using dictionary above
        # use cube.rename to set long name if not
        # a valid standard name
        if species in AURN_2_STDNAME:
            cube.rename(AURN_2_STDNAME[species])
        else:
            raise ValueError(species + " not included in AURN_2_STDNAME")

        cube.units = self.raw_units[abbrev][species]

        #Add time as a dim coord
        cube.add_dim_coord(time_coord, 0)
        cell_method = iris.coords.CellMethod('mean', coords='time',
                                             intervals='1 hour')
        cube.add_cell_method(cell_method)

        #Add site coords
        cube = self.__add_site_coords(cube, site)

        #Add attributes
        cube.attributes['short_name'] = species
        cube.attributes['label'] = self.label
        cube.attributes['source'] = 'AURN'

        return cube


    def __create_blank_cube(self, site, species_cube):
        """
        Create a blank cube for a particular site and species
        where its data is filled with nans.
        This uses an example cube for this species to set up
        times, units and cube name.
        """

        time_coord = species_cube.coord('time')

        data_array = np.full((len(time_coord.points)), np.nan)

        #Set up cube
        cube = iris.cube.Cube(data_array)
        cube.rename(species_cube.name())
        cube.units = species_cube.units

        #Add time coords
        cube.add_dim_coord(time_coord, 0)
        cell_method = iris.coords.CellMethod('mean', coords='time',
                                             intervals='1 hour')
        cube.add_cell_method(cell_method)

        #Add site coords
        cube = self.__add_site_coords(cube, site)

        #Add extra attributes
        cube.attributes['short_name'] = species_cube.attributes['short_name']
        cube.attributes['label'] = self.label
        cube.attributes['source'] = 'AURN'

        return cube


    def __quality_control(self, data):
        """
        Quality control to remove negative data and replace with nan.
        """
        #Quality control - remove negative data
        for name in data.dtype.names:
            #To avoid case sensitivity use lower case version of name
            if name.lower() != 'locationid' and name != 'date_time' \
            and name.lower() != 'time':
                #Only need to check if all the data is not nan
                if data[name].dtype == 'float64':
                    if not np.isnan(data[name]).all():
                        indices = np.where(data[name] < 0)
                        if indices[0].size:
                            data[name][indices[0]] = np.nan

        return data


    def __fix_dates(self, data):
        """
        If a start date and end date are not specified use values
        from data. Create a np array of this length replacing any missing
        data with nan.
        """
        #Get start and end times from self,
        #otherwise use values from data.
        if self.start_datetime is None:
            start_datetime = data['date_time'][0]
        else:
            start_datetime = self.start_datetime
        if self.end_datetime is None:
            end_datetime = data['date_time'][-1]
        else:
            end_datetime = self.end_datetime


        #If a certain period is required (set in self), then may need to
        #extend the range of dates given in data.
        #Note not all files have the same dates in
        #(eg if a station is turned on/off mid-year)
        #To do this, loop over all required dates, and create new
        #ndarray with these dates in and the data from original data array
        #Any dates without data from original nd-data array are set
        #to have nan data values.
        if self.start_datetime is not None and self.end_datetime is not None:
            #Required number of hours
            nsecs = (self.end_datetime-self.start_datetime).total_seconds()
            nhours = np.int(nsecs/3600)+1
            #Set up required data np array
            reqdata = np.full((nhours,), np.nan, dtype=data.dtype)
            dt = start_datetime
            i = 0
            while dt <= end_datetime:
                reqdata[i]['date_time'] = dt

                #Fill in with data where available from data array
                indices = np.where(data['date_time'] == dt)
                if indices[0].size:

                    if not data.shape:
                        #Only a single time contained within data
                        reqdata[[i]] = data
                    else:
                        reqdata[[i]] = data[indices[0]]
                dt += datetime.timedelta(hours=1)
                i += 1

            #Overwrite original data array.
            data = reqdata
        return data


    def __read_ffheader(self, filename):
        """
        Read header from filename and return appropriate header names:
          * _ugm3 and _mgm3 are removed
          * date is replaced with date_time
          * PM2.5 is replaced with PM2p5 as python has problems with the '.'

        :returns: header, dtypes, units, empty_file
                  Where header is is a list of sensible header names;
                  dtypes is a list of data types for each column;
                  units is a dictionary whose keys are species names, and values
                  are the unit for that species;
                  empty_file is a logical which is set to False if no data
                  is contained in file.

        Example header format:
        "locationID date     time    CO_mgm3   NO2_ugm3    O3_ugm3
        NO_ugm3  PM10_ugm3 PM2.5_ugm3   SO2_ugm3  C5H8_ugm3"
        """
        if not os.path.exists(filename):
            raise ValueError("File %s not found" % filename)

        empty_file = False
        with open(filename, "r") as fin:
            line = fin.readline()
            headers = line.split()
            units = {}
            dtypes = []
            for i, header in enumerate(headers):
                if header == 'date':
                    headers[i] = 'date_time'
                    dtype = 'O'
                elif header == 'time':
                    dtype = 'O'
                elif header == 'locationID':
                    dtype = (str, 10)
                elif len(header.split('_')) > 1:
                    #Remove _ugm3/_mgm3
                    headersplit = header.split('_')
                    #Remove '.' as not recognised by python
                    # - replace with p - eg change PM2.5 to PM2p5
                    for ihs, hitem in enumerate(headersplit):
                        headersplit[ihs] = hitem.replace('.', 'p')
                    headers[i] = '_'.join(headersplit[0:-1])
                    units[headers[i]] = headersplit[-1]
                    #Convert units to udunits format
                    if units[headers[i]] == 'ugm3':
                        units[headers[i]] = 'ug/m3'
                    elif units[headers[i]] == 'mgm3':
                        units[headers[i]] = 'mg/m3'
                    dtype = np.float64
                dtypes.append(dtype)
            nextline = fin.readline()
            if not nextline:
                empty_file = True
                print('Warning: no data in file ' + filename)

        return headers, dtypes, units, empty_file

    def __read_fffile(self, filename):
        """
        Read all data from fixed format file,
        only returning required daterange.
        Any negative data is removed (set to nan).
        Note, the code assumes all data should be hourly means.

        Example file format:
        locationID date     time    CO_mgm3   NO2_ugm3    O3_ugm3    \
NO_ugm3  PM10_ugm3 PM2.5_ugm3   SO2_ugm3  C5H8_ugm3
        YW         20140101 01          NaN      0.900     77.000      \
0.100        NaN        NaN        NaN        NaN
        YW         20140101 02          NaN      1.600     74.000      \
0.200        NaN        NaN        NaN        NaN
        YW         20140101 03          NaN      1.200     73.000      \
0.100        NaN        NaN        NaN        NaN

        :param filename: Input filename to read
        :returns: data, units (explained in more detail through example below)

        >>> od = AURNData()
        >>> import config
        >>> filename = config.SAMPLE_DATADIR+'AURN_obs/YW_20140101_20140818.txt'
        >>> data, units = od._AURNData__read_fffile(filename)
        ... # doctest: +ELLIPSIS
        Found obs for  .../YW_20140101_20140818.txt

        Here data is an ndarray, which has names based on the headers
        in the file:

        >>> print(data.dtype.names)
        ('locationID', 'date_time', 'time', 'CO', 'NO2', 'O3', \
'NO', 'PM10', 'PM2p5', 'SO2', 'C5H8')

        Each row is a different time for the same site, eg from
        the above file:

        >>> print(data[0]) #doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        ('YW', datetime.datetime(2014, 1, 1, 1, 0), \
datetime.timedelta(0, 3600),  nan,  0.9,  ......  0.1,  nan,  nan,  nan,  nan)

        All the values from a single header can be accessed in a similar
        way to a dictionary, eg:

        >>> np.set_printoptions(formatter={'float':lambda x: '{:5.2f}'.format(x)})
        >>> print(data['O3']) #doctest: +ELLIPSIS
        [77.00 74.00 73.00 ...... 58.00 57.00 53.00]
        >>> np.set_printoptions()

        The returned value of units is a dictionary of the units used
        by each species, eg for the above file,

        >>> units = {'CO': 'mg/m3', 'PM10': 'ug/m3', 'NO': 'ug/m3',
        ... 'C5H8': 'ug/m3', 'PM2p5': 'ug/m3', 'SO2': 'ug/m3',
        ... 'O3': 'ug/m3', 'NO2': 'ug/m3'}

        """

        if not os.path.exists(filename):
            raise ValueError("File %s not found" % filename)
        headers, dtypes, units, empty_file = self.__read_ffheader(filename)

        if empty_file:
            return np.empty((0)), [] #data, units

        #Read data in - this is a numpy ndarray, which will
        #allow columns to be accessed by their names, eg data['PM10']
        data = np.genfromtxt(filename, dtype=dtypes, names=headers, skip_header=1,
                             converters={'date_time':ffdate2dt,
                                         'time':fftime2tdelta})

        #Change date to be date+time
        data['date_time'] = data['date_time']+data['time']

        data = self.__quality_control(data)
        data = self.__fix_dates(data)


        if not data.size:
            print('Required dates not found for ', filename)
        else:
            print('Found obs for ', filename)

        return data, units

    def __read_csvheader(self, filename):
        """
        Read header from csv filename and return appropriate header names:
          * The header is on the 4th line in the .csv file
          * This header line is split where there is a comma separation
          * Date is replaced with date_time
          * Time is replaced with time
          * If the header entry is 'Status' this is converted to 'Units' and
          then each species is appended e.g. Units_NO

        Example header format - Marga.csv:
        Hourly measurement data supplied by UK-air on 17/11/2016
        All Data GMT hour ending
        Status: V=Verified P=Provisionaly Verified N=Not Verified S=Suspect
        ,,"Harwell",,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
        Date,Time,"Nitric oxide",Status,"Nitrogen dioxide",Status,    \
        "calcium in PM10",Status,"chloride in PM10",Status,    \
        "potassium in PM10",Status,"magnesium in PM10",Status,    \
        "sodium in PM10",Status, "ammonium in PM10",Status,    \
        "nitrate in PM10",Status,"sulphate in PM10",Status,    \
        "calcium in PM2.5",Status,"chloride in PM2.5",Status,    \
        "potassium in PM2.5",Status,"magnesium in PM2.5",Status,    \
        "sodium in PM2.5",Status,"ammonium in PM2.5",Status,    \
        "nitrate in PM2.5",Status,"sulphate in PM2.5",Status,    \
        "gaseous hydrochloric acid",Status,"gaseous nitrous acid",    \
        Status,"gaseous nitric acid",Status,"gaseous ammonia",Status,    \
        "gaseous sulphur dioxide",Status

        Example header format - AURN.csv:
        Data supplied by AEA on 25/3/2016
        All Data GMT hour ending
        Status: R =Ratified P=Provisional,P*=As supplied
        ,,Yarner Wood,,,,,,,,,,,,,,,,,,,,,,,
        Date,time,"Nitric oxide",status,unit,"Nitrogen dioxide",status,    \
        unit,"Nitrogen oxides as nitrogen dioxide",status,unit,"Ozone",    \
        status,unit

        :param filename: Input filename to read
        :returns: newheaders, dtypes, empty_file
                  Where newheaders are a list of headers which will include any
                  renamed headers
                  dtypes is a list of (header, datatype)
                  empty_file is a logical which is set to True if file is empty
                  or only contains header information (ie empty of data)
        """
        with open(filename, "r") as fin:
            reader = csv.reader(fin, delimiter=',')
            for iline, line in enumerate(reader):
                if iline == 4:
                    headers = line
                if iline > 6:
                    break
            if iline <= 6:
                empty_file = True
                print('Warning: no data in file ' + filename)
            else:
                empty_file = False
        #The dytpes list will contain a list of tuples containing
        #(header, datatype) the data types for each header.
        #The newheaders list will include a list of headers
        #which will include any renamed unused duplicates.
        dtypes = []
        newheaders = []
        for iheader, header in enumerate(headers):
            header = header.replace(",", "_")
            header = header.replace("<sub>", "")
            header = header.replace("</sub>", "")
            if header == 'Date':
                header = 'date_time'
                dtype = 'O'
            if header == 'Time':
                header = 'time'
                dtype = 'O'
            if header in AURN_LONGNAME_2_SHORTNAME:
                short_name = AURN_LONGNAME_2_SHORTNAME[header]
                if short_name in newheaders:
                    #There may be duplicate columns in some .csv data
                    #If this occurs the last duplicate column will be used
                    #Therefore, the previous instance of the name is changed
                    #in addition to the matching Units column.
                    indicies = np.where(np.array(newheaders) == short_name)
                    if indicies[0].size:
                        index = indicies[0][0]
                        newheaders[index] += '_unused' + str(index)
                        dtypes[index] = (newheaders[index], np.float64)
                    indicies = np.where(np.array(newheaders) ==
                                        'Units_' + short_name)
                    if indicies[0].size:
                        index = indicies[0][0]
                        newheaders[index] += '_unused' + str(index)
                        dtypes[index] = (newheaders[index], (str, 17))
                header = AURN_LONGNAME_2_SHORTNAME[header]
                dtype = np.float64
            #For Marga .csv the 'Status' column has both Status of the
            #observation and units e.g. 'V ugm-3'
            if header == 'Status':
                header = 'Units_' + newheaders[iheader-1]
                dtype = (str, 17)
            if header == 'status':
                header = header + '_' + str(iheader)
                dtype = (str, 1)
            #For AURN .csv there are separate 'status' and 'units' columns
            #e.g. the units column may include 'ugm-3 (TEOM FDMS)'
            if header == 'unit':
                header = 'Units_' + newheaders[iheader-2]
                dtype = (str, 17)
            dtypes.append((header, dtype))
            newheaders.append(header)

        return newheaders, dtypes, empty_file

    def __fix_csvfile(self, filename):
        """
        This reads in the original .csv file and inserts a comma in the final
        column of each line if it is missing i.e. No data will become No data,
        This modified output_site_file will then be output to the directory:
        self.obsdir_fixed_csv.
        """
        site_filename = os.path.basename(filename) #removes directory component
        output_site_file = self.obsdir_fixed_csv+'/'+site_filename
        if not os.path.isdir(self.obsdir_fixed_csv):
            os.makedirs(self.obsdir_fixed_csv)
        with open(filename, "r") as fin:
            with open(output_site_file, 'w') as fout:
                for line in fin:
                    if line[-8:] == 'No data\n':
                        line = line[:-1] + ',\n'
                    fout.write(line)
        print('created new file: ', output_site_file)

        return output_site_file

    def __read_csvfile(self, filename):
        """
        Read all data from csv format file,
        only returning required daterange.
        A new csv file may be produced if there is no comma at the end of the
        line in the file.
        Note, the code assumes all data should be hourly means.

        Example file format AURN csv:
        Data supplied by AEA on 25/3/2016
        All Data GMT hour ending
        Status: R =Ratified P=Provisional,P*=As supplied
        ,,Yarner Wood,,,,,,,,,,,,,,,,,,,,,,,
        Date,time,"Nitric oxide",status,unit,"Nitrogen dioxide",status,    \
        unit,"Nitrogen oxides as nitrogen dioxide",status,unit,"Ozone",    \
        status,unit

        01-01-2015,01:00,0.19084,R,ugm-3,2.42760,R,    \
        ugm-3,2.72021,R,ugm-3,78.89667,R,ugm-3

        :param filename: Input filename to read
        :returns: data, units (explained in more detail through example below)

        >>> od = AURNData()
        >>> import config
        >>> filename = config.SAMPLE_DATADIR+'AURN_obs/raw_csv/YW_2015.csv'
        >>> data, units = od._AURNData__read_csvfile(filename)
        ... # doctest: +ELLIPSIS
        Found obs for  .../YW_2015.csv

        Here data is an ndarray, which has names based on the headers
        in the file:

        >>> print(data.dtype.names) #doctest: +NORMALIZE_WHITESPACE
        ('date_time', 'time', 'NO', 'status_3', 'Units_NO', 'NO2', \
'status_6', 'Units_NO2', 'NOx_as_NO2', 'status_9', 'Units_NOx_as_NO2', 'O3', \
'status_12', 'Units_O3')

        Each row is a different time for the same site, eg from
        the above file:

        >>> print(data[0]) #doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        (datetime.datetime(2015, 1, 1, 1, 0), datetime.timedelta(0, 3600), \
0.19..., 'R', 'ugm-3', 2.42..., 'R', 'ugm-3', 2.72..., 'R', 'ugm-3', \
78.89..., 'R', 'ugm-3')

        All the values from a single header can be accessed in a similar
        way to a dictionary, eg:

        >>> np.set_printoptions(formatter={'float':lambda x: '{:5.2f}'.format(x)})
        >>> print(data['O3']) #doctest: +ELLIPSIS
        [78.90 77.03 79.58 ...... 66.16 64.41 53.43]

        The returned value of units is a dictionary of the units used
        by each species, eg for the above file,

        >>> units == {'O3': 'ug/m3', 'NOx_as_NO2': 'ug/m3',
        ... 'NO2': 'ug/m3', 'NO': 'ug/m3'}
        True

        Example file format MARGA csv:
        Hourly measurement data supplied by UK-air on 17/11/2016
        All Data GMT hour ending
        Status: V=Verified P=Provisionaly Verified N=Not Verified S=Suspect
        ,,"Harwell",,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
        Date,Time,"Nitric oxide",Status,"Nitric oxide",Status, \
"Nitrogen dioxide",Status,"Nitrogen dioxide",Status,"calcium in PM10", \
Status,"chloride in PM10",Status,"potassium in PM10",Status, \
"magnesium in PM10",Status,"sodium in PM10",Status,"ammonium in PM10", \
Status,"nitrate in PM10",Status,"sulphate in PM10",Status, \
"calcium in PM2.5",Status,"chloride in PM2.5",Status, \
"potassium in PM2.5",Status,"magnesium in PM2.5",Status, \
"sodium in PM2.5",Status,"ammonium in PM2.5",Status,"nitrate in PM2.5", \
Status,"sulphate in PM2.5",Status,"gaseous hydrochloric acid",Status, \
"gaseous nitrous acid",Status,"gaseous nitric acid",Status,"gaseous ammonia", \
Status,"gaseous sulphur dioxide",Status
        2015-01-01,01:00:00,No data,,0.32554,V ugm-3,No data,,2.87258,V ugm-3, \
0.003,V ugm-3,1.878,V ugm-3,0.211,V ugm-3,0.123,V ugm-3,1.013,V ugm-3, 2.7, \
V ugm-3,5.006,V ugm-3,4.126,V ugm-3,0.008,V ugm-3,1.081,V ugm-3,0.126,V ugm-3, \
0.099,V ugm-3,0.429,V ugm-3,2.479,V ugm-3,3.915,V ugm-3,3.857,V ugm-3,0.021, \
V ugm-3,0.332,V ugm-3,0.15,V ugm-3,0.816,V ugm-3, 0.144,V ugm-3

        >>> od = AURNData()
        >>> filename = (config.SAMPLE_DATADIR +
        ... 'AURN_obs/MARGA_obs/raw_csv/HAR_MARGA_2015.csv')
        >>> od.obsdir_fixed_csv = config.TEST_DIR + '/fixed_csv'
        >>> data, units = od._AURNData__read_csvfile(filename)
        ... # doctest: +ELLIPSIS
        created new file:  .../fixed_csv/HAR_MARGA_2015.csv
        Found obs for .../fixed_csv/HAR_MARGA_2015.csv

        Here data is an ndarray, which has names based on the headers
        in the file:

        >>> print(data.dtype.names)
        ('date_time', 'time', 'NO_unused2', 'Units_NO_unused3', 'NO', \
'Units_NO', 'NO2_unused6', 'Units_NO2_unused7', 'NO2', 'Units_NO2', \
'PM10_Ca2', 'Units_PM10_Ca2', 'PM10_Cl', 'Units_PM10_Cl', 'PM10_K', \
'Units_PM10_K', 'magnesium_in_PM10', 'Units_magnesium_in_PM10', 'PM10_NA', \
'Units_PM10_NA', 'PM10_NH4', 'Units_PM10_NH4', 'PM10_NO3', 'Units_PM10_NO3', \
'PM10_SO42', 'Units_PM10_SO42', 'PM2p5_Ca2', 'Units_PM2p5_Ca2', 'PM2p5_Cl', \
'Units_PM2p5_Cl', 'PM2p5_K', 'Units_PM2p5_K', 'PM2p5_Mg2', 'Units_PM2p5_Mg2', \
'PM2p5_NA', 'Units_PM2p5_NA', 'PM2p5_NH4', 'Units_PM2p5_NH4', 'PM2p5_NO3', \
'Units_PM2p5_NO3', 'PM2p5_SO4_2', 'Units_PM2p5_SO4_2', 'HCl', 'Units_HCl', \
'HNO2', 'Units_HNO2', 'HNO3', 'Units_HNO3', 'NH3', 'Units_NH3', 'SO2', \
'Units_SO2')

        Each row is a different time for the same site, eg from
        the above file:

        >>> print(data[0])
        (datetime.datetime(2015, 1, 1, 1, 0), datetime.timedelta(0, 3600), \
  nan, '',  0.33, 'V ugm-3',   nan, '',  2.87, 'V ugm-3',  0.00, 'V ugm-3', \
 1.88, 'V ugm-3',  0.21, 'V ugm-3', '0.123', 'V ugm-3',  1.01, 'V ugm-3', \
 2.70, 'V ugm-3',  5.01, 'V ugm-3',  4.13, 'V ugm-3',  0.01, 'V ugm-3', \
 1.08, 'V ugm-3',  0.13, 'V ugm-3',  0.10, 'V ugm-3',  0.43, 'V ugm-3', \
 2.48, 'V ugm-3',  3.92, 'V ugm-3',  3.86, 'V ugm-3',  0.02, 'V ugm-3', \
 0.33, 'V ugm-3',  0.15, 'V ugm-3',  0.82, 'V ugm-3',  0.14, 'V ugm-3')

        All the values from a single header can be accessed in a similar
        way to a dictionary, eg:

        >>> print(data['NH3']) #doctest: +ELLIPSIS
        [ 0.82  0.82  0.88 ......  0.51  0.50  0.50]
        >>> np.set_printoptions()

        The returned value of units is a dictionary of the units used
        by each species, eg for the above file,

        >>> units == {'PM10_NH4': 'ug/m3', 'PM2p5_NA': 'ug/m3',
        ... 'NO2': 'ug/m3', 'PM10_NO3': 'ug/m3', 'NO': 'ug/m3',
        ... 'PM10_NA': 'ug/m3', 'PM2p5_NH4': 'ug/m3', 'NO2_unused7': 'ug/m3',
        ... 'magnesium_in_PM10': 'ug/m3', 'NO_unused3': 'ug/m3',
        ... 'HNO3': 'ug/m3', 'HNO2': 'ug/m3', 'PM2p5_Ca2': 'ug/m3',
        ... 'PM10_Cl': 'ug/m3', 'NH3': 'ug/m3', 'PM2p5_K': 'ug/m3',
        ... 'PM2p5_SO4_2': 'ug/m3', 'PM2p5_NO3': 'ug/m3', 'PM10_K': 'ug/m3',
        ... 'PM2p5_Cl': 'ug/m3', 'PM10_Ca2': 'ug/m3', 'SO2': 'ug/m3',
        ... 'PM10_SO42': 'ug/m3', 'HCl': 'ug/m3', 'PM2p5_Mg2': 'ug/m3'}
        True
        """

        headers, dtypes, empty_file = self.__read_csvheader(filename)

        if empty_file:
            return np.empty((0)), [] #data, units

        #Read all data in - this is a numpy ndarray, which will
        #allow columns to be accessed by their names, eg data['NH3']
        #A try-except statement is used herebecause of possible problems
        #with commas missing at the end of lines in the file causing
        # errors with reading file in directly
        try:
            data = genfromtxt_csv(filename, headers, dtypes)
        except:
            if self.obsdir_fixed_csv is None:
                self.obsdir_fixed_csv = self.obsdir_list[0] + '/fixed_csv'
            try:
                filename = self.__fix_csvfile(filename)
            except IOError:
                raise IOError('Unable to write to '+ self.obsdir_fixed_csv)
            data = genfromtxt_csv(filename, headers, dtypes)

        ## Quality control - check that all units in a Status column are the
        ## same and return units
        units = {}
        for name in data.dtype.names:
            if "Units_" in name:
                species = name[6:]
                speciesunits = None
                uniqueunits = np.unique(data[name])
                for uniqueunit in uniqueunits:
                    if uniqueunit[-5:] == 'ugm-3': #'V ugm-3'
                        speciesunits = 'ug/m3'
                    elif uniqueunit[-5:] == 'mgm-3': #'V mgm-3'
                        speciesunits = 'mg/m3'
                    elif uniqueunit[:-12] == 'ugm-3': #'ugm-3 (TEOM FDMS)'
                        speciesunits = 'ug/m3'
                    elif uniqueunit[:-12] == 'mgm-3': #'mgm-3 (TEOM FDMS)'
                        speciesunits = 'mg/m3'
                if speciesunits is None:
                    #Assume by default to be ug/m3 for most species
                    speciesunits = 'ug/m3'
                    if species == 'CO':
                        #CO is ususally given in mg/m3 so set this as default
                        speciesunits = 'mg/m3'
                units[species] = speciesunits

        #Change date to be date+time
        data['date_time'] = data['date_time']+data['time']

        data = self.__quality_control(data)
        data = self.__fix_dates(data)

        if data.size:
            print('Found obs for ', filename)
        else:
            print('Required dates not found for ', filename)

        return data, units

    def __read_alldata(self):
        """
        Read all data in for all sites

        :returns: alldata, allunits.

        Where alldata is a dictionary, whose keys are site abbreviations
        and whose values for each site is a list of an ndarrays of data,
        one ndarray for each successfully read in file.
        (see description in the example for __read_fffile)

        Where allunits is a dictionary, whose keys are site abbreviations,
        and whose values for each site is a dictionary of the units
        for each species.
        """

        alldata = {}
        allunits = {}
        for site in self.sites:
            alldata[site['abbrev']] = []
            for obsdir in self.obsdir_list:
                filenames = sorted(glob.glob(obsdir+'/'+site['abbrev']+'_*'))
                if filenames:
                    for filename in filenames:
                        #Get file extension
                        fname, fext = os.path.splitext(filename)
                        if fext == ".txt":
                            data, units = self.__read_fffile(filename)
                        elif fext == ".csv":
                            data, units = self.__read_csvfile(filename)
                        else:
                            warnings.warn('Unknown file extension: '+ fext)
                            data, units = np.empty((0)), []

                        if data.size:
                            #Add to list of data
                            alldata[site['abbrev']].append(data)
                            #Add to units dictionary
                            if site['abbrev'] not in allunits:
                                allunits[site['abbrev']] = units
                            else:
                                #Check that all units are included
                                for key, value in units.items():
                                    if key not in allunits[site['abbrev']]:
                                        allunits[site['abbrev']][key] = value
            if not alldata[site['abbrev']]:
                print("No file found for "+site['abbrev'])

        return alldata, allunits

if __name__ == '__main__':

    import doctest
    doctest.testmod()
