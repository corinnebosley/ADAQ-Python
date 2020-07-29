"""
Functions that work on ADAQData objects
"""
from __future__ import division
from __future__ import print_function

import datetime
import warnings
import os
import copy
import glob
import six

from six.moves.builtins import zip
from six.moves.builtins import range

import iris
import numpy as np

#Import required modules from adaqcode
import inifile
import sites_info
import adaq_data
import config
import timeseries_stats
import cube_statistics
import cube_chemistry
import cube_time
import aq_indices

#Import individual data classes for use in get_obs and get_models
#Observations:
import aurn_data
import camsaqobs_data
import pollenobs_data
#Models:
import ecgrib_data
import maccens_data
import name_data
import nimrod_data
import pp_data
import trajectory_data


def create_example_data(exampletype='1days', save=True):
    """
    Create example data files of gridded_cube_lists and sites_cube_lists.

    :param exampletype: Type and/or number of days of data to retrieve.
                        Currently set up:

                        * **1days** - short 1 day example, only has two species
                        * **5days** - fully 5 day example, with 5 species
                        * **3d** - 3d data, 3 hours of data with one species
                          (ozone) plus pressure and temperature
                        * **stats** - 10 days of 2 species saved as statistics
                          cubes
                        * **10days** - longer 10 day period, with only 2
                          species, spread across two months

    :param save: Save the created data - True, or False. False can be
                 useful if testing.

    :returns: (ini_dict, sites_data, od, md_list) - or for exampletype='stats':
              (ini_dict, stats_cubes)

    .. note:: The doctests here are designed to test that the code has not
              changed the sample data. If these need changing, then the sample
              data files should also be updated!

    To create sample data for use, need to run with each possible exampletype.
    Note running with exampletype='stats' should be done after
    exampletype='full'.
    The resulting files should then be put into config.SAMPLE_DATADIR.

    >>> ini_dict, sites_data, od, md_list = create_example_data(
    ... save=False) # doctest: +ELLIPSIS
    Reading inifile .../example_data_1days.ini
    Number of sites:  5
    Creating observation data at  ...
    Reading obs data files
    Found obs for .../ABD_20140101_20140818.txt
    Found obs for .../ACTH_20140101_20140818.txt
    Found obs for .../AH_20140101_20140818.txt
    Found obs for .../HAR_20140101_20140818.txt
    Found obs for .../YW_20140101_20140818.txt
    Creating obs cubes
    Getting model data for  aqum_oper  at ...
    Getting model data for  aqum_casestudy  at ...
    >>> print(sorted(ini_dict.keys())) # doctest: +NORMALIZE_WHITESPACE
    ['end_date', 'end_datetime',
    'forecast_day', 'models_dir_list', 'models_fmt_list', 'models_list',
    'obs_dir', 'obs_fmt', 'range_days', 'short_name_list',
    'sites_file', 'start_date', 'start_datetime']
    >>> fmtstr = sites_info.format_string(sites_data)
    >>> for site in sites_data:
    ...     print(fmtstr.format(*site))
    GB0001, ABD,  57.158, -2.094,  20,URBAN_BACKGROUND,        Aberdeen
    GB0003,ACTH,  55.883, -3.347, 260,           RURAL,Auchencorth_Moss
    GB0002,  AH,  52.504, -3.034, 370,           RURAL,      Aston_Hill
    GB0045, HAR,  51.571, -1.327, 137,           RURAL,         Harwell
    GB0128,  YW,  50.597, -3.716, 119,           RURAL,     Yarner_Wood
    >>> print(od.sites_cube_list)
    0: mass_concentration_of_ozone_in_air / (ug/m3) \
(site_id: 5; time: 25)
    1: mass_concentration_of_pm2p5_ambient_aerosol_in_air / (ug/m3) \
(site_id: 5; time: 25)
    >>> for md in md_list:
    ...    print(md) # doctest: +ELLIPSIS
    <class '....PPData'> - Subclass of ADAQData - Contains:
    sites_cube_list:
    0: mass_concentration_of_pm2p5_dry_aerosol_in_air / (ug/m3) \
(site_id: 5; time: 25)
    1: mass_concentration_of_ozone_in_air / (ug/m3) \
(site_id: 5; time: 25)
    gridded_cube_list:
    0: mass_concentration_of_pm2p5_dry_aerosol_in_air / (ug/m3) \
(time: 25; grid_latitude: 182; grid_longitude: 146)
    1: mass_concentration_of_ozone_in_air / (ug/m3) \
(time: 25; grid_latitude: 182; grid_longitude: 146)
    trajectory_cube_list:
    < No cubes >
    <class '....PPData'> - Subclass of ADAQData - Contains:
    sites_cube_list:
    0: mass_concentration_of_pm2p5_dry_aerosol_in_air / (ug/m3) \
(site_id: 5; time: 25)
    1: mass_concentration_of_ozone_in_air / (ug/m3) (site_id: 5; time: 25)
    gridded_cube_list:
    0: mass_concentration_of_pm2p5_dry_aerosol_in_air / (ug/m3) \
(time: 25; grid_latitude: 182; grid_longitude: 146)
    1: mass_concentration_of_ozone_in_air / (ug/m3) \
(time: 25; grid_latitude: 182; grid_longitude: 146)
    trajectory_cube_list:
    < No cubes >

    >>> sc = md_list[0].extract(short_name='O3', singlecube=True)
    >>> print(sc) # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    mass_concentration_of_ozone_in_air / (ug/m3) (site_id: 5; time: 25)
         Dimension coordinates:
              site_id                                    x        -
              time                                       -        x
         Auxiliary coordinates:
              abbrev                                     x        -
              grid_latitude                              x        -
              grid_longitude                             x        -
              latitude                                   x        -
              longitude                                  x        -
              site_altitude                              x        -
              site_name                                  x        -
              site_type                                  x        -
              surface_altitude                           x        -
              forecast_period                            -        x
         Scalar coordinates:
              forecast_day: 1.0 Days
              level_height: 20.000... m, bound=(0.0, 49.998...) m
              model_level_number: 1
              sigma: 0.997..., bound=(1.0, 0.994...)
         Attributes:
              STASH: m01s34i001
              label: aqum_oper
              short_name: O3
              source: Data from Met Office Unified Model
         Cell methods:
              mean: time (1 hour)
    >>> print('{:.3f} {:.3f} {:.3f}'.format(sc.data.max(), sc.data.min(), sc.data.mean()))
    89.025 0.200 50.612
    >>> gc = md_list[0].extract(short_name='PM2p5', singlecube=True,
    ... gridded=True)
    >>> print(gc) # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    mass_concentration_of_pm2p5_dry_aerosol_in_air / (ug/m3) \
(time: 25; grid_latitude: 182; grid_longitude: 146)
         Dimension coordinates:
              time                x                  -                    -
              grid_latitude       -                  x                    -
              grid_longitude      -                  -                    x
         Auxiliary coordinates:
              forecast_period     x                  -                    -
              surface_altitude    -                  x                    x
         Derived coordinates:
              altitude            -                  x                    x
         Scalar coordinates:
              forecast_day: 1.0 Days
              level_height: 20.000... m, bound=(0.0, 49.998...) m
              model_level_number: 1
              sigma: 0.997..., bound=(1.0, 0.994...)
         Attributes:
              STASH: m01s17i221
              label: aqum_oper
              short_name: PM2p5
              source: Data from Met Office Unified Model
         Cell methods:
              mean: time (1 hour)
    >>> print('{:.3f} {:.3f} {:.3f}'.format(gc.data.max(), gc.data.min(), gc.data.mean()))
    197.086 0.449 25.949

    Create example stats cubes:

    >>> ini_dict, stats_cubes = create_example_data(exampletype='stats',
    ... save=False) # doctest: +ELLIPSIS
    Reading inifile .../example_data_10days.ini
    Number of sites:  5
    Saved to  .../stats_cube_list/stats_Obs-AURN_Mod-aqum_oper.nc
    ...
    Saved to  .../stats_cube_list/stats_Obs-AURN_Mod-aqum_oper.nc

    Check cubes:

    >>> print(stats_cubes)
    [[<iris 'Cube' of mass_concentration_of_ozone_in_air / (1) \
(istatistic: 34; time: 10)>,
    <iris 'Cube' of mass_concentration_of_pm2p5_ambient_aerosol_in_air / (1) \
(istatistic: 34; time: 10)>]]

    The first item in this list is an iris.cube.CubeList() corresponding to
    the first model in ini_dict['models_list'] which is aqum_oper, the second
    item therefore is a list of stats cubes corresponding to aqum_casestudy.

    >>> o3_oper = stats_cubes[0].extract(
    ... iris.AttributeConstraint(short_name='O3'), strict=True)
    >>> print(o3_oper)  # doctest: +ELLIPSIS
    mass_concentration_of_ozone_in_air / (1) (istatistic: 34; time: 10)
         Dimension coordinates:
              istatistic                                x         -
              time                                      -         x
         Auxiliary coordinates:
              statistic                                 x         -
              statistic_long_name                       x         -
              statistic_units                           x         -
         Attributes:
              ...
              label: aqum_oper
              obs: AURN
              short_name: O3

    Check the statistics included:

    >>> np.set_printoptions(linewidth=74)
    >>> print(o3_oper.coord('statistic').points)
    ['mdi' 'nsites' 'npts' 'correlation' 'bias' 'nmb' 'mnmb' 'mge' 'nmge'
     'fge' 'rmse' 'fac2' 'ioa' 'threshold' 'orss' 'odds_ratio' 'hitrate'
     'falsealarmrate' 'falsealarmratio' 'o>=t_m>=t' 'o<t_m>=t' 'o>=t_m<t'
     'o<t_m<t' 'maxobs' 'maxmod' 'meanobs' 'meanmod' 'sdobs' 'sdmod'
     'perc_correct' 'perc_over' 'perc_under' 'pc95obs' 'pc95mod']
    >>> np.set_printoptions(linewidth=75)
    """


    if exampletype == 'stats':
        #Use pre-created example data as a starting point.
        ini_dict, sites_data, od, md_list = get_exampledata(
            exampletype='10days')
        if save:
            ini_dict['plot_dir'] = './stats_cube_list'
        else:
            ini_dict['plot_dir'] = config.TEST_DIR + '/stats_cube_list'
        #Calculate statistics for each day and save to combined netcdf file.
        dates = [ini_dict['start_datetime'].date() + datetime.timedelta(days=n)
                 for n in range(ini_dict['range_days'])]
        for idate, date in enumerate(dates):
            #Extract observations object for this date only
            od_date = od.extract(time=lambda t: t.point.date() == date)
            #Generate list of model objects, with each object having a single
            #date of model data in.
            md_date_list = []
            for md in md_list:
                md_date = md.extract(time=lambda t: t.point.date() == date)
                md_date_list.append(md_date)
            #Save statistics
            if idate == 0:
                calc_stats(ini_dict, od_date, md_date_list, nc_append=False)
            else:
                calc_stats(ini_dict, od_date, md_date_list, nc_append=True)

        #Load data back in
        stats_cubes = []
        for model in ini_dict['models_list']:
            stats_cubes_model = iris.load(
                ini_dict['plot_dir'] + '/stats_Obs-AURN_Mod-' + model + '.nc')
            #Sort cubes alphabetically to ensure consistent doctests
            stats_cubes_model = iris.cube.CubeList(sorted(list(
                stats_cubes_model), key=lambda cube: cube.name()))
            stats_cubes.append(stats_cubes_model)

        return ini_dict, stats_cubes

    else:

        #Read ini file
        inifilename = 'adaqcode/example_data_' +exampletype + '.ini'
        ini_dict = inifile.get_inidict(defaultfilename=inifilename)
        #Get site data (if reading from file - otherwise None returned)
        sites_data = sites_info.get_siteinfo(ini_dict)

        #Read in observations and model data
        od = get_obs(ini_dict, sites_data)

        #Read in model data
        md_list = get_models(ini_dict, sites_data)

        #Unit conversion
        if exampletype != '3d':
            od, md_list = unit_conversion(od, md_list, chem_units='ug/m3',
                                          aerosol_units='ug/m3')

        if save:

            #Setup output directories
            sites_dir = "./sites_cube_list/"
            if not os.path.isdir(sites_dir):
                os.makedirs(sites_dir)
            gridded_dir = "./gridded_cube_list/"
            if not os.path.isdir(gridded_dir):
                os.makedirs(gridded_dir)

            #Output sites and gridded data
            suffix = '_' + exampletype + '.nc'
            od.save_ts(sites_dir + ini_dict['obs_fmt'] + suffix)
            for label, md in zip(ini_dict['models_list'], md_list):
                md.save_ts(sites_dir + label + suffix)
                print('Saved to '+sites_dir + label + suffix)
                md.save_gridded(gridded_dir + label + suffix)
                print('Saved to '+gridded_dir + label + suffix)

        return ini_dict, sites_data, od, md_list


def get_exampledata(sites_cube_list=True, gridded_cube_list=False,
                    exampletype="short"):
    """
    Get some example data to enable easier testing without reading in lots of
    data files.
    Note only two models are read in so this may not be consistent with
    ini_dict.

    :param sites_cube_list: Logical. If True, loads site cube data.
    :param gridded_cube_list: Logical. If True, loads gridded cube data.
                              Note by default this is False.
    :param exampletype: Type of example data to load
                        Currently set up:

                        * **short** - short 1 day example, only has two species
                        * **full** - fully 5 day example, with 5 species
                        * **3d** - 3d data, 3 hours of data with one species
                          (ozone) plus pressure and temperature
                        * **stats** - 5 days of 2 species saved as statistics
                          cubes
                        * **10days** - longer 10 day period, with only two
                          species, spread across two months

    :return:  (ini_dict, sites_data, od, md_list)
              Although if exampletype='stats', (ini_dict, stats_cubes) are
              returned instead.

     * **ini_dict** - Dictionary of a :class:`inifile` object
     * **sites_data** - numpy ndarray containing site information data
       from a :class:`sites_info.SitesInfo` object
     * **od** - Observation Data -  in :class:`adaq_data` format
     * **md_list** - Model Data list - each item in :class:`adaq_data` format

    .. Note:: Currently this only returns AQUM based data. We may want to
              add some NAME data to this (possibly determined by keyword
              'exampletype').

    >>> ini_dict, sites_data, od, md_list = get_exampledata()
    ... #  doctest: +ELLIPSIS
    Reading inifile .../adaqcode/example_data_1days.ini
    Number of sites:  5

    >>> ini_dict, stats_cubes = get_exampledata(exampletype='stats')
    ... #  doctest: +ELLIPSIS
    Reading inifile .../example_data_10days.ini
    >>> print(stats_cubes[0])
    0: mass_concentration_of_ozone_in_air / (1) \
(istatistic: 32; time: 10)
    1: mass_concentration_of_pm2p5_ambient_aerosol_in_air / (1) \
(istatistic: 32; time: 10)
    >>> np.set_printoptions(linewidth=74)
    >>> print(stats_cubes[0][0].coord('statistic').points)
    ['mdi' 'nsites' 'npts' 'correlation' 'bias' 'nmb' 'mnmb' 'mge' 'nmge'
     'fge' 'rmse' 'fac2' 'ioa' 'threshold' 'orss' 'odds_ratio' 'hitrate'
     'falsealarmrate' 'falsealarmratio' 'o>=t_m>=t' 'o<t_m>=t' 'o>=t_m<t'
     'o<t_m<t' 'maxobs' 'maxmod' 'meanobs' 'meanmod' 'sdobs' 'sdmod'
     'perc_correct' 'perc_over' 'perc_under']
    >>> np.set_printoptions(linewidth=74)
    """

    if exampletype == "short":
        inifilename = 'adaqcode/example_data_1days.ini'
    elif exampletype == "full":
        inifilename = 'adaqcode/example_data_5days.ini'
    elif exampletype == '3d':
        inifilename = 'adaqcode/example_data_3d.ini'
    elif exampletype == '10days':
        inifilename = 'adaqcode/example_data_10days.ini'
    elif exampletype == 'stats':
        inifilename = 'adaqcode/example_data_10days.ini'
    else:
        raise ValueError("Unknown exampletype:" + exampletype)

    ini_dict = inifile.get_inidict(defaultfilename=inifilename)

    if exampletype == 'stats':
        samplepath = config.SAMPLE_DATADIR + 'stats_cube_list/'
        stats_cubes = []
        for model in ini_dict['models_list']:
            stats_cubes_model = iris.load(
                samplepath + '/stats_Obs-AURN_Mod-' + model + '.nc')
            #Sort cubes alphabetically to ensure consistent doctests
            stats_cubes_model = iris.cube.CubeList(sorted(list(
                stats_cubes_model), key=lambda cube: cube.name()))
            stats_cubes.append(stats_cubes_model)
        return ini_dict, stats_cubes

    sites_data = None
    od = adaq_data.ADAQData()
    md1 = adaq_data.ADAQData()
    md2 = adaq_data.ADAQData()

    if sites_cube_list:
        sites_data = sites_info.get_siteinfo(ini_dict)
        samplepath = config.SAMPLE_DATADIR+'sites_cube_list/'
        if exampletype == "short":
            scl = od.load_ts(samplepath+'aurn_1days.nc')
            scl = md1.load_ts(samplepath+'aqum_oper_1days.nc')
            scl = md2.load_ts(samplepath+'aqum_casestudy_1days.nc')
        elif exampletype == "full":
            scl = od.load_ts(samplepath+'aurn_5days.nc')
            scl = md1.load_ts(samplepath+'aqum_oper_5days.nc')
            scl = md2.load_ts(samplepath+'aqum_casestudy_5days.nc')
        elif exampletype == '10days':
            scl = od.load_ts(samplepath+'aurn_10days.nc')
            scl = md1.load_ts(samplepath+'aqum_oper_10days.nc')
        elif exampletype == '3d':
            scl = md1.load_ts(samplepath+'aqum_oper_3d.nc')
    if gridded_cube_list:
        samplepath = config.SAMPLE_DATADIR+'gridded_cube_list/'
        if exampletype == "short":
            gcl = md1.load_gridded(samplepath+'aqum_oper_1days.nc')
            gcl = md2.load_gridded(samplepath+'aqum_casestudy_1days.nc')
        elif exampletype == "full":
            gcl = md1.load_gridded(samplepath+'aqum_oper_5days.nc')
            gcl = md2.load_gridded(samplepath+'aqum_casestudy_5days.nc')
        elif exampletype == "10days":
            gcl = md1.load_gridded(samplepath+'aqum_oper_10days.nc')
        elif exampletype == '3d':
            gcl = md1.load_gridded(samplepath+'aqum_oper_3d.nc')

    if exampletype == '3d' or exampletype == '10days':
        md_list = [md1]
    else:
        md_list = [md1, md2]

    return ini_dict, sites_data, od, md_list


def get_obs(ini_dict, sites_data):
    """
    Read in observations in :class:`adaq_data.ADAQData` format,
    for sites defined in sites_data and start/end dates and short_names
    as defined in ini_dict.

    :param ini_dict: Dictionary of a :class:`inifile` object
                     Should contain:

                      * 'obs_fmt' - currently only 'aurn','camsaqobs'
                        and 'pollen' allowed
                      * 'obs_dir' - directory contain obs files
                      * 'short_name_list'
                      * 'start_datetime'
                      * 'end_datetime'

    :param sites_data: numpy ndarray containing site information data
                       from a :class:`sites_info.SitesInfo` object

    :return: Observation data in the form of an :class:`adaq_data.ADAQData`
             object

    .. note:: Currently, only a single observation type is returned.

    .. note:: If 'obs_dir' is not in ini_dict (or if obs_fmt is not a
              currently valid format type),
              then no observations are read in and the observations are
              returned as a single ADAQData object.

    >>> ini_dict = inifile.get_inidict(
    ... defaultfilename='adaqcode/example_data_1days.ini') # doctest: +ELLIPSIS
    Reading inifile .../adaqcode/example_data_1days.ini
    >>> sites_data = sites_info.get_siteinfo(ini_dict)
    Number of sites:  5
    >>> od = get_obs(ini_dict, sites_data) # doctest: +ELLIPSIS
    Creating observation data at ...
    Reading obs data files
    Found obs for  .../AURN_obs/ABD_20140101_20140818.txt
    Found obs for  .../AURN_obs/ACTH_20140101_20140818.txt
    Found obs for  .../AURN_obs/AH_20140101_20140818.txt
    Found obs for  .../AURN_obs/HAR_20140101_20140818.txt
    Found obs for  .../AURN_obs/YW_20140101_20140818.txt
    Creating obs cubes
    >>> print(od) # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    <class '...AURNData'> - Subclass of ADAQData - Contains:
    sites_cube_list:
    0: mass_concentration_of_ozone_in_air / (ug/m3) \
(site_id: 5; time: 25)
    1: mass_concentration_of_pm2p5_ambient_aerosol_in_air / (ug/m3) \
(site_id: 5; time: 25)
    gridded_cube_list:
    < No cubes >
    trajectory_cube_list:
    < No cubes >
    >>> cube = od.extract(short_name='O3',singlecube=True)
    >>> print(cube.data.max(), cube.data.min())
    96.0 2.0
    """


    #--------------------------------------------------
    #Get observations

    #Set up list of obs directories - use from obs_dir_list if possible,
    #otherwise from obs_dir
    obs_dir_list = ini_dict.get('obs_dir_list', None)
    if not obs_dir_list:
        if 'obs_dir' in ini_dict:
            obs_dir_list = [ini_dict['obs_dir']]

        else:
            #No observation directory set - can't use obs
            #Set observations as blank ADAQData object
            #this is so can still be compared to model data
            #(or to only plot model data)
            od = adaq_data.ADAQData()
            return od

    print("Creating observation data at ", datetime.datetime.now())

    obs_fmt = ini_dict.get('obs_fmt', None)
    if obs_fmt == 'aurn':
        od = aurn_data.AURNData()
        od.readdata(obsdir_list=obs_dir_list,
                    short_name_list=ini_dict['short_name_list'],
                    sites=sites_data,
                    start_datetime=ini_dict['start_datetime'],
                    end_datetime=ini_dict['end_datetime'],
                    obsdir_fixed_csv=ini_dict.get('obsdir_fixed_csv', None))
    elif obs_fmt == 'camsaqobs':
        #Only accepts a single directory (could contain wildcards)
        obs_dir = obs_dir_list[0]
        od = camsaqobs_data.CAMSAQObsData()
        od.readdata(filenames=obs_dir,
                    short_name_list=ini_dict['short_name_list'],
                    sites_data=sites_data,
                    site_types=ini_dict.get('site_types_list', None),
                    start_datetime=ini_dict['start_datetime'],
                    end_datetime=ini_dict['end_datetime'])
    elif obs_fmt == 'pollen':
        #Only accepts a single directory (could contain wildcards)
        od = pollenobs_data.PollenObsData()
        od.readdata(file_pattern=obs_dir_list[0],
                    short_name_list=ini_dict['short_name_list'],
                    sites_data=sites_data,
                    start_datetime=ini_dict['start_datetime'],
                    end_datetime=ini_dict['end_datetime'])

    else:
        warnings.warn('Unknown observations fmt: ' + obs_fmt)
        #No observation fmt set - can't use obs
        #Set observations as blank ADAQData object
        #this is so can still be compared to model data
        #(or to only plot model data)
        #also set sites_cube_list and gridded_cube_list to
        # a cubelist (rather than None) for the same reason.
        od = adaq_data.ADAQData()

    return od

def get_models(ini_dict, sites_data=None, delete_gridded=False,
               surface_only=False):
    """
    Read in models in :class:`adaq_data.ADAQData` format,
    for sites defined in sites_data and start/end dates and short_names
    as defined in ini_dict.

    :param ini_dict: Dictionary of a :class:`inifile` object.
                     Should contain:

                      * 'models_list'
                      * 'models_dir_list'
                      * 'models_fmt_list' (currently supported: 'macc_ens',
                        'name','nimrod','pp')
                      * 'short_name_list'
                      * 'start_datetime'
                      * 'end_datetime'
                      * 'forecast_day' (depending on data type)
                      * 'field_attribute_dict' (for NAME format only)

    :param sites_data: numpy ndarray containing site information data
                       from a :class:`sites_info.SitesInfo` object
    :param delete_gridded: Set to True to delete the gridded cubes as soon
                         as sites cubes are created to save memory.
    :param surface_only: Set to True to return surface data only.
                         Set to False to get all model levels.
                         Note: not set up for NAMEData.
    :return: List of model data objects in the form of
             :class:`adaq_data.ADAQData` objects.

    .. note:: If 'models_dir_list' is not in ini_dict,
              then no models are read in and md_list is returned as
              a single ADAQData object in a list.

    >>> ini_dict = inifile.get_inidict(
    ... defaultfilename='adaqcode/example_data_1days.ini') # doctest: +ELLIPSIS
    Reading inifile .../adaqcode/example_data_1days.ini
    >>> sites_data = sites_info.get_siteinfo(ini_dict)
    Number of sites:  5
    >>> md_list = get_models(ini_dict, sites_data) # doctest: +ELLIPSIS
    Getting model data for  aqum_oper  at  ...
    Getting model data for  aqum_casestudy  at  ...
    >>> print(md_list) # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    [<...pp_data.PPData object at ...>,
    <...pp_data.PPData object at ...]
    >>> cube = md_list[0].extract(short_name='O3',singlecube=True)
    >>> print('{:.3e}'.format(cube.data.max()))
    7.390e-08
    """

    verbose = 1 #Set to 2 for more detailed print output for debugging

    if 'models_dir_list' not in ini_dict:
        #No models listed
        #Set up blank md_list
        #this is so can still be compared to obs data
        #(or to only plot obs data)
        md = adaq_data.ADAQData()
        md_list = [md]
        return md_list

    md_list = []
    for imodel, model in enumerate(ini_dict['models_list']):
        datadir = ini_dict['models_dir_list'][imodel]
        fmt = ini_dict['models_fmt_list'][imodel]

        if model is None:
            continue

        if verbose > 0:
            print("Getting model data for ", model, \
                  " at ", datetime.datetime.now())

        if fmt == 'ecgrib':

            md = ecgrib_data.ECGribData()
            md.get_filenames(datadir,
                             start_datetime=ini_dict['start_datetime'],
                             end_datetime=ini_dict['end_datetime'],
                             forecast_day=ini_dict['forecast_day'],
                             short_name_list=ini_dict['short_name_list'])
            md.readdata(label=model,
                        surface_only=surface_only)

        elif fmt == 'maccens':

            md = maccens_data.MaccEnsData()
            md.get_filenames(datadir,
                             start_datetime=ini_dict['start_datetime'],
                             end_datetime=ini_dict['end_datetime'],
                             forecast_day=ini_dict['forecast_day'])
            if surface_only:
                level = 0
            else:
                level = None
            md.readdata(short_name_list=ini_dict['short_name_list'],
                        start_datetime=ini_dict['start_datetime'],
                        end_datetime=ini_dict['end_datetime'],
                        label=model,
                        level=level)

        elif fmt == 'name':

            z_levels = ini_dict.get('z_level_list', None)
            z_leveltype = ini_dict.get('z_leveltype', None)
            back = ini_dict.get('back', None)
            short_name_list = ini_dict.get('short_name_list', None)
            field_attributes = ini_dict.get('field_attribute_dict', None)
            start_datetime = ini_dict.get('start_datetime', None)
            end_datetime = ini_dict.get('end_datetime', None)

            md = name_data.NAMEData()
            md.get_filenames(datadir,
                             start_datetime=start_datetime,
                             end_datetime=end_datetime)
            md.readdata(field_attributes=field_attributes,
                        z_levels=z_levels,
                        z_leveltype=z_leveltype,
                        back=back,
                        label=model,
                        short_name_list=short_name_list)

        elif fmt == 'ens_name':

            z_levels = ini_dict.get('z_level_list', None)
            z_leveltype = ini_dict.get('z_leveltype', None)
            back = ini_dict.get('back', None)
            short_name_list = ini_dict.get('short_name_list', None)
            field_attributes = ini_dict.get('field_attribute_dict', None)
            start_datetime = ini_dict.get('start_datetime', None)
            end_datetime = ini_dict.get('end_datetime', None)
            nmembers_ini = int(ini_dict.get('num_members', 18))

            md = name_data.NAMEData()
            md.get_filenames(datadir,
                             start_datetime=start_datetime,
                             end_datetime=end_datetime)

            directory = os.path.dirname(datadir)
            workdir = os.path.dirname(directory)
            member_dir = os.path.basename(directory).replace('*', '{0}')
            nmembers = len(glob.glob(directory))

            if int(nmembers) != int(nmembers_ini):
                raise ValueError('Expecting ' + str(nmembers_ini) +
                                 ' ensemble members but found: ' +
                                 str(nmembers))

            ensemble_dict = {}
            for imem in range(nmembers):
                fullpath_member_dir = os.path.join(workdir, member_dir.format(imem))
                member_label = member_dir.format(imem)
                ensemble_dict[fullpath_member_dir + '/'] = [member_label, imem]

            md.readdata(ensemble_dict=ensemble_dict,
                        field_attributes=field_attributes,
                        z_levels=z_levels,
                        z_leveltype=z_leveltype,
                        back=back,
                        label=model,
                        short_name_list=short_name_list)

        elif fmt == 'nimrod':

            #Get data_type - default to sppo if not set up.
            data_type = None
            data_type_list = ini_dict.get('nimrod_data_type_list', None)
            if data_type_list is not None:
                data_type = data_type_list[imodel]
            if data_type is None:
                data_type = 'sppo'

            md = nimrod_data.NimrodData()
            md.get_filenames(datadir,
                             start_datetime=ini_dict['start_datetime'],
                             end_datetime=ini_dict['end_datetime'],
                             forecast_day=ini_dict['forecast_day'])
            md.readdata(data_type=data_type,
                        short_name_list=ini_dict['short_name_list'],
                        start_datetime=ini_dict['start_datetime'],
                        end_datetime=ini_dict['end_datetime'],
                        label=model)

        elif fmt == 'pp':

            #--------------------------------------------------
            # Get a list of pp filenames between the specified times

            ppfiles = pp_data.AQUMppFiles()
            filenames = ppfiles.get_filenames(
                datadir,
                start_datetime=ini_dict['start_datetime'],
                end_datetime=ini_dict['end_datetime'],
                forecast_day=ini_dict['forecast_day'])

            #--------------------------------------------------
            # Read model data
            md = pp_data.PPData()
            md.readdata(filenames=filenames,
                        short_name_list=ini_dict['short_name_list'],
                        start_datetime=ini_dict['start_datetime'],
                        end_datetime=ini_dict['end_datetime'],
                        forecast_day=ini_dict['forecast_day'],
                        label=model,
                        surface_only=surface_only)

        elif fmt == 'traj':
            md = trajectory_data.TrajData()
            md.readdata(datadir)

        else:
            raise ValueError('Unknown model format ' + fmt)

        #--------------------------------------------------
        # Extract site data from model

        if sites_data is not None:
            if not md.gridded_cube_list:
                warnings.warn('No gridded data read in for ' + model)
            else:
                if verbose > 1:
                    print(("Creating model site cubes at ",
                           datetime.datetime.now()))
                md.extract_sites(sites_data=sites_data)

        if delete_gridded:
            md.gridded_cube_list = iris.cube.CubeList()

        md_list.append(md)

        if verbose > 1:
            print("Finished getting model data at ", datetime.datetime.now())

    return md_list


def add_missing_times(od, md_list):
    """
    For plotting time-series, or calculation of statistics which require
    knowledge of previous hours, the time-coordinate needs to have all
    possible times, even if the data values are set to nan.
    This may occur if for example, a day of datafiles don't exist.
    This routine adds the missing times to sites_cube_list and gridded_cube_list
    for observations and model data objects.

    :param od: observation data in the form of an
               :class:`adaq_data.ADAQData` object
    :param md_list: list of model data objects in the form of
                    :class:`adaq_data.ADAQData` objects
    :return: od, md_list in the same formats as input, but with the missing
             times added to the time-coordinate and the data values for these
             points set to nan.

    >>> ini_dict, sites_data, od, md_list = get_exampledata(
    ... exampletype="full") # doctest: +ELLIPSIS
    Reading inifile .../adaqcode/example_data_5days.ini
    Number of sites:  5

    Remove the centre times of a cube for example purposes:

    >>> cube1 = md_list[0].sites_cube_list[0][:,:48]
    >>> cube2 = md_list[0].sites_cube_list[0][:,72:]
    >>> cube = iris.cube.CubeList([cube1, cube2]).concatenate_cube()
    >>> md_list[0].sites_cube_list[0] = cube
    >>> print(md_list[0].sites_cube_list[0].summary(shorten=True))
    mass_concentration_of_nitrogen_dioxide_in_air / (ug/m3) \
(site_id: 5; time: 97)
    >>> od, md_list = add_missing_times(od, md_list)
    Adding missing times
    >>> print(md_list[0].sites_cube_list[0].summary(shorten=True))
    mass_concentration_of_nitrogen_dioxide_in_air / (ug/m3) \
(site_id: 5; time: 121)
    """

    print('Adding missing times')

    if od.sites_cube_list:
        for i, cube in enumerate(od.sites_cube_list):
            if len(cube.coord('time').points) > 1:
                od.sites_cube_list[i] = \
                    cube_time.cube_add_missing_times(cube)
    if od.gridded_cube_list:
        for i, cube in enumerate(od.gridded_cube_list):
            if len(cube.coord('time').points) > 1:
                od.gridded_cube_list[i] = \
                    cube_time.cube_add_missing_times(cube)
    for md in md_list:
        if md.sites_cube_list:
            for i, cube in enumerate(md.sites_cube_list):
                if len(cube.coord('time').points) > 1:
                    md.sites_cube_list[i] = \
                        cube_time.cube_add_missing_times(cube)
        if md.gridded_cube_list:
            for i, cube in enumerate(md.gridded_cube_list):
                if len(cube.coord('time').points) > 1:
                    md.gridded_cube_list[i] = \
                        cube_time.cube_add_missing_times(cube)

    return od, md_list


def match_data_for_stats(ini_dict, od, md_list):
    """
    For fair calculation of statistics, observations and model data
    statistics should only be calculated where data exist for all datasets.
    This routine therefore ensures the times and sites from observations
    and all model data match for each short-name. Furthermore, whereever
    there is missing (NaN) data for any observation or model point, NaNs are
    set across all observations and model datasets at this point.
    However if any model or observations are entirely missing, this is ignored
    and data is only matched amongst available model/observations

    :param ini_dict: Dictionary of a :class:`inifile` object. Should contain:

                       * 'short_name_list'
                       * 'strict_statistics' (optional) - if False, returns
                         original data and not matched.

    :param od: observation data in the form of an
               :class:`adaq_data.ADAQData` object.
    :param md_list: list of model data objects in the form of
                    :class:`adaq_data.ADAQData` objects/
    :return: od, md_list in the same formats as input, but with extra
             times/sites removed and missing data matched across all
             observations/model datasets.

    .. note:: This only acts on the sites_cube_list (not gridded_cube_list).

    First get some data, and simplify example by only using one set of
    model data:

    >>> ini_dict, sites_data, od, md_list = get_exampledata(exampletype="full")
    ... # doctest: +ELLIPSIS
    Reading inifile .../adaqcode/example_data_5days.ini
    Number of sites:  5
    >>> md_list = [ md_list[0] ]

    Extract a single site and short name, and for the purposes of this example,
    remove some times from the start of the model:

    >>> obs_no2_yw = od.extract(short_name='NO2',abbrev='YW',singlecube=True)
    >>> mod_no2_yw = md_list[0].extract(short_name='NO2',abbrev='YW',
    ... singlecube=True)

    Note that there is some missing data in the observations which are present
    in the models.

    >>> np.set_printoptions(formatter={'float':lambda x: '{:5.2f}'.format(x)})
    >>> print(obs_no2_yw.data[45:48], mod_no2_yw.data[45:48])
    [ 3.20   nan  0.10] [13.04 11.79 13.88]

    Run the routine to correct these:

    >>> od_strict, md_list_strict = match_data_for_stats(ini_dict, od, md_list)

    And check again the data:

    >>> obs_strict_no2_yw = od_strict.extract(short_name='NO2',
    ... abbrev='YW',singlecube=True)
    >>> mod_strict_no2_yw = md_list_strict[0].extract(short_name='NO2',
    ... abbrev='YW',singlecube=True)
    >>> print(obs_strict_no2_yw.data[45:48], mod_strict_no2_yw.data[45:48])
    [ 3.20   nan  0.10] [13.04   nan 13.88]
    >>> np.set_printoptions()

    These both now have the same start times and
    model data is missing where the obs are.

    Note that if one of the datasets has extra times with valid data
    at the beginnning of end which the other dataset does not, the data
    at these extra times will be removed this routine. For more information
    on this, see :func:`timeseries_stats.match_cubes`

    Note that if there are no observations for a particular short_name but
    there are multiple models then the models are only matched amongst
    themselves and not all set to missing data.

    """

    #Check value of strict_statistics in ini_dict
    #Assume default of True if not set.
    if not ini_dict.get('strict_statistics', True):
        #Don't match - return a copy of all original data
        #Use copy to ensure not pointing at input versions.
        return copy.deepcopy(od), copy.deepcopy(md_list)

    else:
        #Loop over cubes in cubelists. Only need to do sites cubes??

        if not od.sites_cube_list:
            #Can't match - return a copy of all original data
            #Use copy to ensure not pointing at input versions.
            warnings.warn("Can not match for statistics: "+\
                          "no observations sites_cube_list")
            return copy.deepcopy(od), copy.deepcopy(md_list)

        #Make a new version of the observation and model data objects,
        #set the sites_cube_list as a blank cubelist
        od_strict = adaq_data.ADAQData()
        od_strict.sites_cube_list = iris.cube.CubeList()
        md_list_strict = [adaq_data.ADAQData() for md in md_list]
        for md in md_list_strict:
            md.sites_cube_list = iris.cube.CubeList()

        #Loop over all short names, extracting singlecubes,
        # then matching data/sites etc
        for short_name in ini_dict['short_name_list']:
            cubes = []

            #Extract this short_name cube from observations
            obs_cube_strict_tmp = od.extract(short_name=short_name,
                                             singlecube=True)
            if obs_cube_strict_tmp is not None:
                obs_cube_strict = obs_cube_strict_tmp.copy()
                cubes.append(obs_cube_strict)
                #Take note that matched_cubes will contain obs
                use_obs = True
            else:
                #Obs don't contain any cubes for this short_name
                #Can only match between different models
                warnings.warn("No observations cube found for "+\
                              short_name)
                #Take note that matched_cubes will not contain obs
                use_obs = False

            use_model = []
            for md in md_list:
                #Extract this short_name cube from model data
                mod_cube_strict_tmp = md.extract(short_name=short_name,
                                                 singlecube=True)
                if mod_cube_strict_tmp is not None:
                    mod_cube_strict = mod_cube_strict_tmp.copy()
                    cubes.append(mod_cube_strict)
                    use_model.append(True)
                else:
                    warnings.warn("No model cube found for "+\
                                  short_name)
                    use_model.append(False)

            #Match cubes against times/sites etc
            matched_cubes = cube_statistics.match_cubes(cubes)

            #Add into sites_cube_list for strict observations...
            i_cube = 0
            if use_obs:
                #Can only add if obs available for this short_name
                od_strict.sites_cube_list.append(matched_cubes[i_cube])
                i_cube += 1
            for imd, md in enumerate(md_list_strict):
                #...and for strict model data
                if use_model[imd]:
                    #Can only add if md available for this short_name
                    md.sites_cube_list.append(matched_cubes[i_cube])
                    i_cube += 1

    return od_strict, md_list_strict


def unit_conversion(od, md_list, chem_units=None, aerosol_units=None,
                    pollen_units=None, at_stp=True):
    """
    Sets up cubelists of ADAQData objects prior to unit and phenomenon
    conversion. This function takes as inputs an observation ADAQData object
    and a list of model ADAQData objects. For non-stp conversions these
    must contain pressure and temperature cubes.

    Currently only conversions to mass concentrations in ug/m3 and to
    volume mixing ratios in ppb are allowed.

    :param od: observation data in the form of an
               :class:`adaq_data.ADAQData` object.
    :param md_list: list of model data objects in the form of
                    :class:`adaq_data.ADAQData` objects.
    :param chem_units: the units that chemical species are converted to.
                       By default, if aerosol_units are used, then chem_units
                       are set to 'ug/m3'.
    :param aerosol_units: the units that aerosol species are converted to.
                          By default, if chem_units are used, then
                          aerosol_units are set to 'ug/m3'.
    :param pollen_units: the units that pollen species are converted to.
    :param at_stp: Set to True for stp assumption, otherwise set False and
                   supply the p and T required to carry out the conversion.
    :return: od, md_list: observation and model data list ADAQData objects
                          with chemical species in transformed units.


    This example is for 3D data unit conversion. However the routine can also
    be used for converting lists of observation and/or model surface data,
    either with or without the assumption of stp.

    First read in some 3D data:

    >>> ini_dict, sites_data, od, md_list = get_exampledata(
    ...     exampletype='3d', gridded_cube_list=True) # doctest: +ELLIPSIS
    Reading inifile .../example_data_3d.ini
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

    Check the current spread of data, at both the top and bottom levels to
    observe the effect temperature and pressure will have

    >>> ozone = md_list[0].extract(short_name='O3', singlecube=True, gridded=True)
    >>> print('Surface mean = {:.3e} {}'.format(ozone.data[0,0,:,:].mean(), ozone.units))
    Surface mean = 5.455e-08 kg kg-1
    >>> print('Top layer mean = {:.3e} {}'.format(ozone.data[0,-1,:,:].mean(), ozone.units))
    Top layer mean = 6.300e-06 kg kg-1

    Now convert units:

    >>> od, md_list = unit_conversion(od, md_list, chem_units='ug/m3',
    ... aerosol_units='ug/m3', at_stp=False)
    >>> print(md_list[0]) # doctest: +ELLIPSIS
    <class '....ADAQData'> - Subclass of ADAQData - Contains:
    sites_cube_list:
    0: air_pressure / (Pa)                 (site_id: 5; time: 2; \
model_level_number: 38)
    1: air_temperature / (K)               (site_id: 5; time: 2; \
model_level_number: 38)
    2: mass_concentration_of_ozone_in_air / (ug/m3) (site_id: 5; time: 2; \
model_level_number: 38)
    3: specific_humidity / (kg kg-1)       (site_id: 5; time: 2; \
model_level_number: 38)
    gridded_cube_list:
    0: air_pressure / (Pa)                 (time: 2; model_level_number: 38; \
grid_latitude: 182; grid_longitude: 146)
    1: air_temperature / (K)               (time: 2; model_level_number: 38; \
grid_latitude: 182; grid_longitude: 146)
    2: mass_concentration_of_ozone_in_air / (ug/m3) (time: 2; \
model_level_number: 38; grid_latitude: 182; grid_longitude: 146)
    3: specific_humidity / (kg kg-1)       (time: 2; model_level_number: 38; \
grid_latitude: 182; grid_longitude: 146)
    trajectory_cube_list:
    < No cubes >

    Confirm the new data has sensible values.

    >>> ozone = md_list[0].extract(short_name='O3', singlecube=True, gridded=True)
    >>> print('Surface mean = {:.3f} {}'.format(ozone.data[0,0,:,:].mean(), ozone.units))
    Surface mean = 65.953 ug/m3
    >>> print('Top layer mean = {:.3f} {}'.format(ozone.data[0,-1,:,:].min(), ozone.units))
    Top layer mean = 27.490 ug/m3


    Similarly we can convert to ppb.

    >>> ini_dict, sites_data, od, md_list = get_exampledata(
    ...     exampletype='short', gridded_cube_list=True) # doctest: +ELLIPSIS
    Reading inifile .../example_data_1days.ini
    Number of sites:  5
    >>> ozone = md_list[0].extract(short_name='O3', singlecube=True,
    ...     gridded=True)
    >>> print('Mean ozone = {:.3f} {}'.format(ozone.data.mean(), ozone.units))
    Mean ozone = 78.087 ug/m3
    >>> pm = md_list[0].extract(short_name='PM2p5', singlecube=True,
    ...     gridded=True)
    >>> print('Mean PM 2.5 = {:.3f} {}'.format(pm.data.mean(), pm.units))
    Mean PM 2.5 = 25.949 ug/m3

    Here we use chem_units as aerosols cannot be converted to ppb.

    >>> od, md_list = unit_conversion(od, md_list, chem_units='ppb',
    ... at_stp=True)
    >>> ozone = md_list[0].extract(short_name='O3', singlecube=True,
    ...     gridded=True)
    >>> print('Mean ozone = {:.3f} {}'.format(ozone.data.mean(), ozone.units))
    Mean ozone = 39.124 ppb
    >>> pm = md_list[0].extract(short_name='PM2p5', singlecube=True,
    ...     gridded=True)
    >>> print('Mean PM 2.5 = {:.3f} {}'.format(pm.data.mean(), pm.units))
    Mean PM 2.5 = 25.949 ug/m3

    """

    if od is not None:
        if od.sites_cube_list:
            if chem_units or aerosol_units:
                od.sites_cube_list = cube_chemistry.convert_cubelist_units(
                    od.sites_cube_list, chem_units=chem_units,
                    aerosol_units=aerosol_units, at_stp=at_stp)
            if pollen_units:
                od.sites_cube_list = pollenobs_data.convert_cubelist_units(
                    od.sites_cube_list, pollen_units=pollen_units)

        if od.gridded_cube_list:
            if chem_units or aerosol_units:
                od.gridded_cube_list = cube_chemistry.convert_cubelist_units(
                    od.gridded_cube_list, chem_units=chem_units,
                    aerosol_units=aerosol_units, at_stp=at_stp)
            if pollen_units:
                od.gridded_cube_list = pollenobs_data.convert_cubelist_units(
                    od.gridded_cube_list, pollen_units=pollen_units)


    for md in md_list:
        if md is not None:
            if md.sites_cube_list:
                if chem_units or aerosol_units:
                    md.sites_cube_list = cube_chemistry.convert_cubelist_units(
                        md.sites_cube_list, chem_units=chem_units,
                        aerosol_units=aerosol_units, at_stp=at_stp)
                if pollen_units:
                    md.sites_cube_list = pollenobs_data.convert_cubelist_units(
                        md.sites_cube_list, pollen_units=pollen_units)
            if md.gridded_cube_list:
                if chem_units or aerosol_units:
                    md.gridded_cube_list = cube_chemistry.convert_cubelist_units(
                        md.gridded_cube_list, chem_units=chem_units,
                        aerosol_units=aerosol_units, at_stp=at_stp)
                if pollen_units:
                    md.gridded_cube_list = pollenobs_data.convert_cubelist_units(
                        md.gridded_cube_list, pollen_units=pollen_units)


    return od, md_list

def scaling_factors(ini_dict, adaq_list):
    """
    Apply scaling factors to adaq data objects.
    Each scaling factor in the ini_dict['scaling_factors_list'] will be applied
    to :class:`adaq_data.ADAQData` objects in the order given in adaq_list.
    The same scaling factor will be applied to all cubes (all short_names, gridded,
    site and trajectory cubes).

    :param ini_dict: Dictionary of a :class:`inifile` object. Should contain
                     'scaling_factors_list'. If this is not available in the
                     ini_dict then no factors are applied. If an individual
                     scaling factor is None, then no factor is applied for this
                     particular adaq_data object in the list.
    :param adaq_list: List of :class:`adaq_data.ADAQData` objects.

    :returns: List of :class:`adaq_data.ADAQData` objects. Note the original
              input list will have also been modified.

    >>> ini_dict, sites_data, od, md_list = get_exampledata()
    ... # doctest: +ELLIPSIS
    Reading inifile .../example_data_1days.ini
    Number of sites:  5

    Set up scaling factors and list of :class:`adaq_data.ADAQData` objects:

    >>> ini_dict['scaling_factors_list'] = [None, 1.0, 10.]
    >>> adaq_list = [od] + md_list

    Check raw data:

    >>> print(np.nanmax(od.sites_cube_list[0].data),
    ... np.nanmax(md_list[0].sites_cube_list[0].data),
    ... np.nanmax(md_list[1].sites_cube_list[0].data))
    96.0 89.02489 92.46645

    >>> adaq_list = scaling_factors(ini_dict, adaq_list)
    Applying scaling factors

    Check scaled data - only the last one has been modified by a factor of 10,
    as requested:

    >>> print(np.nanmax(od.sites_cube_list[0].data),
    ... np.nanmax(md_list[0].sites_cube_list[0].data),
    ... np.nanmax(md_list[1].sites_cube_list[0].data))
    96.0 89.02489 924.66455


    """

    if 'scaling_factors_list' in ini_dict:
        print('Applying scaling factors')
        for scaling_factor, adaq_data in zip(ini_dict['scaling_factors_list'], adaq_list):
            if scaling_factor is not None:
                scaling_factor = float(scaling_factor)
                for cube in adaq_data.sites_cube_list:
                    cube.data *= scaling_factor
                for cube in adaq_data.gridded_cube_list:
                    cube.data *= scaling_factor
                for cube in adaq_data.trajectory_cube_list:
                    cube.data *= scaling_factor

    return adaq_list



def convert_to_daqi(ini_dict, od, md_list, od_stats=None, md_list_stats=None):
    """
    Convert all cubes to DAQI where possible (see :mod:`aq_indices`)
    Original cubes are left in place, with new cubes added to ends of cubelists.
    Total DAQI cube is also added to the cubelists.
    Short name list is modified to reflect these changes.

    :param ini_dict: Dictionary of a :class:`inifile` object. Should contain:

                       * 'short_name_list'

    :param od: observation data in the form of an
               :class:`adaq_data.ADAQData` object.
    :param md_list: list of model data objects in the form of
                    :class:`adaq_data.ADAQData` objects/
    :param od_stats: observation data in the form of
                     :class:`adaq_data.ADAQData` object, which has been matched
                     with times/sites where model data is available.
    :param md_list: list of model data objects in the form of
                    :class:`adaq_data.ADAQData` objects, which has been matched
                    with times/sites where model data is available.
    :return: ini_dict, od, md_list, od_stats, md_list_stats in the same formats
             as input, but with extra cubes added to sites_cube_list
             and gridded_cube_list to represent DAQI.

    >>> ini_dict, sites_data, od, md_list = get_exampledata(exampletype="full")
    ... # doctest: +ELLIPSIS
    Reading inifile .../adaqcode/example_data_5days.ini
    Number of sites:  5

    >>> print(od.sites_cube_list)
    0: mass_concentration_of_nitrogen_dioxide_in_air / \
(ug/m3) (site_id: 5; time: 121)
    1: mass_concentration_of_ozone_in_air / (ug/m3) \
(site_id: 5; time: 121)
    2: mass_concentration_of_pm10_ambient_aerosol_in_air / \
(ug/m3) (site_id: 5; time: 121)
    3: mass_concentration_of_pm2p5_ambient_aerosol_in_air / \
(ug/m3) (site_id: 5; time: 121)
    4: mass_concentration_of_sulphur_dioxide_in_air / \
(ug/m3) (site_id: 5; time: 121)

    >>> ini_dict, od, md_list, _, _ = convert_to_daqi(
    ... ini_dict, od, md_list)
    Converting to DAQI
    >>> print(od.sites_cube_list)
    0: mass_concentration_of_nitrogen_dioxide_in_air / \
(ug/m3) (site_id: 5; time: 121)
    1: mass_concentration_of_ozone_in_air / (ug/m3) \
(site_id: 5; time: 121)
    2: mass_concentration_of_pm10_ambient_aerosol_in_air / \
(ug/m3) (site_id: 5; time: 121)
    3: mass_concentration_of_pm2p5_ambient_aerosol_in_air / \
(ug/m3) (site_id: 5; time: 121)
    4: mass_concentration_of_sulphur_dioxide_in_air / (ug/m3) \
(site_id: 5; time: 121)
    5: daily_air_quality_index_of_nitrogen_dioxide / (1) \
(site_id: 5; time: 6)
    6: daily_air_quality_index_of_ozone / (1) \
(site_id: 5; time: 6)
    7: daily_air_quality_index_of_pm10_ambient_aerosol / (1) \
(site_id: 5; time: 6)
    8: daily_air_quality_index_of_pm2p5_ambient_aerosol / (1) \
(site_id: 5; time: 6)
    9: daily_air_quality_index_of_sulphur_dioxide / (1) \
(site_id: 5; time: 6)
    10: daily_air_quality_index / (1)       \
(site_id: 5; time: 6)

    """

    print('Converting to DAQI')

    #Modify short name list to add _DAQI where appropriate and
    #also add DAQI in
    short_name_list_with_daqi = []
    for short_name in ini_dict['short_name_list']:
        if short_name in aq_indices.DAQI_LEVELS:
            short_name_list_with_daqi.append(short_name + '_DAQI')
    short_name_list_with_daqi.append('DAQI')
    ini_dict['short_name_list'] += short_name_list_with_daqi

    if od.sites_cube_list:
        od.sites_cube_list = aq_indices.convert_all_daqi(
            od.sites_cube_list)

    if od.gridded_cube_list:
        od.gridded_cube_list = aq_indices.convert_all_daqi(
            od.gridded_cube_list)

    for md in md_list:

        if md.sites_cube_list:
            md.sites_cube_list = aq_indices.convert_all_daqi(
                md.sites_cube_list)

        if md.gridded_cube_list:
            md.gridded_cube_list = aq_indices.convert_all_daqi(
                md.gridded_cube_list)

    if od_stats is not None:
        if od_stats.sites_cube_list:
            od_stats.sites_cube_list = aq_indices.convert_all_daqi(
                od_stats.sites_cube_list)

        if od_stats.gridded_cube_list:
            od_stats.gridded_cube_list = aq_indices.convert_all_daqi(
                od_stats.gridded_cube_list)

    if md_list_stats is not None:

        for md in md_list_stats:

            if md.sites_cube_list:
                md.sites_cube_list = aq_indices.convert_all_daqi(
                    md.sites_cube_list)

            if md.gridded_cube_list:
                md.gridded_cube_list = aq_indices.convert_all_daqi(
                    md.gridded_cube_list)

    return ini_dict, od, md_list, od_stats, md_list_stats


def calc_stats(ini_dict, od, md_list, statsfile_prefix=None, thresholds=None,
               nc_append=False, stats_tcoord=None):
    '''
    Given an ObservationData object and list of ModelData objects,
    calculates statistics, comparing each model to the observations,
    outputting to statistics file if set.

    :param ini_dict: Dictionary of a :class:`inifile` object.
                     Could contain 'plot_dir' (optional)
                     - if not set, uses statsfile, or else doesn't write a file.

    :param od: observation data - an :class:`adaq_data.ADAQData` object
               with data in sites_cube_list
    :param md_list: model data - a list of :class:`adaq_data.ADAQData` objects
                    with data in sites_cube_list
    :param statsfile_prefix: string, output file to write stats data to.
                      If set to None, will write to 'plot_dir'/stats.csv
                      if 'plot_dir' is given in ini_dict.
                      If not available, no file will be written.
    :param thresholds: dictionary, or constant value.
                       If constant value (eg float, int), then the same
                       threshold is used for all phenomenon.
                       If dictionary, the keys should be either cube names
                       (first priority), or else short_name.
                       The values matched to these keys are then the thresholds
                       which are used for that key.
    :param nc_append: If writing to netcdf file, then if True, append to
                      file of the same name if it already exists. If False,
                      then overwrites the file instead.
    :param stats_tcoord: iris DimCoord representing time, which can be used
                         to overwrite default setting of the time coordinate
                         which will be generated for creating netcdf statistics
                         files from.

    :return: List of :class:`timeseries_stats.TimeSeriesStats` objects
             (dictionaries of statistics).

    >>> ini_dict, sites_data, od, md_list = get_exampledata()
    ... # doctest: +ELLIPSIS
    Reading inifile .../example_data_1days.ini
    Number of sites:  5
    >>> ini_dict = {} #Only used for plot_dir if want to write to file
    >>> thresholds = {'O3': 100.5, 'PM2p5': 35.5}
    >>> tsstats_list = calc_stats(ini_dict, od, md_list, thresholds=thresholds)
    >>> print(tsstats_list[0])  # doctest: +ELLIPSIS
    <...timeseries_stats.TimeSeriesStats object at ...>
    >>> print(sorted(tsstats_list[0].statsdict.keys()))
    ... # doctest: +NORMALIZE_WHITESPACE
    ['bias', 'correlation', 'fac2', 'falsealarmrate', 'falsealarmratio',
    'fge', 'hitrate', 'ioa', 'maxmod', 'maxobs', 'mdi', 'meanmod', 'meanobs',
    'mge', 'mnmb', 'nmb', 'nmge', 'npts', 'nsites', 'o<t_m<t', 'o<t_m>=t',
    'o>=t_m<t', 'o>=t_m>=t', 'odds_ratio', 'orss', 'pc95mod', 'pc95obs',
    'perc_correct', 'perc_over', 'perc_under', 'rmse', 'sdmod', 'sdobs', 'threshold', 'units']
    >>> print('{:.3f}'.format(tsstats_list[0].statsdict['bias']))
    -6.204
    '''

    percentile_value = ini_dict.get('percentile', 95)

    #If a single modeldata class passed in, convert this first to a list
    if not isinstance(md_list, list):
        md_list = list(md_list)

    tsstats_list = []

    if not isinstance(od.sites_cube_list, iris.cube.CubeList):
        #No cubelist available to work from
        return tsstats_list

    phenomenon_list = [cube.attributes['short_name']
                       for cube in od.sites_cube_list]
    phenomenon_list.sort()

    startdate, enddate = None, None
    for phenomenon in phenomenon_list:

        #Extract single species cubes for obs:
        obs_cube_full = od.extract(short_name=phenomenon, singlecube=True)
        if obs_cube_full is None:
            #Can not calculate statistics
            continue

        #Set up thresholds for statistics calculation if
        # a dictionary of thresholds is given
        if thresholds is None:
            threshold = None
        elif isinstance(thresholds, dict):
            if obs_cube_full.name() in thresholds:
                threshold = thresholds[obs_cube_full.name()]
            elif phenomenon in thresholds:
                threshold = thresholds[phenomenon]
            else:
                threshold = None
        elif len(thresholds) == 1:
            #Otherwise set the threshold as a constant value
            threshold = thresholds
        else:
            raise ValueError("thresholds argument should be None,\
                              constant value or a dictionary")

        #Loop over model data:
        # - extract single species cubes for model:
        # - check they have the same times
        # - get information about start and end times
        # - calculate statistics
        # - add extra information to the stats dictionary

        for md in md_list:

            mod_cube = md.extract(short_name=phenomenon, singlecube=True)
            if mod_cube is None:
                #Can not calculate statistics
                continue
            obs_cube, mod_cube = cube_time.intersect_cubetime([obs_cube_full,
                                                               mod_cube])
            startdt, enddt = cube_time.get_startenddt(cube=obs_cube)
            if startdate is None:
                startdate = startdt
            if startdt < startdate:
                startdt = startdate
            if enddate is None:
                enddate = enddt
            if enddt > enddate:
                enddate = enddt



            stats = timeseries_stats.TimeSeriesStats(obs_cube, mod_cube,
                                                     phenomenon=phenomenon,
                                                     threshold=threshold,
                                                     stats_tcoord=stats_tcoord)
            stats.get_basic_stats(percentile_value=percentile_value)
            stats.convert_to_cube()

            tsstats_list.append(stats)


    #Write out to file if required
    statsdir = ini_dict.get('plot_dir', None)
    if statsfile_prefix is not None or statsdir is not None:
        if tsstats_list:
            #Some stats in list, so can be written to file
            header = 'Statistics for '+\
                     startdate.strftime("%a %d-%b-%Y %H:%M")+\
                     ' - '+enddate.strftime("%a %d-%b-%Y %H:%M")

            format_list = ini_dict.get('calc_stats_format_list', ['csv'])

            timeseries_stats.save_stats(tsstats_list,
                                        percentile_value, filename_prefix=statsfile_prefix,
                                        format_list=format_list, directory=statsdir,
                                        header=header, nc_append=nc_append)

    return tsstats_list

def percent_contribution(md_list):
    '''
    Compute percentage contributions to integrated quantities
    such as air concentrations and deposits. This is done by
    dividing the value at each point by the total value over all
    points in the horizontal grid.

    .. note:: There is currently no support for calculating percent
              contributions for averaged quantities

    As an example, load in some model data and compute percentage
    contributions

    Use the sample data path to locate data

    >>> sample_data_path = config.SAMPLE_DATADIR+'name/'

    Read in the data

    >>> name = name_data.NAMEData()
    >>> name.readdata(sample_data_path + 'Fields_grid2_201304172000.txt',
    ... field_attributes = {'Quantity': 'Wet deposition'})
    [<iris 'Cube' of TRACER_WET_DEPOSITION / (g/m2) \
(latitude: 200; longitude: 200)>]

    >>> [name] = percent_contribution([name])
    >>> print(name.gridded_cube_list)
    0: Percent_Contribution / (1)          (latitude: 200; longitude: 200)
    '''

    for md in md_list:
        for mind, ngc in enumerate(md.gridded_cube_list):

            # Determine whether quantity is an average or integral
            if 'Time Av or Int' in ngc.attributes:
                av_or_int = ngc.attributes['Time Av or Int']
            else:
                av_or_int = ngc.attributes['Time av/int info']

            if 'average' in av_or_int:
                ve_message = 'Cannot compute percent contributions'
                ve_message += ' for average quantities'
                raise ValueError(ve_message)
                # For average quantities we need release rate and duration
                # which will require some work

            else:

                total_ngc = ngc.collapsed(['latitude', 'longitude'],
                                          iris.analysis.SUM)
                ngc_percent = 100. * (ngc / total_ngc)
                # Remove nans form the percent contribution cube
                xdata = np.nan_to_num(ngc_percent.data)
                ngc_percent.data = xdata

                for att in ['End of release', 'Forecast duration',
                            'Met data', 'NAME Version', 'Release height',
                            'Release location', 'Release rate', 'Run time',
                            'Start of release', 'Time Av or Int',
                            'Time av/int info', 'Run name',
                            'Title', 'label', 'Source strength']:
                    if att in ngc.attributes:
                        ngc_percent.attributes[att] = ngc.attributes[att]

            ngc_percent.attributes['short_name'] = 'Percent_Contribution'
            ngc_percent.attributes['Quantity'] = 'Percent Contribution'
            ngc_percent.rename('Percent_Contribution')
            md.gridded_cube_list[mind] = ngc_percent

    return md_list

if __name__ == "__main__":

    #To create new sample data files, uncomment the lines below and run.
    #The resulting netcdf files need to be copied into config.SAMPLE_DATADIR

    #create_example_data(exampletype='1days', save=True)
    #create_example_data(exampletype='5days', save=True)
    #create_example_data(exampletype='10days', save=True)
    #create_example_data(exampletype='3d', save=True)
    #create_example_data(exampletype='stats', save=True)

    warnings.filterwarnings("ignore") #ignore warnings when running
    import doctest
    doctest.testmod()
