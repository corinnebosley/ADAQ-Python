#!/usr/bin/env python
"""
Generalised script to generate AQ verification plots from model and obs data.

Copy the file 'adaqscripts/aq_plot.ini' to your working directory and modify it.

From the 'adaqscripts' directory, run via:

.. code-block:: ksh

    ./aq_plot.py [inifilename]

where inifilename includes the full path to the ini file.

Alternatively, from within other python code, run using:

.. code-block:: python

    from adaqscripts import aq_plot
    aq_plot.aq_plotting(inifilename='inifilename')

"""
from __future__ import print_function

import os
import sys
import datetime
import iris.cube
#Find adaq directory -
# should be in directory above the one containing this file
adaq_path = os.path.dirname(os.path.realpath(__file__))+'/../'
sys.path.append(adaq_path)
import adaqcode.adaq_functions as adaq_functions
import adaqcode.adaq_plotting as adaq_plotting
import adaqcode.aq_indices
import adaqcode.config
import adaqcode.inifile
import adaqcode.plotting_functions
import adaqcode.sites_info
import adaqcode.website

#: Thresholds based on start of DAQI index 4 (maximum of index 3)
STATS_THRESHOLDS = {k: v[3]
                    for k, v in adaqcode.aq_indices.DAQI_LEVELS.items()}

#: Site classifications (used for soccer plots and webpages)
SITE_CLASSIFICATIONS = {
    'Rural': {'site_type': ['REMOTE', 'RURAL', 'RURAL_BACKGROUND'],
              'marker': 'o'},
    'Urban': {'site_type': ['SUBURBAN', 'URBAN_BACKGROUND',
                            'SUBURBAN_BACKGROUND', 'URBAN_CENTRE'],
              'marker': 's'},
    'Traffic': {'site_type': ['ROADSIDE', 'KERBSIDE',
                              'RURAL_TRAFFIC', 'SUBURBAN_TRAFFIC',
                              'URBAN_TRAFFIC'],
                'marker': 'D'},
    'Other': {'site_type': ['RURAL_INDUSTRIAL', 'SUBURBAN_INDUSTRIAL',
                            'SUBURBAN_UNKNOWN', 'UNKNOWN_BACKGROUND',
                            'URBAN_INDUSTRIAL', 'URBAN_UNKNOWN'],
              'marker': '^'}
    }


def aq_get_data(inifilename=None, ini_dict=None, sites_data=None,
                keep_all_data=False, surface_only=True,
                defaultinifile='adaqscripts/aq_plot.ini'):
    """
    Function to get all the required data into memory:

    * Reads ini file and sets any extra required defaults
    * Calculates from ini file whether statistics will be produced
    * Reads in sites data
    * Reads in observations
    * Reads in model data
    * Deletes any data which won't be used further

    :param inifilename: String giving filename of ini file.
    :param ini_dict: Dictionary of a :class:`inifile` object. If this is given,
                     then it is used instead of reading in from inifilename.
    :param sites_data: numpy ndarray containing site information data
                       from a :class:`sites_info.SitesInfo` object.
                       If this is given, then it is used instead of reading in
                       from a file/cube.
    :param keep_all_data: If set to False, gridded_cube_lists are set to empty
                          iris cubelists after extracting sites_cube_list if no
                          longer required (eg for contours), and sites_cube_list
                          are set to empty iris cubelists if not required.
                          This leads to improvements in speed and memory as
                          routines are not being run on this unused data.
                          Set to True to keep all data - useful for debugging.
    :param surface_only: If set to True then only loads surface data, if False
                         then will also load 3D data.

    :returns: (ini_dict, sites_data, od, md_list)

    """

    #----
    #Read ini file
    #----

    if ini_dict is None:
        ini_dict = adaqcode.inifile.get_inidict(
            inifilename=inifilename,
            defaultfilename=defaultinifile)

    #----
    #Set up variables
    #----

    #Ensure models_list is setup if not already given in inifile
    #Note - requires models_dir_list to be given instead.
    #If models_dir_list is also not given, models_dir_list is set to []
    ini_dict['models_dir_list'] = ini_dict.get('models_dir_list', [])
    if 'models_list' not in ini_dict:
        ini_dict['models_list'] = [
            os.path.basename(path) if path is not None else None
            for path in ini_dict['models_dir_list']]

    #For ease of use later, figure out whether any statistics/timeseries/
    #fieldplotting will be produced.
    #Store these in ini_dict to allow easy passing around of variable.
    ini_dict['statistics'] = ini_dict.get('calc_stats', False) or \
       ini_dict.get('histograms', False) or \
       ini_dict.get('qqplots', False) or \
       ini_dict.get('soccer_plots', False) or \
       ini_dict.get('diurnal', False) or \
       ini_dict.get('timeseries_of_stats', False) or \
       ini_dict.get('timeseries_of_aggregate_stats', False)
    ini_dict['tseries'] = ini_dict.get('timeseries', False) or \
                          ini_dict.get('timeseries_multiple_short_names', False)
    #Turn off daqi_site_maps if daqi=False
    if 'daqi_site_maps' in ini_dict:
        ini_dict['daqi_site_maps'] = ini_dict.get('daqi', False) and \
                                     ini_dict.get('daqi_site_maps', False)
    ini_dict['fieldplot'] = ini_dict.get('contours', False) or \
                            ini_dict.get('daqi_site_maps', False)

    #Set up some extra defaults if not set in ini_dict
    if ini_dict.get('fieldplot', False):
        #Default to 11 levels (boundary points) to match DAQI 10 levels
        ini_dict['nlevels'] = ini_dict.get('nlevels', 11)
        #Default to DAQI colour map
        ini_dict['cmap'] = ini_dict.get('cmap', 'DAQI')
        #Default to vertical colour bar
        ini_dict['cbar_orientation'] = ini_dict.get('cbar_orientation',
                                                    'vertical')
        #Allow a fixed levels list to be used for all plots
        ini_dict['fixed_levels_list'] = ini_dict.get('levels_list', None)

    #Set up default location for corrected csv files for MARGA data if required
    #will write to plot_dir if set in ini_dict otherwise defaults to TEST_DIR
    if 'obsdir_fixed_csv' not in ini_dict:
        ini_dict['obsdir_fixed_csv'] = ini_dict.get(
            'plot_dir', adaqcode.config.TEST_DIR)+'/fixed_csv'

    #Determine which data can be deleted early
    #Don't need to keep the gridded data if only using site-specific
    # ie don't need to plot contours
    delete_gridded = not (keep_all_data
                          or ini_dict.get('contours', False)
                          or ini_dict.get('save_cubes_nc', False))
    #And delete the sites data if only gridded data is required for contours
    # ie don't need to plot any site-specific timeseries or statistics.
    delete_sites = not (keep_all_data
                        or ini_dict.get('statistics', False)
                        or ini_dict.get('tseries', False)
                        or ini_dict.get('daqi_site_maps', False))

    #----
    #Read data
    #----

    #Get site data (if reading from file - otherwise None returned)
    if sites_data is None:
        sites_data = adaqcode.sites_info.get_siteinfo(ini_dict)

    #Read in observations
    od = adaq_functions.get_obs(ini_dict, sites_data)
    if sites_data is None and od.sites_cube_list:
        #Get site data from observation site cube
        sites_data = adaqcode.sites_info.get_siteinfo(
            ini_dict,
            sites_cube=od.sites_cube_list[0])

    #Read in model data (deleting any gridded data if not required)
    md_list = adaq_functions.get_models(ini_dict, sites_data,
                                        delete_gridded=delete_gridded,
                                        surface_only=surface_only)

    #----
    #Remove unneeded data
    #----

    #Reset gridded_cube_list:
    if delete_gridded:
        print('Removing gridded_cube_list')
        od.gridded_cube_list = iris.cube.CubeList()
        for md in md_list:
            md.gridded_cube_list = iris.cube.CubeList()
    #Reset sites_cube_list:
    if delete_sites:
        print('Removing sites_cube_list')
        od.sites_cube_list = iris.cube.CubeList()
        for md in md_list:
            md.sites_cube_list = iris.cube.CubeList()


    return ini_dict, sites_data, od, md_list

def aq_prepare_data(ini_dict, od, md_list):
    """
    Function to do various preparatory functions on the data:

    * Apply any scaling factors
    * Unit conversion to ug/m3
    * Add missing times
    * Match data ready for statistics
    * Conversion to DAQI

    :param ini_dict: Dictionary of a :class:`inifile` object.
    :param od: observation data - an :class:`adaq_data.ADAQData` object
               with data in sites_cube_list
    :param md_list: model data - a list of :class:`adaq_data.ADAQData` objects
                    with data in sites_cube_list

    :returns: (ini_dict, od, md_list, od_stats, md_list_stats)

    """

    #Apply scaling factors
    adaq_functions.scaling_factors(ini_dict, [od]+md_list)

    #Unit conversion
    od, md_list = adaq_functions.unit_conversion(
        od, md_list,
        chem_units=ini_dict.get('chem_units', 'ug/m3'),
        aerosol_units=ini_dict.get('aerosol_units', 'ug/m3'),
        pollen_units=ini_dict.get('pollen_units', 'grains m-3'))

    #Add missing times in to get regular times in arrays.
    od, md_list = adaq_functions.add_missing_times(od, md_list)

    #If using forecast_day='all', only works upto here
    #- match data for stats expects to be able to just get a single cube
    # for each species, but there could be upto 6 cubes per species.

    #Match the data ready for calculating statistics.
    if ini_dict.get('statistics', False):
        #Only need to calculate if actually using statistics
        od_stats, md_list_stats = \
                  adaq_functions.match_data_for_stats(ini_dict, od, md_list)
    else:
        od_stats = None
        md_list_stats = None

    #Convert to DAQI if required
    if ini_dict.get('daqi', False):
        ini_dict, od, md_list, od_stats, md_list_stats = \
                  adaq_functions.convert_to_daqi(
                      ini_dict, od, md_list, od_stats, md_list_stats)

    #Save prepared data as netcdf if requested
    if ini_dict.get('save_cubes_nc', False):
        output_dir = ini_dict.get('plot_dir', './')
        try:
            #note: when we only need to support python 3, can use
            # keyword argument exist_ok=True instead of try..except
            os.makedirs(output_dir)
        except OSError:
            pass

        suffix = ini_dict.get('save_cubes_label', '')
        if suffix and not suffix.startswith("_"):
            suffix = "_"+suffix

        if od:
            od.save_all(output_dir, prefix="obs_", suffix=suffix)
        if md_list:
            for data, label in zip(md_list, ini_dict['models_list']):
                data.save_all(output_dir, prefix=label+"_", suffix=suffix)
        suffix = "_stats"+suffix
        if od_stats:
            od_stats.save_all(output_dir, prefix="obs_", suffix=suffix)
        if md_list_stats:
            for data, label in zip(md_list_stats, ini_dict['models_list']):
                data.save_all(output_dir, prefix=label+"_", suffix=suffix)


    return ini_dict, od, md_list, od_stats, md_list_stats

def aq_stats_and_plots(ini_dict, sites_data, od, md_list,
                       od_stats, md_list_stats, verbose=1):
    """
    Function to call all required statistics and plotting routines,
    depending on values set in ini_dict.

    :param ini_dict: Dictionary of a :class:`inifile` object.
    :param sites_data: numpy ndarray containing site information data
                       from a :class:`sites_info.SitesInfo` object.
    :param od: observation data - an :class:`adaq_data.ADAQData` object
               with data in sites_cube_list
    :param md_list: model data - a list of :class:`adaq_data.ADAQData` objects
                    with data in sites_cube_list
    :param od_stats: observation data which has been matched to model data
                     ready for statistic calculations - an
                     :class:`adaq_data.ADAQData` object with data in
                     sites_cube_list
    :param md_list_stats: model data which has been matched to observation data-
                    a list of :class:`adaq_data.ADAQData` objects with data in
                    sites_cube_list
    :param verbose: Level of print output:

                    * 0 = No extra printing
                    * 1 = Standard level (recommended)
                    * 2 = Extra level for debugging

    """
    # Calculate statistics and save to csv file
    #If 'calc_stats' not in ini_dict then set to a default of False
    #Only calculate stats if 'calc_stats'=True
    if ini_dict.get('calc_stats', False):
        adaq_functions.calc_stats(ini_dict, od_stats, md_list_stats,
                                  thresholds=STATS_THRESHOLDS,
                                  nc_append=True)



    # Produce species-specific plots
    for short_name in ini_dict['short_name_list']:

        print('Short_name:', short_name)

        #--- Field plot data

        if ini_dict['fieldplot']:

            #Set up contour levels
            if ini_dict['fixed_levels_list'] is None:
                if ini_dict['cmap'] == 'CAMS':
                    #Get levels from CAMS_LEVELS if available
                    contour_levels = adaqcode.aq_indices.CAMS_LEVELS.get(
                        short_name, None)
                else:
                    #By default take from DAQI_LEVELS if available
                    contour_levels = adaqcode.aq_indices.DAQI_LEVELS.get(
                        short_name, None)
                # Check units - for CAMS or DAQI levels, should be in ug/m3
                chem_units = ini_dict.get('chem_units', 'ug/m3')
                aerosol_units = ini_dict.get('aerosol_units', 'ug/m3')
                if chem_units != 'ug/m3' or aerosol_units != 'ug/m3':
                    cube = md_list[0].extract(short_name=short_name,
                                              gridded=True, singlecube=True)
                    if cube and cube.units != 'ug/m3':
                        contour_levels = None

            else:
                contour_levels = ini_dict['fixed_levels_list']
            ini_dict['levels_list'] = contour_levels

        if ini_dict.get('contours', False):
            adaq_plotting.plot_md_gridded_fields(ini_dict, md_list, short_name,
                                                 defaults='AQ',
                                                 plot_subdir='gridded_fields',
                                                 verbose=verbose)

        # Maps of DAQI sites data
        if 'DAQI' in short_name and ini_dict.get('daqi', False) and \
           ini_dict.get('daqi_site_maps', False):
            #Increase default marker size to 40 (instead of 20)
            ini_dict['marker_size'] = ini_dict.get('marker_size', 40)
            adaq_plotting.plot_sitescube_maps(
                ini_dict, od, md_list, short_name,
                SITE_CLASSIFICATIONS,
                plot_subdir='gridded_fields',
                tight=True,
                verbose=verbose)


        #--- Site specific data

        # Histogram plots - only if set to True in inifile
        if ini_dict.get('histograms', False):
            #Force binsize to 2 for main AQ species
            if short_name in ['O3', 'NO2', 'SO2', 'PM10', 'PM2.5']:
                binsize = 2.
                if short_name in ['O3', 'NO2', 'SO2'] and \
                   ini_dict.get('chem_units', 'ug/m3') == 'cm-3':
                    #Automatically calculate binsize
                    binsize = None
            elif short_name[-4:] == 'DAQI':
                binsize = 1
            else:
                #Automatically calculate binsize
                binsize = None
            adaq_plotting.plot_histogram(ini_dict, od_stats, md_list_stats,
                                         short_name, binsize=binsize,
                                         maxperc=99.5)

        # Quantile-Quantile plots - only if set to True in inifile
        if ini_dict.get('qqplots', False):
            adaq_plotting.plot_qq(ini_dict, od_stats, md_list_stats, short_name)


        # Soccer plots - only if set to True in inifile
        if ini_dict.get('soccer_plots', False):
            adaq_plotting.plot_soccer_plots(ini_dict, od_stats, md_list_stats,
                                            short_name, SITE_CLASSIFICATIONS)

        # Diurnal plots
        if ini_dict.get('diurnal', False):
            #Doesn't make sense to produce hourly plots for daily values
            if short_name[-4:] != 'DAQI':
                adaq_plotting.plot_diurnal(ini_dict, od_stats, md_list_stats,
                                           short_name)

        # Timeseries of statistics
        if ini_dict.get('timeseries_of_stats', False):
            threshold = STATS_THRESHOLDS.get(short_name, None)
            adaq_plotting.plot_timeseries_of_stats(ini_dict, od_stats,
                                                   md_list_stats,
                                                   short_name,
                                                   threshold=threshold)

        #Timeseries of aggregate statistics across sites
        if ini_dict.get('timeseries_of_aggregate_stats', False):
            for stat in ini_dict.get('timeseries_of_aggregate_stats_list', []):
                adaq_plotting.plot_timeseries_aggregate_stats(
                    ini_dict,
                    [od_stats]+md_list_stats,
                    short_name,
                    stat=stat,
                    linestyles=['--']+['-']*len(md_list_stats))

    # Produce site-specific timeseries plots
    if ini_dict.get('tseries', False):
        #Only loop through all sites if 'tseries' set to True
        #- set up in aq_get_data()
        for site in sites_data:
            print('Site:', site['abbrev'])

            #Plot time-series - only if set to True in inifile
            if ini_dict.get('timeseries', False):
                for short_name in ini_dict['short_name_list']:
                    print('Short_name:', short_name)

                    adaq_plotting.plot_timeseries(ini_dict, od, md_list,
                                                  short_name, site,
                                                  verbose=verbose)

            #Produce multiple-species timeseries
            if ini_dict.get('timeseries_multiple_short_names', False):
                tseries_multiple_sn_dict = ini_dict.get(
                    'timeseries_multiple_short_names_dict', {})

                for group, sname_list in tseries_multiple_sn_dict.items():
                    print('Multiple short_names group:', group)
                    adaq_plotting.plot_tseries_multiple_snames(
                        ini_dict, od, md_list,
                        site, sname_list, group,
                        verbose=verbose)


def aq_write_websites(ini_dict, sites_data):
    """
    Function to write the required websites

    :param ini_dict: Dictionary of a :class:`inifile` object. If this is given,
                     then it is used instead of reading in from inifilename.
    :param sites_data: numpy ndarray containing site information data
                       from a :class:`sites_info.SitesInfo` object.
                       If this is given, then it is used instead of reading in
                       from a file/cube.
    """

    if ini_dict.get('html_dir', None) is not None:
        if ini_dict.get('tseries', False):

            for classification in SITE_CLASSIFICATIONS:
                site_types = SITE_CLASSIFICATIONS[classification]['site_type']
                if ini_dict.get('timeseries', False):
                    filename = 'aq_plot_tseries_'+classification+'.html'
                    adaqcode.website.write_aq_plot_tseries_website(
                        ini_dict,
                        sites_data,
                        site_types=site_types,
                        filename=filename)
                if ini_dict.get('timeseries_multiple_short_names', False):
                    filename = ('aq_plot_tseries_multiple_species_'
                                +classification+'.html')
                    #Set up short_names using the keys from the
                    #timeseries_multiple_short_names_dict
                    #so only these are included on webpage
                    short_names = [
                        sn + '*' for sn in
                        list(ini_dict['timeseries_multiple_short_names_dict'].keys())]
                    adaqcode.website.write_aq_plot_tseries_website(
                        ini_dict,
                        sites_data,
                        site_types=site_types,
                        short_names=short_names,
                        filename=filename)

        if ini_dict.get('statistics', False):
            adaqcode.website.write_aq_plot_stats_website(ini_dict)

        if ini_dict.get('contours', False):
            adaqcode.website.write_aq_plot_contours_website(ini_dict)

        if ini_dict.get('daqi_site_maps', False) or \
           (ini_dict.get('contours', False) and ini_dict.get('daqi', False)):
            #Note also produce these pages if DAQI contour maps are plotted
            #This will redo the pages if contours are plotted separately to
            #site maps (contours will probably take longer to run).
            adaqcode.website.write_aq_plot_sites_website(ini_dict,
                                                         daqi_only=True)


def aq_plotting(inifilename=None, ini_dict=None, sites_data=None,
                keep_all_data=False, verbose=1):
    """
    Top-level plotting routine to incorporate all possible types of plotting.

    :param inifilename: String giving filename of ini file.
    :param ini_dict: Dictionary of a :class:`inifile` object. If this is given,
                     then it is used instead of reading in from inifilename.
    :param sites_data: numpy ndarray containing site information data
                       from a :class:`sites_info.SitesInfo` object.
                       If this is given, then it is used instead of reading in
                       from a file/cube.
    :param keep_all_data: If set to False, gridded_cube_lists are set to empty
                          iris cubelists after extracting sites_cube_list if no
                          longer required (eg for contours), and sites_cube_list
                          are set to empty iris cubelists if not required.
                          This leads to improvements in speed and memory as
                          routines are not being run on this unused data.
                          Set to True to keep all data - useful for debugging.
    :param verbose: Level of print output:

                    * 0 = No extra printing
                    * 1 = Standard level (recommended)
                    * 2 = Extra level for debugging

    :returns: (ini_dict, sites_data, od, md_list)

     * **ini_dict** - Dictionary of an :class:`inifile` object.
     * **sites_data** - numpy ndarray containing site information data
       from a :class:`sites_info.SitesInfo` object.
     * **od** - Observation Data -
       Either AURN observations from :class:`aurn_data`
       or CAMSAQ observations from :class:`camsaqobs_data`
       (both in :class:`adaq_data` format).
     * **md_list** - Model Data list - each item could be pp or nimrod in
       :class:`adaq_data` format.

    **Method:**
     * Reads in inifile
     * Reads in sites file (as defined in inifile)
     * Reads in observation data
     * Reads in model data
     * Does any required phenomenon and unit conversion.
     * Matches up data points between model and observations
     * Then where requested from inifile, produces output of
       statistics and plots.

    >>> ini_dict, sites_data, od, md_list = aq_plotting(verbose=0)
    ... # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    Current version of python: ...
    [...]
    Current version of iris: ...
    Beginning at  ...
    Reading inifile .../aq_plot.ini
    Number of sites:  5
    Creating observation data at ...
    Reading obs data files
    Found obs for .../AURN_obs/ABD_20140101_20140818.txt
    Found obs for .../AURN_obs/ACTH_20140101_20140818.txt
    Found obs for .../AURN_obs/AH_20140101_20140818.txt
    Found obs for .../AURN_obs/HAR_20140101_20140818.txt
    Found obs for .../AURN_obs/YW_20140101_20140818.txt
    Creating obs cubes
    Getting model data for  oper  at ...
    Getting model data for  name_casestudy  at ...
    Adding missing times
    Converting to DAQI
    Saving to .../obs_sites_cube_list.nc
    Saving to .../oper_sites_cube_list.nc
    Saving to .../oper_gridded_cube_list.nc
    Saving to .../name_casestudy_sites_cube_list.nc
    Saving to .../name_casestudy_gridded_cube_list.nc
    Saving to .../obs_sites_cube_list_stats.nc
    Saving to .../oper_sites_cube_list_stats.nc
    Saving to .../name_casestudy_sites_cube_list_stats.nc
    Statistics saved to  .../stats.csv
    Statistics saved to  .../stats.wiki
    Short_name: O3
    Plotting gridded fields
    Plotting histogram
    Saved figure  .../Histogram_O3.png
    Plotting quantile-quantile plot
    Saved figure  .../Quantile-Quantile_O3.png
    Plotting Soccer Plot
    Saved figure  .../Soccer_Plot_O3.png
    Plotting diurnal
    Saved figure  .../Diurnal_O3.png
    Plotting timeseries of statistics
    Saved figure  .../Timeseries_of_bias_O3.png
    Saved figure  .../Timeseries_of_rmse_O3.png
    Plotting timeseries of aggregated cube statistics
    Saved figure  .../Timeseries_of_sites_nanmean_O3.png
    Plotting timeseries of aggregated cube statistics
    Saved figure  .../Timeseries_of_sites_nanmax_O3.png
    Short_name: PM2p5
    Plotting gridded fields
    Plotting histogram
    Saved figure  .../Histogram_PM2p5.png
    Plotting quantile-quantile plot
    Saved figure  .../Quantile-Quantile_PM2p5.png
    Plotting Soccer Plot
    Saved figure .../Soccer_Plot_PM2p5.png
    Plotting diurnal
    Saved figure  .../Diurnal_PM2p5.png
    Plotting timeseries of statistics
    Saved figure  .../Timeseries_of_bias_PM2p5.png
    Saved figure  .../Timeseries_of_rmse_PM2p5.png
    Plotting timeseries of aggregated cube statistics
    Saved figure  .../Timeseries_of_sites_nanmean_PM2p5.png
    Plotting timeseries of aggregated cube statistics
    Saved figure  .../Timeseries_of_sites_nanmax_PM2p5.png
    Short_name: O3_DAQI
    Plotting gridded fields
    Plotting sitescube maps
    Plotting histogram
    Saved figure  .../Histogram_O3_DAQI.png
    Plotting quantile-quantile plot
    Saved figure  .../Quantile-Quantile_O3_DAQI.png
    Plotting Soccer Plot
    Saved figure  .../Soccer_Plot_O3_DAQI.png
    Plotting timeseries of statistics
    Plotting timeseries of aggregated cube statistics
    Saved figure  .../Timeseries_of_sites_nanmean_O3_DAQI.png
    Plotting timeseries of aggregated cube statistics
    Saved figure  .../Timeseries_of_sites_nanmax_O3_DAQI.png
    Short_name: PM2p5_DAQI
    ...
    Short_name: DAQI
    ...
    Site: ABD
    Short_name: O3
    Plotting timeseries
    Short_name: PM2p5
    Plotting timeseries
    Short_name: O3_DAQI
    Plotting timeseries
    Short_name: PM2p5_DAQI
    Plotting timeseries
    Short_name: DAQI
    Plotting timeseries
    Multiple short_names group: O3+PM2p5
    Plotting timeseries for multiple short_names
    Site: ACTH
    Short_name: O3
    Plotting timeseries
    ...
    Multiple short_names group: O3+PM2p5
    Plotting timeseries for multiple short_names
    Site: AH
    Short_name: O3
    Plotting timeseries
    ...
    Saved website  .../aq_plot_tseries_Rural.html
    Saved website  .../aq_plot_tseries_multiple_species_Rural.html
    Saved website  .../aq_plot_tseries_Urban.html
    Saved website  .../aq_plot_tseries_multiple_species_Urban.html
    Saved website  .../aq_plot_stats.html
    Saved website  .../aq_plot_gridded_fields_oper.html
    Saved website  .../aq_plot_gridded_fields_name_casestudy.html
    Saved website  .../aq_plot_site_fields_O3_DAQI.html
    Saved website  .../aq_plot_site_fields_PM2p5_DAQI.html
    Saved website  .../aq_plot_site_fields_DAQI.html
    Finished at  ...
    """

    print('Current version of python:', sys.version)
    print('Current version of iris:', iris.__version__)
    print('Beginning at ', datetime.datetime.now())

    #Read in all required data
    ini_dict, sites_data, od, md_list = aq_get_data(
        inifilename=inifilename,
        ini_dict=ini_dict,
        sites_data=sites_data,
        keep_all_data=keep_all_data)

    #Prepare all the data as required
    ini_dict, od, md_list, od_stats, md_list_stats = aq_prepare_data(
        ini_dict, od, md_list)

    #Produce statistics and all plots
    aq_stats_and_plots(ini_dict, sites_data,
                       od, md_list,
                       od_stats, md_list_stats,
                       verbose=verbose)

    #Write websites
    aq_write_websites(ini_dict, sites_data)

    print('Finished at ', datetime.datetime.now())

    return ini_dict, sites_data, od, md_list


if __name__ == '__main__':

    aq_plotting()

    #ini_dict, sites_data, od, md_list = aq_plotting(inifilename='aq_plot.ini')

    #import doctest
    #doctest.testmod()

    #doctest.run_docstring_examples(aq_plotting, globals())
