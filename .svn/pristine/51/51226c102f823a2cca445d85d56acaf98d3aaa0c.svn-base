#!/usr/bin/env python
"""
Generalised script to generate AQ verification plots, by weather regime,
from PP and obs data.

Run using ./aq_plot_regimes.py [inifilename]
"""
from __future__ import print_function

from six.moves.builtins import str
from six.moves.builtins import range

import datetime
import os
import sys

#Find adaq directory -
# should be in directory above the one containing this file
adaq_path = os.path.dirname(os.path.realpath(__file__))+'/../'
sys.path.append(adaq_path)

from adaqcode.adaq_plotting import plot_md_gridded_fields
import adaqcode.weather_regime as weather_regime
import adaqcode.adaq_plotting as adaq_plotting
import adaqcode
from aq_plot import STATS_THRESHOLDS, SITE_CLASSIFICATIONS

def plot_gridded(ini_dict, md_list, short_name, filesuffix='', title='',
                 verbose=1):
    """
    Produce gridded field plots of regimes,
    plotting MEAN and MAX fields only.

    :param ini_dict: Dictionary of an :class:`inifile` object.
    :param od: observation data in the form of
               an :class:`adaq_data.ADAQData` object
    :param md_list: list of model data objects in the form of an
                    :class:`adaq_data.ADAQData` objects.
    :param short_name: string to match to short_name attribute in cubes.
    :param filesuffix: String to add to end of filename for specific
                       naming purposes. By default adds nothing.
    :param title: String to prefix main title with.
                  If prefix ends with ``\\n`` this will put it
                  on the line above the automatic title.
    :param verbose: Level of print output:

                    * 0 = No extra printing
                    * 1 = Standard level (recommended)
                    * 2 = Extra level for debugging
    """

    contour_levels = adaqcode.aq_indices.DAQI_LEVELS.get(short_name, None)
    ini_dict['levels_list'] = contour_levels

    for md in md_list:

        #Develop aggregated md
        cube = md.extract(short_name=short_name, gridded=True, singlecube=True)

        for agg_name in ["MEAN", "MAX"]:

            aggregator = adaqcode.cube_statistics.CUBE_AGGREGATORS[agg_name]
            aggregated_cube = cube.collapsed('time', aggregator)
            suffix_agg = filesuffix + '_' + agg_name.lower()
            title_agg = title + ' ' + agg_name.capitalize() + '\n'

            plot_md_gridded_fields(ini_dict, None, short_name,
                                   griddedcube=aggregated_cube,
                                   verbose=verbose,
                                   defaults='AQ',
                                   filesuffix=suffix_agg,
                                   titleprefix=title_agg)



def regime_plotting(inifilename=None, verbose=1):
    """
    Top-level plotting routine to incorporate all possible types of regime
    specific plotting.

    :param inifilename: String giving filename of ini file.
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
     * Adds regime data to cubes
     * Then where requested from inifile, produces regime-specific output of
       statistics and plots.

    >>> inifilename = adaq_path+'adaqscripts/aq_plot_regimes.ini'
    >>> ini_dict, sites_data, od, md_list = regime_plotting(inifilename,
    ... verbose=0) # doctest: +ELLIPSIS
    Beginning at  ...
    Reading inifile .../aq_plot_regimes.ini
    Number of sites:  5
    Getting regime data from  .../daily_30regimes_since2010.txt
    Creating observation data at ...
    Reading obs data files
    Found obs for  .../AURN_obs/ABD_20140101_20140818.txt
    Found obs for  .../AURN_obs/ACTH_20140101_20140818.txt
    Found obs for  .../AURN_obs/AH_20140101_20140818.txt
    Found obs for  .../AURN_obs/HAR_20140101_20140818.txt
    Found obs for  .../AURN_obs/YW_20140101_20140818.txt
    Creating obs cubes
    Getting model data for  oper  at ...
    Plotting regime bar chart
    Saved figure .../Regime_Bar.png
    Plotting regime time-series
    Saved figure .../Regime_Timeseries.png
    Site: ABD
    Plotting timeseries
    Site: ACTH
    Plotting timeseries
    ...
    Site: ABD
    Plotting timeseries
    ...
    Site: YW
    Plotting timeseries
    Plotting histogram
    Saved figure .../Histogram_O3_regime6.png
    Plotting Soccer Plot
    Saved figure .../Soccer_Plot_O3_regime6.png
    Plotting gridded fields
    Plotting gridded fields
    Plotting histogram
    Saved figure .../Histogram_PM2p5_regime6.png
    ...
    Plotting histogram
    Saved figure .../Histogram_O3_regime27.png
    ...
    Statistics saved to .../stats.csv
    Statistics found at .../stats.csv
    Saved figure .../meanobs_vs_mnmb_O3.png
    Saved figure .../meanobs_vs_fge_O3.png
    Saved figure .../meanobs_vs_bias_O3.png
    Saved figure .../meanobs_vs_meanmod_O3.png
    Saved figure .../meanobs_vs_mnmb_PM2p5.png
    Saved figure .../meanobs_vs_fge_PM2p5.png
    Saved figure .../meanobs_vs_bias_PM2p5.png
    Saved figure .../meanobs_vs_meanmod_PM2p5.png
    Saved figure .../meanobs_vs_mnmb_PM10.png
    Saved figure .../meanobs_vs_fge_PM10.png
    Saved figure .../meanobs_vs_bias_PM10.png
    Saved figure .../meanobs_vs_meanmod_PM10.png
    Finished at ...
    """

    print('Beginning at ', datetime.datetime.now())

    #Read ini file
    ini_dict = adaqcode.inifile.get_inidict(
        inifilename=inifilename,
        defaultfilename='adaqscripts/aq_plot_regimes.ini')

    #Ensure models_list is setup if not already given in inifile
    #Note - requires models_dir_list to be given instead.
    if 'models_list' not in ini_dict:
        ini_dict['models_list'] = [
            os.path.basename(path) if path is not None else None
            for path in ini_dict['models_dir_list']]

    #Set up some extra defaults if not set in ini_dict
    if ini_dict.get('contours', False):
        ini_dict['nlevels'] = ini_dict.get('nlevels', 11)
        ini_dict['cmap'] = ini_dict.get('cmap', adaqcode.aq_indices.daqi_cmap())
        ini_dict['cbar_orientation'] = ini_dict.get('cbar_orientation',
                                                    'vertical')
        ini_dict['extent_list'] = ini_dict.get('extent_list', None)

    #Get site data
    sites_data = adaqcode.sites_info.get_siteinfo(ini_dict)

    #Read in regime data
    regime_dates = adaqcode.weather_regime.read_regime_txt(
        ini_dict['regimes_file'])

    #Read in observations
    od = adaqcode.adaq_functions.get_obs(ini_dict, sites_data)
    #Read in model data
    md_list = adaqcode.adaq_functions.get_models(ini_dict, sites_data)

    #Unit conversion
    od, md_list = adaqcode.adaq_functions.unit_conversion(
        od, md_list, chem_units='ug/m3', aerosol_units='ug/m3')

    #Add weather regime as a coordinate to obs/model data cubes
    for cube in od.sites_cube_list:
        cube = weather_regime.add_regime_coord(cube, regime_dates, 1)
    for md in md_list:
        for cube in md.sites_cube_list:
            cube = weather_regime.add_regime_coord(cube, regime_dates, 1)
        for cube in md.gridded_cube_list:
            cube = weather_regime.add_regime_coord(cube, regime_dates, 0)

    #If 'regime_plots' not in ini_dict then set to a default of False
    #Only plot if 'regime_plots'=True
    if ini_dict.get('regime_plots', False):
        #Just pick first cube in list, first site in list - regime will be the
        #same for all of these.
        cube = od.sites_cube_list[0][0]

        #Produce bar chart showing frequency of regimes
        adaqcode.weather_regime.plot_regime_bar(cube, ini_dict['plot_dir'])

        #Timeseries of regime classification
        adaqcode.weather_regime.plot_regime_timeseries(cube, ini_dict['plot_dir'])

    for short_name in ini_dict['short_name_list']:
        #Produce site-specific plots
        if ini_dict.get('timeseries', False):
            for site in sites_data:
                print('Site:', site['abbrev'])

                #Plot time-series - only if set to True in inifile
                adaq_plotting.plot_timeseries(
                    ini_dict, od, md_list, short_name,
                    site, verbose=verbose)


    for md in md_list:

        md_newlist = []

        for iregime in range(1, 31):

            #Segregate cubes by regime
            od_regime = od.extract(regime=iregime)
            for cube in od_regime.sites_cube_list:
                cube.attributes['label'] = cube.attributes['label']+str(iregime)

            md_regime = md.extract(regime=iregime, gridded=False)
            md_regime = md_regime.extract(regime=iregime, gridded=True)

            for cube in md_regime.sites_cube_list:
                cube.attributes['label'] = cube.attributes['label']+str(iregime)

            #Produce species-specific and regime-specific plots
            filesuffix = '_regime'+str(iregime)
            title = 'Regime ' + str(iregime)
            for short_name in ini_dict['short_name_list']:
                if md_regime.sites_cube_list:
                    md_newlist.append(md_regime)


                    #Histogram plots - only if set to True in inifile
                    if ini_dict.get('histograms', False):
                        adaq_plotting.plot_histogram(ini_dict, od_regime,
                                                     [md_regime], short_name,
                                                     binsize=2.,
                                                     filesuffix=filesuffix)

                    #Soccer plots - only if set to True in inifile
                    if ini_dict.get('soccer_plots', False):
                        adaq_plotting.plot_soccer_plots(ini_dict, od_regime,
                                                        [md_regime], short_name,
                                                        SITE_CLASSIFICATIONS,
                                                        filesuffix=filesuffix)

                #Contour plots - only if set to True in inifile
                if md_regime.gridded_cube_list:
                    if ini_dict.get('contours', False):
                        plot_gridded(ini_dict, [md_regime], short_name,
                                     filesuffix, title, verbose=verbose)


        #Calculate regime-specific statistics and save to csv file
        # - only if set to True in inifile
        if ini_dict.get('calc_stats', False):
            #Get thresholds based on start of index 4 (maximum of index 3)
            adaqcode.adaq_functions.calc_stats(
                ini_dict, od, md_newlist,
                thresholds=STATS_THRESHOLDS)

        #Scatter plots of key statistics - only if set to True in inifile
        if ini_dict.get('regime_stats_scatter', False):
            weather_regime.plot_regime_stats(ini_dict)

    print('Finished at ', datetime.datetime.now())

    return ini_dict, sites_data, od, md_list


if __name__ == '__main__':

    #regime_plotting()

    #ini_dict, sites_data, od, md_list = regime_plotting('aq_plot_regimes.ini')

    import doctest
    doctest.testmod()
