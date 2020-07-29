'''
This function loads 3D model data (either aqeur or aqcops) for
temperature, pressure and species mass mixing ratio (mmr).
Dobson Units are calculated from the species mmr. Contour
plots are generated for Dobson Units and surface species
concentration.
'''

from __future__ import print_function
import os
import sys
import datetime
import numpy as np
import iris

adaq_path = os.path.dirname(os.path.realpath(__file__))+'/../'
sys.path.append(adaq_path)

import adaqscripts
import adaqcode
import adaqcode.inifile
import adaqcode.constants
from adaqcode.cube_functions import get_latlon_cube


def mmr2dobsonunits(inifilename=None):

    """
    This function loads 3D model data (either aqeur or aqcops) for
    temperature, pressure and species concentration. It loads this data
    as a model data list and calls the ADAQ_Python function
    calc_dobson_units to determine the Dobson units for the required
    species within a user specified vertical range that is defined by
    the surface and p_target_top (upper limit). The cubes returned from
    calc_dobson_units contain gridded Dobson unit data and species
    concentration which are used for diagnostic plotting. \
    Note (1): If Dobson units are being calculated for more than one
    species the pressure and temperature are partially overwritten with
    NaN values depending on the value of p_target_top. \
    Note (2): Due to the size of the 3D array for AQCOPS output, if a
    longer time series is being processed (approx. one month or greater),
    only one satellite pass time can be computed when the model data is
    passed to calc_dobson_units.

    :param inifilename: filename of inifile. If set to None, taken from \
    command line argument if passed in. Required

    :returns ini_dict, md_time_extr, du_time_ave_cubelist: ini_dict is a \
    dictionary of an inifile object; md_time_extr is a list of 4D cubes of \
    temperature, pressure and species where the time coordinate is the \
    specified satellite overpass time and du_time_ave_cubelist is a list \
    of cubes containing the calculated, time averaged Dobson units for the \
    required species.


    Examples:

    >>> ini_dict, md_time_extr, du_time_ave_cubelist = mmr2dobsonunits(\
inifilename=adaq_path + 'adaqscripts/mmr2dobsonunits.ini') #doctest: +ELLIPSIS
    Reading inifile .../adaqscripts/mmr2dobsonunits.ini
    Getting model data for  3d  at  ...
    Extracting temperature cube
    Extracting pressure cube
    Extracting species cube
    Calculating Dobson units for species
    Saved figure .../Fieldplot_pp_3d_O3_surf_201707020900.png
    Saved figure .../Fieldplot_pp_3d_O3_DU_201707020900.png
    Saved figure .../Fieldplot_pp_3d_O3_surf_201707021200.png
    Saved figure .../Fieldplot_pp_3d_O3_DU_201707021200.png

    >>> print(md_time_extr.gridded_cube_list)
    0: air_pressure / (Pa)\
                 (time: 2; model_level_number: 63; grid_latitude: 182; grid_longitude: 146)
    1: air_temperature / (K)\
               (time: 2; model_level_number: 63; grid_latitude: 182; grid_longitude: 146)
    2: mass_fraction_of_ozone_in_air / (kg kg-1)\
 (time: 2; model_level_number: 63; grid_latitude: 182; grid_longitude: 146)

    >>> print(du_time_ave_cubelist[0].name(),
    ... du_time_ave_cubelist[0].attributes['short_name'],
    ... du_time_ave_cubelist[0].units)
    Tropospheric_ozone_column_abundance_in_Dobson_Units O3_DU DU

    >>> print(du_time_ave_cubelist) #doctest: +ELLIPSIS
    0: Tropospheric_ozone_column_abundance_in_Dobson_Units / \
(DU) (grid_latitude: 182; grid_longitude: 146)
    1: Tropospheric_ozone_column_abundance_in_Dobson_Units / \
(DU) (grid_latitude: 182; grid_longitude: 146)

    >>> print(du_time_ave_cubelist[0]) #doctest: +ELLIPSIS
    Tropospheric_ozone_column_abundance_in_Dobson_Units / \
(DU) (grid_latitude: 182; grid_longitude: 146)
         Dimension coordinates:
              grid_latitude                                                  x                    -
              grid_longitude                                                 -                    x
         Auxiliary coordinates:
              surface_altitude                                               x                    x
         Derived coordinates:
              altitude                                                       x                    x
         Scalar coordinates:
              forecast_day: 1.0 Days
              level_height: 21324.0... m, bound=(0.0, 42648...) m
              sigma: 0.5, bound=(0.0, 1.0)
              time: 2017-07-02 09:00:00, bound=(2017-07-02 08:00:00, 2017-07-02 09:00:00)
         Attributes:
              STASH: m01s34i001
              label: 3d
              short_name: O3_DU
              source: Data from Met Office Unified Model
         Cell methods:
              mean: time (1 hour)
              sum: model_level_number

    >>> print('{:.2f} {:.2f}'.format(np.max(du_time_ave_cubelist[0].data),
    ... np.min(du_time_ave_cubelist[0].data)))
    27.72 10.94

    """


    ## Load an inidict
    ini_dict = adaqcode.inifile.get_inidict(
        inifilename=inifilename,
        defaultfilename='adaqscripts/mmr2dobsonunits.ini')

    ## Extract the short_name_list
    short_name_list = ini_dict['short_name_list']

    ## If required, retrieve files from MASS
    if ini_dict['mass_retrieve'] is True:

        print("Retrieving model output from MASS")

        ## Initially get the orography from prods files
        ini_dict['short_name_list'] = ["orog"]
        adaqscripts.aq_mass_retrieve.mass_retrieval(ini_dict=ini_dict)

        ## Then get the rest of the model data
        ini_dict['forecast_day'] = '1'
        ini_dict['short_name_list'] = short_name_list
        ini_dict['mass_filenames'] = "*prodm*"
        adaqscripts.aq_mass_retrieve.mass_retrieval(ini_dict=ini_dict)


    ## Load the model data
    if 'models_list' not in ini_dict:
        ini_dict['models_list'] = [
            os.path.basename(path) if path is not None else None
            for path in ini_dict['models_dir_list']]
    md = adaqcode.adaq_functions.get_models(ini_dict, None,
                                            delete_gridded=False,
                                            surface_only=False)[0]


    ## Extract the satellite pass time/s
    sat_pass_times = [int(t) for t in ini_dict['sat_pass_time_list']]
    md_time_extr = md.extract(time=lambda cell: cell.point.hour in sat_pass_times, gridded=True)


    ## Remove the original model list to free memory
    md = None


    ## Calculate Dobson units for chemical species

    ## Set up a list for the cubes of time averaged
    ## Dobson units that will be generated for the
    ## required species
    du_time_ave_list = []

    ## Loop over the required species
    for short_name in short_name_list:

        if not short_name in ["p", "T"]:
            du_cube, sp_cube = adaqcode.cube_chemistry.calc_dobson_units(
                md_time_extr,
                np.float(ini_dict['p_target_top']),
                short_name)


        ## Plotting - Contour plots for species concentration at the
        ## surface and Dobson Units as calculated above.
        ## If more than one satellite pass time has been extracted,
        ## data can be plotted for each time point or averaged over
        ## all time points.

            if ini_dict['plot_single_time'] is True:

                for time_point in sat_pass_times:
                    tconst = iris.Constraint(time=lambda cell: cell.point.hour == time_point)

                    ## Plot surface concentration
                    sp_surf = sp_cube.extract(iris.Constraint(model_level_number=1))
                    sp_surf_tp = sp_surf.extract(tconst)

                    ## Find mean across time points if there is more than one time point
                    if len(sp_surf_tp.coord('time').points) > 1:
                        sp_surf_time_ave = sp_surf_tp.collapsed('time', iris.analysis.MEAN)
                    else:
                        sp_surf_time_ave = sp_surf_tp

                    sp_surf_short_name = short_name + '_surf'
                    sp_surf_time_ave.attributes['short_name'] = sp_surf_short_name


                    ## Plot Dobson Units
                    du_cube_tp = du_cube.extract(iris.Constraint(
                        time=lambda cell: cell.point.hour == time_point))

                    ## Find mean across time points if there is more than one time point
                    if len(du_cube_tp.coord('time').points) > 1:
                        du_time_ave = du_cube_tp.collapsed('time', iris.analysis.MEAN)
                    else:
                        du_time_ave = du_cube_tp


                    ## Append time averaged Dobson unit data to a cube list
                    du_time_ave_list.append(du_time_ave)


                    ## Plot contours if required
                    if ini_dict['contours'] is True:
                        __plot_contours(sp_surf_time_ave, du_time_ave, ini_dict)


            else:

                print('plot_single_time is False')

                sp_surf = sp_cube.extract(iris.Constraint(model_level_number=1))

                ## Find mean across time points if there is more than one time point
                if len(sp_surf.coord('time').points) > 1:
                    sp_surf_time_ave = sp_surf.collapsed('time', iris.analysis.MEAN)
                else:
                    sp_surf_time_ave = sp_surf

                sp_surf_short_name = short_name + '_surf'
                sp_surf_time_ave.attributes['short_name'] = sp_surf_short_name

                ## Find mean across time points if there is more than one time point
                if len(du_cube.coord('time').points) > 1:
                    du_time_ave = du_cube.collapsed('time', iris.analysis.MEAN)
                else:
                    du_time_ave = du_cube


                ## Plot contours if required
                if ini_dict['contours'] is True:
                    __plot_contours(sp_surf_time_ave, du_time_ave, ini_dict)


                ## Append time averaged Dobson unit data to a cube list
                du_time_ave_list.append(du_time_ave)


    ## Convert the list to a cube list
    du_time_ave_cubelist = iris.cube.CubeList(du_time_ave_list)


    return ini_dict, md_time_extr, du_time_ave_cubelist



def __plot_contours(sp_surf_time_ave, du_time_ave, ini_dict):

    """
    This function plots contour plots of Dobson Units
    and surface concentration for the same species. The
    Dobson Units and surface concentration can be
    averaged for individual satellite pass times (if
    more than one time point is given), or averaged over
    all time points in the cube.
    This function can only be called from mmr2dobsonunits.
    """

    # Set up plotting extent depending on model domain
    if ini_dict['models_op_list'] == ['aqeur']:
        map_extent = [-12.0, 4.0, 48.0, 60.0]

    if ini_dict['models_op_list'] == ['aqcops']:
        map_extent = [-25.0, 45.0, 30.0, 70.0]

    # Convert to ug/m3
    if ini_dict.get('units', 'ug/m3'):
        # Try to convert units to ug/m3 if possible
        try:
            sp_surf_time_ave = adaqcode.cube_chemistry.\
                cube_elemental_to_molecular_mass(sp_surf_time_ave)
        except ValueError:
            pass
        sp_surf_time_ave = adaqcode.cube_chemistry.cube_to_mass_conc_in_air(
            sp_surf_time_ave, at_stp=True, ounits='ug/m3')

        # Use sensible contour levels based on DAQI/other sensible values
        sp_surf_levels = adaqcode.aq_indices.DAQI_LEVELS.get(
            sp_surf_time_ave.attributes['short_name'].split('_')[0], None)
        if sp_surf_time_ave.attributes['short_name'] == 'NO2_surf':
            sp_surf_levels = np.array(([0., 2., 4., 6., 8., 10., 15., 20.,
                                        30., 40., 70.]))
    else:
        sp_surf_levels = None

    # Plot surface ozone
    if sp_surf_levels is None:
        nlevels = 10
    else:
        nlevels = len(sp_surf_levels)
    fl = adaqcode.field_layer.FieldLayer(sp_surf_time_ave)
    fl.set_layerstyle(nlevels=nlevels,
                      levels=sp_surf_levels,
                      plottype='pcolormesh',
                      cmap='DAQI')

    fp = adaqcode.field_plot.FieldPlot(ini_dict)
    fp.title = ("Z+".join(ini_dict['sat_pass_time_list'])+'Z\n'
                + ini_dict['start_datetime'].strftime("%d/%m/%d %HZ")
                + ' - ' + ini_dict['end_datetime'].strftime("%d/%m/%d %HZ"))
    fp.add_layer(fl)
    fp.setup_mapping(extent=map_extent,
                     mapping='coastlines',
                     gridlines=False)

    fp.plot()
    fp.save_fig(plotdir=ini_dict['plot_dir'], fileprefix='Fieldplot_pp_',)


    ## Plot ozone Dobson units
    if du_time_ave.attributes['short_name'] == 'O3_DU':
        levels = np.arange(0.0, 44.0, 2.0)
        nlevels = len(levels)

    if du_time_ave.attributes['short_name'] == 'NO2_DU':
        levels = np.arange(0.0, 0.55, 0.025)
        nlevels = len(levels)


    fl = adaqcode.field_layer.FieldLayer(du_time_ave)
    fl.set_layerstyle(nlevels=nlevels,
                      levels=levels,
                      plottype='pcolormesh',
                      cmap='nipy_spectral')

    fp = adaqcode.field_plot.FieldPlot(ini_dict)
    fp.title = ("Z+".join(ini_dict['sat_pass_time_list'])+'Z\n'
                + ini_dict['start_datetime'].strftime("%d/%m/%d %HZ")
                + ' - ' + ini_dict['end_datetime'].strftime("%d/%m/%d %HZ"))
    fp.add_layer(fl)
    fp.setup_mapping(extent=map_extent,
                     mapping='coastlines',
                     gridlines=False)
    fp.plot()
    fp.save_fig(plotdir=ini_dict['plot_dir'], fileprefix='Fieldplot_pp_',)


def convert_to_reg_latlon_nc(inifilename=None):
    """
    Convert data to a regular lat-lon grid and save to netcdf
    in format required for RAL.

    >>> inifilename=adaq_path + 'adaqscripts/mmr2dobsonunits.ini'
    >>> ini_dict, md_ll = convert_to_reg_latlon_nc(inifilename)
    ... #doctest: +ELLIPSIS
    Reading inifile .../mmr2dobsonunits.ini
    Getting model data for  3d  at  ...
    Saved to file: /.../20170702.nc

    """

    #Load inifile
    ini_dict = adaqcode.inifile.get_inidict(
        inifilename=inifilename)

    # Load the model data
    if 'models_list' not in ini_dict:
        ini_dict['models_list'] = [
            os.path.basename(path) if path is not None else None
            for path in ini_dict['models_dir_list']]
    md = adaqcode.adaq_functions.get_models(ini_dict, None,
                                            delete_gridded=False,
                                            surface_only=False)[0]
    # Extract the satellite pass time/s
    sat_pass_times = [int(t) for t in ini_dict['sat_pass_time_list']]
    md = md.extract(time=lambda cell: cell.point.hour
                    in sat_pass_times, gridded=True)

    #Get regular lat-lon cube
    llgrid = get_latlon_cube(lonmin=-12.0, lonmax=9.5,
                             latmin=44.0, latmax=62.0,
                             dlon=0.1, dlat=0.1)

    #Regrid
    md_ll = adaqcode.adaq_data.ADAQData()
    for cube in md.gridded_cube_list:
        llcube = cube.regrid(llgrid, iris.analysis.Linear())
        #Remove altitude - gives an error when saving to netcdf
        llcube.remove_aux_factory(llcube.aux_factories[0])
        md_ll.gridded_cube_list.append(llcube)

    #Save to file - one day per file
    ncdir = ini_dict['plot_dir'] + '/' + 'nc'
    if not os.path.isdir(ncdir):
        os.makedirs(ncdir)

    dt = ini_dict['start_datetime']
    while dt <= ini_dict['end_datetime']:
        tconstraint = iris.Constraint(
            time=lambda t: dt <= t.point < dt+datetime.timedelta(hours=24))
        cubes = md_ll.gridded_cube_list.extract(tconstraint)
        if cubes:
            date = dt.strftime('%Y%m%d')
            filename = ncdir + '/' + date + '.nc'
            iris.save(cubes, filename)
            print('Saved to file:', filename)

        dt += datetime.timedelta(hours=24)

    return ini_dict, md_ll

if __name__ == '__main__':

#    ini_dict, md_time_extr, du_time_ave_cubelist = mmr2dobsonunits()

    import doctest
    doctest.testmod()
