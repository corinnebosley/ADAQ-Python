Gallery 
========

.. contents:: Gallery Contents
    :local:

Line plots
----------

.. |Lineplot| image:: figures/lineplot.png
    :width: 400px

.. |timeseries| image:: figures/Harwell_O3.png
    :width: 400px
  
.. |Soccer| image:: figures/Soccer_Plot_O3.png
    :width: 400px
        
.. |Diurnal| image:: figures/adaq_plotting/Diurnal_O3.png
    :width: 400px
    
.. |timeseries_stats| image:: figures/adaq_plotting/Timeseries_of_bias_O3.png
    :width: 400px

.. |site_mean_timeseries| image:: figures/adaq_plotting/Timeseries_of_sites_nanmean_O3.png
    :width: 400px
    
.. |multi_species_timeseries| image:: figures/adaq_plotting/timeseries/Aberdeen_PM_Components_aqum_casestudy.png
    :width: 400px

.. |QQ| image:: figures/Quantile-Quantile_O3.png
    :width: 400px

+---------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------+
| |Lineplot|                                                                                        | |timeseries|                                                                                                         |
+---------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------+
| Example basic line plot - see :class:`line_plot.LinePlot`                                         | Example time-series - see :class:`timeseries_plot.TimeSeriesPlot`                                                    |
+---------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------+
| |Soccer|                                                                                          | |Diurnal|                                                                                                            |
+---------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------+
| Example soccer plot - see :class:`statistics_plotting.SoccerPlot`                                 | Example diurnal plot - see :func:`adaq_plotting.plot_diurnal`                                                        |
+---------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------+
| |timeseries_stats|                                                                                | |site_mean_timeseries|                                                                                               |
+---------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------+
| Example timeseries of statistics - see :func:`adaq_plotting.plot_timeseries_of_stats`             | Example of mean timeseries across all sites - see :func:`adaq_plotting.plot_timeseries_aggregate_stats`              |
+---------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------+
| |QQ|                                                                                              | |multi_species_timeseries|                                                                                           |
+---------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------+
| Example quantile-quantile plot - see :class:`statistics_plotting.QQPlot`                          | multiple species on same timeseries plot - see :func:`adaq_plotting.plot_tseries_multiple_snames`                    |
+---------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------+


Field plots
-----------

.. |Cs137_total_dep| image:: figures/Fieldplot_NAME_CAESIUM-137_TOTAL_DEPOSITION_Boundarylayer_201511291800.png
    :width: 400px
    
.. |VA_concen| image:: figures/Fieldplot_name_VOLCANIC_ASH_AIR_CONCENTRATION_FromFL025toFL050_201105230000.png
    :width: 700px
    
.. |veg_plot| image:: figures/Leafy_green_veg_plot.png
    :width: 400px
    
.. |o3_basic_NAME| image:: figures/adaq_plotting/gridded_fields/Fieldplot_aqum_oper_O3_201404020000_basicNAME.png
    :width: 400px
    
.. |o3_basic_aq| image:: figures/adaq_plotting/gridded_fields/Fieldplot_aqum_oper_O3_201404020000_basicAQ.png
    :width: 500px
    
.. |o3_full_aq| image:: figures/adaq_plotting/gridded_fields/Fieldplot_aqum_oper_O3_201404020000_fullAQ.png
    :width: 300px
    
.. |site_locs| image:: figures/site_locations.png
    :width: 300px
    
.. |lim_site_locs| image:: figures/limited_site_locations.png
    :width: 300px
    
.. |cross_section| image:: figures/adaq_plotting/CrossSection_aqum_T_201707020600.png
    :width: 500px
    :height: 250px
    
.. |cross_section_waypoints| image:: figures/adaq_plotting/CrossSectionSamplePoints_aqum_oper.png
    :width: 300px

+--------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------+
| |Cs137_total_dep|                                                                                      | |VA_concen|                                                                                                          |
+--------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------+
| Example field plot - see :class:`field_plot.FieldPlot`                                                 | Example field plot with annotation - see :class:`field_plot.FieldPlot`                                               |
+--------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------+
| |veg_plot|                                                                                             | |o3_basic_NAME|                                                                                                      |
+--------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------+
| Example veg plot see :class:`field_plot.FieldPlot`                                                     | Example field plot from :func:`adaq_plotting.plot_md_gridded_fields` using defaults='NAME' option.                   |
+--------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------+
| |o3_basic_aq|                                                                                          | |o3_full_aq|                                                                                                         |
+--------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------+
| Example field plot from :func:`adaq_plotting.plot_md_gridded_fields` using defaults='AQ' option.       | Example field plot from :func:`adaq_plotting.plot_md_gridded_fields` using defaults='AQ' option with extra settings. |
+--------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------+
| |site_locs|                                                                                            | |lim_site_locs|                                                                                                      |
+--------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------+
| Map of site locations from SitesInfo class - see :func:`sites_info.SitesInfo.plot_location_map`        | Map of site locations limited to 3 per country from CAMSAQObsData class - see :func:`camsaqobs_data.limit_sites`     |
+--------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------+
| |cross_section|                                                                                        | |cross_section_waypoints|                                                                                            |
+--------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------+
| Example cross section along given waypoints - see                                                      | Location of waypoints used in cross section - see :func:`adaq_vertical_plotting.plot_md_cross_section`               |
| :func:`adaq_vertical_plotting.plot_cube_cross_section`                                                 | or :func:`adaq_vertical_plotting.plot_section_sample_locations`                                                      |
| or :func:`adaq_vertical_plotting.plot_md_cross_section`                                                |                                                                                                                      |
+--------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------+



RSMC plots
^^^^^^^^^^

.. |rsmc_air_conc| image:: figures/adaq_plotting/rsmc_plots/Fieldplot_name_CAESIUM-137_DOSAGE_From0_0to500_0magl_201702261200_COL.png
    :width: 500px

.. |rsmc_dep| image:: figures/adaq_plotting/rsmc_plots/Fieldplot_name_CAESIUM-137_TOTAL_DEPOSITION_Boundarylayer_201702261200_COL.png
    :width: 500px

.. |rsmc_grey| image:: figures/adaq_plotting/rsmc_plots/Fieldplot_name_CAESIUM-137_TOTAL_DEPOSITION_Boundarylayer_201702261200_BW.png
    :width: 500px

.. |rsmc_toa| image:: figures/adaq_plotting/rsmc_plots/Fieldplot_name_Time_of_arrival_of_CAESIUM-137_02.png
    :width: 500px


+------------------------------------------+-----------------------------------+
| |rsmc_air_conc|                          | |rsmc_dep|                        |
|                                          |                                   |  
| RSMC color air concentration plot        | RSMC color deposition plot        |
| see :mod:`rsmc_plot`                     | see :mod:`rsmc_plot`              |
+------------------------------------------+-----------------------------------+
| |rsmc_grey|                              | |rsmc_toa|                        |
|                                          |                                   |
| RSMC greyscale deposition plot           | RSMC time of arrival plot         |
| see :mod:`rsmc_plot`                     | see :mod:`rsmc_plot`              |
+------------------------------------------+-----------------------------------+

Ensemble Plots
^^^^^^^^^^^^^^

.. |postage_stamp_air_conc| image:: figures/adaq_plotting/ensemble/Fieldplot_mogreps-g_name_CAESIUM-137_AIR_CONCENTRATION_Boundarylayeraverage_201908290200.png
    :width: 500px
    
.. |ens_exceedance_plot| image:: figures/adaq_plotting/ensemble/Fieldplot_mogreps-g_name_MEM_EXCEED_VOLCANIC_ASH_AIR_CONCENTRATION_OF_200_FromFL100toFL125_201905191600.png
    :width: 500px

.. |ens_percentile_plot| image:: figures/adaq_plotting/ensemble/Fieldplot_mogreps-g_name_95TH_PERCENTILE_OF_VOLCANIC_ASH_AIR_CONCENTRATION_FromFL100toFL125_201905191600.png
    :width: 500px

+---------------------------------------------------------------------+------------------------------------------------------------+
| |postage_stamp_air_conc|                                            | |ens_exceedance_plot|                                      |
+---------------------------------------------------------------------+------------------------------------------------------------+
| Example postage stamp plot - see :mod:`name_postage_stamp_plot`     | Example exceedance plot - see :mod:`ens_exceedance_plot`   |
+---------------------------------------------------------------------+------------------------------------------------------------+
| |ens_percentile_plot|                                               |                                                            |
+---------------------------------------------------------------------+------------------------------------------------------------+
| Example percentile plot - see :class:`ens_percentile_plot`          |                                                            |
+---------------------------------------------------------------------+------------------------------------------------------------+


Trajectory Plots
----------------

.. |TrajectoryPlot| image:: figures/adaq_plotting/TrajectoryPlot.png
    :width: 300px

+----------------------------+
| |TrajectoryPlot|           |
+----------------------------+
| Forward Trajectory Plot    |
+----------------------------+


Comparison Plots
----------------

.. |Montage| image:: figures/adaq_plotting/montage/Fieldplot_global-name_NO2_EULERIAN_CONCENTRATION_From0_0to100_0magl_201506140000.png
    :width: 1000px

+---------------------------------------------------------------------+
| |Montage|                                                           |
+---------------------------------------------------------------------+
| Example Montage plot - see :class:`name_field_plot_diff`            |
+---------------------------------------------------------------------+


Colour Maps
-----------

.. |YlGnBu| image:: figures/cmap.png
    :width: 300px
        
.. |daqi| image:: figures/daqi_cmap.png
    :width: 300px
    
.. |cams| image:: figures/cams_cmap.png
    :width: 300px

.. |rsmcdep| image:: figures/rsmcdep_cmap.png
    :width: 300px

.. |rsmcac| image:: figures/rsmcac_cmap.png
    :width: 300px

.. |rsmctoa| image:: figures/rsmctoa_cmap.png
    :width: 300px

.. |rsmcgrey| image:: figures/rsmcgrey_cmap.png
    :width: 300px
    
+------------------------------------------+-----------------------------------+----------------------------------+
| |YlGnBu|                                 | |daqi|                            | |cams|                           |
|                                          |                                   |                                  |  
| colormap 'YlGnBu' -                      | DAQI colormap -                   | CAMS colormap -                  |
| see :func:`plotting_functions.plot_cmap` | see :func:`aq_indices.daqi_cmap`  | see :func:`aq_indices.cams_cmap` |
+------------------------------------------+-----------------------------------+----------------------------------+
| |rsmcdep|                                | |rsmcac|                          |                                  |
|                                          |                                   |                                  |
| RSMC deposition colormap                 | RSMC air concentration colormap   |                                  |
| to match RSMC Toulouse -                 | to match RSMC Toulouse -          |                                  |
| see :func:`adaq_cmap.rsmcdep_cmap`       | see :func:`adaq_cmap.rsmcac_cmap` |                                  |
+------------------------------------------+-----------------------------------+----------------------------------+
| |rsmctoa|                                | |rsmcgrey|                        |                                  |
|                                          |                                   |                                  |
| RSMC time of arrival colormap -          | RSMC grey colormap for faxing -   |                                  |
| see :func:`adaq_cmap.rsmcdep_cmap`       | see :func:`adaq_cmap.rsmcac_cmap` |                                  |
+------------------------------------------+-----------------------------------+----------------------------------+

Histogram
---------

.. |Regime_bar| image:: figures/Regime_Bar.png
    :width: 400px
    
.. |Hist| image:: figures/Histogram_O3.png
    :width: 400px

+------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------+
| |Regime_bar|                                                                             | |Hist|                                                                                            |
+------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------+
| Weather regime bar chart - see :func:`weather_regime.plot_regime_bar`                    | Example histogram plot - see :class:`statistics_plotting.Histogram`                               |
+------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------+

Sandpit
-------

Tephi Plots
^^^^^^^^^^^

.. |Tephi_full| image:: figures/VAAC_Tephi_Barbs_plot_full.png
    :width: 400px

.. |Tephi_focus| image:: figures/VAAC_Tephi_Barbs_plot_focus.png
    :width: 400px

+------------------------------------------------------------------------------------------+---------------------------------------------------------------------+
| |Tephi_full|                                                                             | |Tephi_focus|                                                       | 
+------------------------------------------------------------------------------------------+---------------------------------------------------------------------+
| Example full tephigram - see :mod:`tephibarbs`                                           | Example zoomed tephigram - see :mod:`tephibarbs`                    |
+------------------------------------------------------------------------------------------+---------------------------------------------------------------------+

Trajectory Plots
^^^^^^^^^^^^^^^^

.. |trajplume| image:: figures/Trajectory_plume.png
    :width: 600px

+---------------------------------------------------------------------+
| |trajplume|                                                         |
+---------------------------------------------------------------------+
|  Example of trajectory plume plotting - see :mod:`trajectory_plume` |
+---------------------------------------------------------------------+

Ensemble Plots
^^^^^^^^^^^^^^^^

.. |wind_metogram| image:: figures/example_wind_metogram.png
    :width: 700px
   
.. |met_plume| image:: figures/example_met_plume.png
    :width: 400px

+------------------------------------------------------------------------------------------+---------------------------------------------------------------------+
| |met_plume|                                                                              | |wind_metogram|                                                     |
+------------------------------------------------------------------------------------------+---------------------------------------------------------------------+
| wind speed ensemble met plume - see :mod:`ensemble_met_plume`                            | wind direction metogram - see :mod:`wind_direction_metogram`        |
+------------------------------------------------------------------------------------------+---------------------------------------------------------------------+


Other
^^^^^

.. |windrose| image:: figures/StAthan_windrose_070714to130714.png
    :width: 400px

.. |pie_map| image:: figures/plot_pie_map_example.png
    :width: 400px

.. |rings| image:: figures/rings_example.png
    :width: 700px

.. |check_domain| image:: figures/check_domain_thetagrid.png
    :width: 400px

.. |cross_sec| image:: figures/SO2_cross_section.png
    :width: 400px

+------------------------------------------------------------------------------------------+---------------------------------------------------------------------+
| |windrose|                                                                               | |pie_map|                                                           |
+------------------------------------------------------------------------------------------+---------------------------------------------------------------------+
| Example windrose - see :mod:`windrose_plot`                                              | Example of pie charts plotted on a map - see :mod:`plot_pie_map`    |
+------------------------------------------------------------------------------------------+---------------------------------------------------------------------+
| |rings|                                                                                  | |check_domain|                                                      |
+------------------------------------------------------------------------------------------+---------------------------------------------------------------------+
| Example of distance rings - see :mod:`Rings`                                             | Example of domain checking - see :mod:`check_domain`                |
+------------------------------------------------------------------------------------------+---------------------------------------------------------------------+
| |cross_sec|                                                                              |                                                                     |
+------------------------------------------------------------------------------------------+---------------------------------------------------------------------+
| Plume cross section - see :mod:`plume_cross_section_plot`                                |                                                                     |
+------------------------------------------------------------------------------------------+---------------------------------------------------------------------+
