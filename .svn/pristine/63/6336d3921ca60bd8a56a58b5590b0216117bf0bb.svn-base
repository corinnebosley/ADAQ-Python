Answers to Exercises
====================

..
   Note to developers. If you change any of the questions in this document
   you will need to change the questions in the relevant section of the
   user guide.

.. _intro_answers_ref:

Introduction
------------
(:ref:`intro_questions_ref`)


  1. Print the short name and the label of the 5th gridded cube in the
     sample data file called *aqum_casestudy_5days.nc*

     a. Load the data file *aqum_casestudy_5days.nc* which is in the
        *gridded_cube_list* directory in the sample data directory.
        (Hint: the sample data directory is *adaqcode.SAMPLE_DATADIR*)
        **Note: We begin by ensuring the adaqcode is imported. You will need to
        do this anytime you exit python and restart the exercises but for brevity
        we have only included it in Exercise 1, Question 1.**

        >>> import sys
        >>> adaq_path='../../' #Full address of checked out code
        >>> sys.path.append(adaq_path)
        >>> import adaqcode

        >>> md = adaqcode.adaq_data.ADAQData()
        >>> md.load_gridded(adaqcode.SAMPLE_DATADIR+'gridded_cube_list/aqum_casestudy_5days.nc')
        [<iris 'Cube' of mass_concentration_of_nitrogen_dioxide_in_air / (ug/m3) (time: 121; grid_latitude: 182; grid_longitude: 146)>,
        <iris 'Cube' of mass_concentration_of_ozone_in_air / (ug/m3) (time: 121; grid_latitude: 182; grid_longitude: 146)>,
        <iris 'Cube' of mass_concentration_of_pm10_dry_aerosol_in_air / (ug/m3) (time: 121; grid_latitude: 182; grid_longitude: 146)>,
        <iris 'Cube' of mass_concentration_of_pm2p5_dry_aerosol_in_air / (ug/m3) (time: 121; grid_latitude: 182; grid_longitude: 146)>,
        <iris 'Cube' of mass_concentration_of_sulfur_dioxide_in_air / (ug/m3) (time: 121; grid_latitude: 182; grid_longitude: 146)>]

     b. Select the 5th gridded cube and print the short name and the label

        >>> gcl = md.gridded_cube_list
        >>> print(gcl[4].attributes['short_name'])
        SO2
        >>> print(gcl[4].attributes['label'])
        aqum_casestudy

  2. Now load in the sites_cube_list *aqum_oper_1days.nc* from the
     *sites_cube_list* directory in the sample data directory and print
     the abbreviated names of all the sites in the first sites cube list.
     (Hint: The coordinate of abbreviated names is called *abbrev*).

     >>> md_scl = md.load_ts(adaqcode.SAMPLE_DATADIR+'sites_cube_list/aqum_oper_1days.nc')
     >>> scl = md.sites_cube_list
     >>> sites_cube = scl[0]
     >>> print(sites_cube.coord('abbrev').points)
     ['YW' 'ACTH' 'AH' 'ABD' 'HAR']

     Note that to get only the site abbreviations it is necessary to print
     the points attribute of the coordinate not the whole coordinate.

  3. Read in the sites files *ctbto_sites.txt* in the *ctbto_obs* subfolder
     in the sample data directory and print out the sites data.

     >>> sitesfile = adaqcode.SAMPLE_DATADIR+'ctbto_obs/ctbto_sites.txt'
     >>> si = adaqcode.sites_info.SitesInfo()
     >>> sites_data = si.read_from_file(sitesfile, allsites=True)
     Number of sites:  2
     >>> fmtstr = adaqcode.sites_info.format_string(sites_data)
     >>> for site in sites_data:
     ...     print(fmtstr.format(*site))
     JPP38,  139.080,  36.300,-999,  Takasaki,CTBTO
     USP70, -121.360,  38.670,-999,Sacramento,CTBTO


  4. Read in the inifile *name_field_plot.ini* which is in the *adaqscripts*
     directory in the ADAQ code. Print the *models_fmt_list* and the
     *field_attribute_dict*. (Hint: the code directory is *adaqcode.CODE_DIR*)

     >>> inifilename = adaqcode.CODE_DIR + 'adaqscripts/name_field_plot.ini'
     >>> ini_dict = adaqcode.inifile.Inifile(inifilename) # doctest: +ELLIPSIS
     Reading inifile .../adaqscripts/name_field_plot.ini

     or

     >>> ini_dict = adaqcode.inifile.get_inidict(inifilename) # doctest: +ELLIPSIS
     Reading inifile .../adaqscripts/name_field_plot.ini

     then

     >>> print(ini_dict['models_fmt_list'])
     ['name']
     >>> print(ini_dict['field_attribute_dict'])
     {'Species': 'VOLCANIC_ASH'}

.. _loadingdata_answers_ref:

Loading Data
------------
(:ref:`loadingdata_questions_ref`)

  1. Load in data with the short name *PM10*, a validity date
     between 0Z 01/04/2014 and 0Z 02/04/2014 from the first full day of
     forecast. This data can be found in the *aqum_output/oper* directory
     in the sample data directory and should be given the label *Oper*
     (short for operational).

     a. Create a short_name, start_datetime and end_datetime variable to
        limit the data loading with.

        >>> short_name_list = ['PM10']
        >>> import datetime
        >>> start_datetime = datetime.datetime(2014,4,1,0)
        >>> end_datetime = datetime.datetime(2014,4,2,0)

     b. Setup a AQUMppfiles object and use it obtain the filenames for day
        1 forecasts.

        >>> directory = adaqcode.SAMPLE_DATADIR+'aqum_output/oper'
        >>> ppfiles = adaqcode.pp_data.AQUMppFiles()
        >>> filenames = ppfiles.get_filenames(directory,
        ...    start_datetime = start_datetime, end_datetime = end_datetime,
        ...    forecast_day = 1)
        >>> print(filenames) #  doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
        ['.../oper/aqum_20140330_T+024_QM18.pp',
        '.../oper/aqum_20140331_T+000_QM18.pp',
        '.../oper/aqum_20140331_T+012_QM18.pp',
        '.../oper/aqum_20140331_T+024_QM18.pp',
         '.../oper/aqum_20140401_T+000_QM18.pp']

     c. Set up a PPData object and pass all the variables (e.g. short_name,
        filenames) to it to load the data for day 1 forecasts.

        >>> md = adaqcode.pp_data.PPData()
        >>> md.readdata(filenames = filenames,
        ...    short_name_list = short_name_list,
        ...    start_datetime = start_datetime,
        ...    end_datetime = end_datetime,
        ...    forecast_day = 1,
        ...    label='Oper')
        [<iris 'Cube' of mass_concentration_of_pm10_dry_aerosol_in_air / (ug m-3) (time: 25; grid_latitude: 182; grid_longitude: 146)>]

  2. Load in the XENON-133 data from the Fields_grid1*.txt files stored
     in the *name* directory in the sample data directory. Then using the sites
     data stored in the *ctbto_sites.txt* file in the *ctbto_obs* directory in
     the sample data directory extract the Xenon-133 data at Sacramento
     and Takasaki.

     a. Setup the directory and some field_attributes to restrict loading

        >>> directory = adaqcode.SAMPLE_DATADIR+'name/'
        >>> field_attributes = {'Species':'XENON-133'}

     b. Setup a NAMEData object and use it to read some files in

        >>> name = adaqcode.name_data.NAMEData()
        >>> filenames = [directory + 'Fields_grid1*.txt']
        >>> gcl = name.readdata(filenames, field_attributes=field_attributes)
        >>> print(gcl)
        0: XENON-133_AIR_CONCENTRATION / (Bq / m^3) (time: 9; latitude: 90; longitude: 180)

     c. Setup a SitesInfo object and pass the sites file location to it to read

        >>> si = adaqcode.sites_info.SitesInfo()
        >>> sitesfile = adaqcode.SAMPLE_DATADIR+'ctbto_obs/ctbto_sites.txt'
        >>> sites_data = si.read_from_file(sitesfile, allsites=True)
        Number of sites:  2

     d. Use the sites_data to extract the sites_cubes from the NAMEData object

        >>> name_scl = name.extract_sites(sites_data)
        >>> print(name) # doctest: +ELLIPSIS
        <class '...name_data.NAMEData'> - Subclass of ADAQData - Contains:
        sites_cube_list:
        0: XENON-133_AIR_CONCENTRATION / (Bq / m^3) (site_id: 2; time: 9)
        gridded_cube_list:
        0: XENON-133_AIR_CONCENTRATION / (Bq / m^3) (time: 9; latitude: 90; longitude: 180)
        trajectory_cube_list:
        < No cubes >

  3. Use the information in the inifile example_data_5days.ini to load some
     observations data.

     a. Read the information in adaqcode/example_data_5days.ini into an ini_dict

        >>> inifilename = adaqcode.CODE_DIR + 'adaqcode/example_data_5days.ini'
        >>> ini_dict = adaqcode.inifile.get_inidict(inifilename) # doctest: +ELLIPSIS
        Reading inifile .../adaqcode/example_data_5days.ini

     b. Use the ini_dict to retrieve the sites information

        >>> sites_data = adaqcode.sites_info.get_siteinfo(ini_dict)
        Number of sites:  5

     c. Load in the observations

        >>> od = adaqcode.adaq_functions.get_obs(ini_dict, sites_data)  # doctest: +ELLIPSIS
        Creating observation data at  ...
        Reading obs data files
        Found obs for  .../AURN_obs/ABD_20140101_20140818.txt
        Found obs for  .../AURN_obs/ACTH_20140101_20140818.txt
        Found obs for  .../AURN_obs/AH_20140101_20140818.txt
        Found obs for  .../AURN_obs/HAR_20140101_20140818.txt
        Found obs for  .../AURN_obs/YW_20140101_20140818.txt
        Creating obs cubes

     d. Print the returned ADAQ object

        >>> print(od) # doctest: +ELLIPSIS
        <class '...aurn_data.AURNData'> - Subclass of ADAQData - Contains:
        sites_cube_list:
        0: mass_concentration_of_nitrogen_dioxide_in_air / (ug/m3) (site_id: 5; time: 121)
        1: mass_concentration_of_ozone_in_air / (ug/m3) (site_id: 5; time: 121)
        2: mass_concentration_of_pm10_ambient_aerosol_in_air / (ug/m3) (site_id: 5; time: 121)
        3: mass_concentration_of_pm2p5_ambient_aerosol_in_air / (ug/m3) (site_id: 5; time: 121)
        4: mass_concentration_of_sulfur_dioxide_in_air / (ug/m3) (site_id: 5; time: 121)
        gridded_cube_list:
        < No cubes >
        trajectory_cube_list:
        < No cubes >


  4. Use the information in the inifile name_field_plot.ini to load some NAME
     model data.

     a. Read the information in adaq_scripts/name_field_plot.ini into an ini_dict

        >>> inifilename = adaqcode.CODE_DIR + 'adaqscripts/name_field_plot.ini'
        >>> ini_dict = adaqcode.inifile.get_inidict(inifilename) # doctest: +ELLIPSIS
        Reading inifile .../adaqscripts/name_field_plot.ini

     b. Now use the ini_dict to load in some NAME model data

        >>> md_list = adaqcode.adaq_functions.get_models(ini_dict) # doctest: +ELLIPSIS
        Getting model data for  name  at ...

     c. Print out the contents of md_list[0]

        >>> print(md_list[0]) # doctest: +ELLIPSIS
        <class '...name_data.NAMEData'> - Subclass of ADAQData - Contains:
        sites_cube_list:
        < No cubes >
        gridded_cube_list:
        0: VOLCANIC_ASH_AIR_CONCENTRATION / (g/m3) (latitude: 480; longitude: 640)
        trajectory_cube_list:
        < No cubes >

.. _plottinggridded_answers_ref:

Plotting Gridded Data
---------------------
(:ref:`plottinggridded_questions_ref`)

  1. Plotting gridded fields automatically. Load in some NAME data using the inifile
     *name_field_plot.ini* in the adaqscripts directory and plot the field with the
     short name *VOLCANIC_ASH_AIR_CONCENTRATION* using the *NAME* defaults. You may
     also find it helpful to take a look at the inifile with a text editor or print
     out the contents of inifile so you can see what variables are set.

     a. Load in the NAME data using the inifile *name_field_plot.ini*. (Hint: You did
        this in question 4 in :ref:`loadingdata_questions_ref`)

        >>> inifilename = adaqcode.CODE_DIR + 'adaqscripts/name_field_plot.ini'
        >>> ini_dict = adaqcode.inifile.get_inidict(inifilename) # doctest: +ELLIPSIS
        Reading inifile .../adaqscripts/name_field_plot.ini
        >>> md_list = adaqcode.adaq_functions.get_models(ini_dict) # doctest: +ELLIPSIS
        Getting model data for  name  at  ...

     b. Change the location of the plot_dir to somewhere in your home space. Here the
        user_guide directory is only used for the purpose of the doctests.

        >>> ini_dict['plot_dir'] = adaqcode.CODE_DIR + "/adaqdocs/figures/user_guide/"

     c. Plot the data setting the short name to *VOLCANIC_ASH_AIR_CONCENTRATION* and
        using the *NAME* defaults for plotting.

        >>> adaqcode.adaq_plotting.plot_md_gridded_fields(ini_dict, md_list,
        ... "VOLCANIC_ASH_AIR_CONCENTRATION", defaults='NAME') # doctest: +ELLIPSIS
        Plotting gridded fields
        Saved figure  .../Fieldplot_name_VOLCANIC_ASH_AIR_CONCENTRATION_FromFL025toFL050_201105230000.png

  2. If you look at the inifile *name_field_plot.ini* you will see that only the extent
     and the levels have been set. Try replotting the *VOLCANIC_ASH_AIR_CONCENTRATION*
     with a different colourmap (*cmap*) and changing the extent to focus on Iceland.

     a. Change the colourmap - some alternative colourmaps can be seen here:
        http://matplotlib.org/1.3.1/examples/color/colormaps_reference.html

        >>> ini_dict['cmap'] = 'YlOrBr'

     b. Change the extent_list. Note that the numbers for the extent list can be provided
        as strings or floats. They are converted into floats (if necessary) in
        *plot_md_gridded_fields*

        >>> ini_dict['extent_list'] = [-30, -10, 60, 70]


     c. Re-run *plot_md_gridded_fields*. Don't forget to add a *filesuffix* so you
        don't overwrite your first plot.

        >>> adaqcode.adaq_plotting.plot_md_gridded_fields(ini_dict, md_list, "VOLCANIC_ASH_AIR_CONCENTRATION",
        ... defaults='NAME', filesuffix='_zoom') # doctest: +ELLIPSIS
        Plotting gridded fields
        Saved figure  .../Fieldplot_name_VOLCANIC_ASH_AIR_CONCENTRATION_FromFL025toFL050_201105230000_zoom.png

.. _plottingsites_answers_ref:

Plotting Sites Data
-------------------
(:ref:`plottingsites_questions_ref`)

  1. Produce a time series plot of Volcanic ash air concentration at Bristol
     using the NAME time series file loaded in the examples in section 2.

    a. First load in the data in the name time series file in the *name*
       directory in the sample data directory

       >>> directory = adaqcode.SAMPLE_DATADIR+'name/'
       >>> filenames = [directory + 'Time_series_grid1.txt']
       >>> name_ts = adaqcode.name_data.NAMEData()
       >>> scl = name_ts.readdata(filenames)

    b. Then use extract to select the data from Bristol

       >>> scl_cube = name_ts.extract(site_name='Bristol', singlecube=True)

    c. Set up a time series plot

       >>> tsp = adaqcode.timeseries_plot.TimeSeriesPlot()

    d. Add the Bristol cube as a line

       >>> tsp.add_line(scl_cube)

    e. Plot the time series plot

       >>> fig = tsp.plot()
       >>> import matplotlib.pyplot as plt
       >>> #plt.show()
       >>> plt.close()

    f. Now create a new time series plot object and plot the same data again
       but change the colour of the line and the title
       (Hint: to change the title use *tsp.title=*)

       >>> tsp = adaqcode.timeseries_plot.TimeSeriesPlot()
       >>> tsp.add_line(scl_cube, colour='b')
       >>> tsp.title = "Air Concentration of Volcanic Ash at Bristol"
       >>> fig = tsp.plot()

       >>> import matplotlib.pyplot as plt
       >>> plt.close()

  2. Now use the example data to produce a quantile-quantile plot
     (see :class:`statistics_plotting.QQPlot`)

    a. First you will need to load in the full 5 day set of site example data as
       shown at the top of this page.

       >>> ini_dict, sites_data, od, md_list = adaqcode.adaq_functions.get_exampledata(
       ... gridded_cube_list=False, sites_cube_list=True) # doctest: +ELLIPSIS
       Reading inifile .../adaqcode/example_data_1days.ini
       Number of sites:  5

    b. Select a species to plot (this could be *O3* or *PM2p5*) and extract
       these from both the model and observations site cube. (Hint: the example
       data contains a model data list so you will need to choose which item in
       the list to extract the data from or you can extract both separately)

       >>> obs_cube = od.extract(short_name='PM2p5', singlecube=True)
       >>> mod_cube = md_list[0].extract(short_name='PM2p5', singlecube=True)
       >>> mod_cube2 = md_list[1].extract(short_name='PM2p5', singlecube=True)

    c. Set up a quantile-quantile plot object (qq = adaqcode.statistics_plotting.QQplot())

       >>> qq = adaqcode.statistics_plotting.QQPlot()

    d. Add the observations cube to the x-axis (qq.xcube = my_obs_cube)

       >>> qq.xcube = obs_cube

    e. Add the model cube(s) using the add_line function

       >>> qq.add_line(mod_cube)
       >>> qq.add_line(mod_cube2, marker='+')

    f. Plot the quantile-quantile plot

       >>> fig = qq.plot()
       >>> import matplotlib.pyplot as plt
       >>> #plt.show()
       >>> plt.close()
       
       
  3. Use the standard short set of example data to produce site data maps for
     both observations for PM2.5.
     
     a. First load in the short set of example data (hint no keywords are 
        required for this) - you should get a 1 day example:
        
        >>> ini_dict, sites_data, od, md_list = adaqcode.adaq_functions.get_exampledata() 
        ... # doctest: +ELLIPSIS
        Reading inifile .../example_data_1days.ini
        Number of sites:  5
     
     b. Set up a new ini_dict and just give a plot directory . Here the
        user_guide directory is only used for the purpose of the doctests. 
      
      	>>> ini_dict = {}
        >>> ini_dict['plot_dir'] = adaqcode.CODE_DIR + "/adaqdocs/figures/user_guide/"
        
     c. Produce site maps for all times. Hint: use the keyword 'od' rather 
        than using 'sitescubelist'. Also need to give short_name='PM2p5'
        
        >>> adaqcode.adaq_plotting.plot_sitescube_maps(ini_dict, od=od, short_name='PM2p5') 
        ... # doctest: +ELLIPSIS
        Plotting sitescube maps
        Saved figure  ...user_guide/Fieldplot_sites_Obs_PM2p5_201404020000.png
	...
        Saved figure  .../user_guide/Fieldplot_sites_Obs_PM2p5_201404030000.png
     
     d. Repeat the above, but using an appropriate domain 'extent_list' and
        'levels_list' and a different colour map ('cmap') and marker_size. You 
        to make a prettier looking plot. Note, setting 'filesuffix' will stop
        these plots overwritting the first ones.
        
        >>> ini_dict['extent_list'] =  [-12, 4.0, 48., 60.]
        >>> ini_dict['levels_list'] = [0., 10., 20., 30., 40., 50., 60.]
        >>> ini_dict['cmap'] = 'jet'
        >>> ini_dict['marker_size'] = 50
        >>> adaqcode.adaq_plotting.plot_sitescube_maps(ini_dict, od=od, short_name='PM2p5',
        ... filesuffix='_pretty') # doctest: +ELLIPSIS
        Plotting sitescube maps
        Saved figure  .../user_guide/Fieldplot_sites_Obs_PM2p5_201404020000_pretty.png
	...
        Saved figure  .../user_guide/Fieldplot_sites_Obs_PM2p5_201404030000_pretty.png
        
       

.. _comparing_answers_ref:

Comparing Model and Observations
--------------------------------
(:ref:`comparing_questions_ref`)

  1. Using the example data extract the data with short_name *PM2p5* at Aberdeen
     compute the bias and the correlation.

    a. Load the example data

       >>> ini_dict, sites_data, od, md_list = adaqcode.adaq_functions.get_exampledata()
       ... # doctest: +ELLIPSIS
       Reading inifile .../adaqcode/example_data_1days.ini
       Number of sites:  5

    b. Run *match_data_for_stats* to ensure that the data has the same
       points/times

       >>> od_stats, md_stats_list = adaqcode.adaq_functions.match_data_for_stats(ini_dict, od, md_list)

    c. Extract the observations for PM2p5 at Aberdeen

       >>> pm2p5_obs = od_stats.extract(short_name='PM2p5', abbrev='ABD', singlecube=True)

    d. Extract some model data for PM2p5 at Aberdeen

       >>> pm2p5_mod = md_stats_list[0].extract(short_name='PM2p5', abbrev='ABD',
       ... singlecube=True)

    e. Initialise the statistics using the extracted observations and model data

       >>> stats = adaqcode.timeseries_stats.TimeSeriesStats(pm2p5_obs, pm2p5_mod)

    f. Compute and print the bias and the correlation

       >>> bias = stats.bias()
       >>> print('{:.3f}'.format(bias))
       -6.244

       >>> correlation = stats.correlation()
       >>> print('{:.3f}'.format(correlation))
       -0.216

  2. Following the method shown in :ref:`plotting_model_obs_ref` reproduce the
     quantile-quantile plot you produced in Exercise 2 in
     :ref:`plottingsites_questions_ref`. (See :func:`adaq_plotting.plot_qq`)

     a. Load in the example data

        >>> ini_dict, sites_data, od, md_list = adaqcode.adaq_functions.get_exampledata()
        ... # doctest: +ELLIPSIS
        Reading inifile .../adaqcode/example_data_1days.ini
        Number of sites:  5

     b. Change the location of the plotting directory 'plot_dir'. Here the
        user_guide directory is only used for the purpose of the doctests. You
        should be setting the plot_dir to a location in your own home or
        $DATADIR space.

        >>> ini_dict['plot_dir'] = adaqcode.CODE_DIR+'/adaqdocs/figures/user_guide/'

     c. Produce the quantile-quantile plot
        (Hint: the function is called plot_qq)

        >>> dp = adaqcode.adaq_plotting.plot_qq(ini_dict, od, md_list, 'PM2p5')
        ... # doctest: +ELLIPSIS
        Plotting quantile-quantile plot
        Saved figure  .../Quantile-Quantile_PM2p5.png

