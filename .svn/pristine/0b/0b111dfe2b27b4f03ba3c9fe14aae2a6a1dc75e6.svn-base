Plotting sites cubes
====================

In this section, we will explain how to plot site-specific cubes, for example
as time-series or for statistical plots.

We begin by ensuring the adaqcode is imported. This is important both for
automated testing of this user guide page, but also if you want to follow
this page yourself.

>>> import sys
>>> adaq_path = '../../' #Full address of checked out code
>>> sys.path.append(adaq_path)
>>> import adaqcode

Now load some example sites data. Let us try using a full 5 day
set of data.

>>> ini_dict, sites_data, od, md_list = adaqcode.adaq_functions.get_exampledata(
... gridded_cube_list=False, sites_cube_list=True, exampletype='full') # doctest: +ELLIPSIS
Reading inifile .../adaqcode/example_data_5days.ini
Number of sites:  5

Plotting time-series
--------------------

A simple time-series can be produced using the
:class:`timeseries_plot.TimeSeriesPlot` class.

Firstly get a 1D cube with only a time-coordinate. For this
we will take the site Harwell and the short_name O3 from our sites_cube_list:

>>> print(md_list[0]) # doctest: +ELLIPSIS
<class '...adaq_data.ADAQData'> - Subclass of ADAQData - Contains:
sites_cube_list:
0: mass_concentration_of_nitrogen_dioxide_in_air / (ug/m3) (site_id: 5; time: 121)
1: mass_concentration_of_ozone_in_air / (ug/m3) (site_id: 5; time: 121)
2: mass_concentration_of_pm10_dry_aerosol_in_air / (ug/m3) (site_id: 5; time: 121)
3: mass_concentration_of_pm2p5_dry_aerosol_in_air / (ug/m3) (site_id: 5; time: 121)
4: mass_concentration_of_sulfur_dioxide_in_air / (ug/m3) (site_id: 5; time: 121)
gridded_cube_list:
< No cubes >
trajectory_cube_list:
< No cubes >

>>> har_o3 = md_list[0].extract(
... site_name='Harwell', short_name='O3', singlecube=True)
>>> print(har_o3.summary(shorten=True))
mass_concentration_of_ozone_in_air / (ug/m3) (time: 121)
>>> print(har_o3.attributes['label'])
aqum_oper

We could also get some observation data at this same site:

>>> har_o3_obs = od.extract(
... site_name='Harwell', short_name='O3', singlecube=True)
>>> print(har_o3_obs.summary(shorten=True))
mass_concentration_of_ozone_in_air / (ug/m3) (time: 121)
>>> print(har_o3_obs.attributes['label'])
Obs

Now start the time-series plot, by initialising the class:

>>> tsp = adaqcode.timeseries_plot.TimeSeriesPlot()

The cube can then be added as a line on this plot:

>>> tsp.add_line(har_o3)

The observation data could also be added as a line.
The :func:`line_plot.LinePlot.add_line` function also takes input about what
the line should look like, so we could ask for a black dashed line:

>>> tsp.add_line(har_o3_obs,
... colour='k', linestyle='--')

Other lines could also be added, modifying the label which will appear on the
plot:

>>> tsp.add_line(har_o3*0.5, label='mod*0.5')

The plot is then produced very simply:

>>> fig = tsp.plot()

The figure object is returned for further use or adjustment.
The plot could be saved, using tsp.save_fig(plotdir="outputdir/"),
shown to screen using plt.show(), or just closed:

>>> import matplotlib.pyplot as plt
>>> #plt.show()
>>> tsp.save_fig(adaqcode.CODE_DIR+'/adaqdocs/figures/user_guide/')  # doctest: +ELLIPSIS
Saved figure  .../user_guide/Harwell_O3.png

>>> plt.close()

.. image:: ../figures/Harwell_O3.png
    :scale: 50%

.. _user_guide_plothist:

Plotting a histogram
--------------------

A histogram of data can be easily plotted using
:class:`statistics_plotting.Histogram`.
Note other statisitical line plots can also be done in a similar way,
for example :class:`statistics_plotting.SoccerPlot`, and
:class:`statistics_plotting.QQPlot`.

Firstly we get some example cubes, in this case for O3,
for both model and observations:

>>> o3_obs = od.extract(short_name='O3', singlecube=True)
>>> o3_mod = md_list[0].extract(short_name='O3', singlecube=True)

Now initialise the Histogram class:

>>> hist = adaqcode.statistics_plotting.Histogram()

Add the cubes as individual lines:

>>> hist.add_line(o3_mod)
>>> hist.add_line(o3_obs, label='obs', colour='k')

The binsize used for the histogram can be manually adjusted:

>>> hist.binsize = 2.

The plot can then be easily produced:

>>> fig = hist.plot()

The figure object is returned for further use or adjustment.
The plot could be saved, using hist.save_fig(plotdir="outputdir/"),
shown to screen using plt.show(), or just closed:

>>> #plt.show()
>>> hist.save_fig(adaqcode.CODE_DIR+'/adaqdocs/figures/user_guide/')  # doctest: +ELLIPSIS
Saved figure  .../user_guide/Histogram_O3.png
>>> plt.close()

.. image:: ../figures/user_guide/Histogram_O3.png
    :scale: 50%


.. Note to developers. If you change any of the questions below then you will
   also need to change the answers in answers.rst

.. _plottingsites_questions_ref:

Plotting sites data on a map
----------------------------

To automatically plot all data from observations and model (or any list of sitescubes)
on a map (a different plot for each model, observation and time), we can use
:func:`adaq_plotting.plot_sitescube_maps`. 

The routine will either take both od and md_list as inputs, or you could give
it a list of site cubes of your own. Here we will give just a single cube, taken
from the first time of the observations sites cube extracted in the example above:

>>> sitescubelist = [o3_obs[:,0]]
>>> print(sitescubelist)
[<iris 'Cube' of mass_concentration_of_ozone_in_air / (ug/m3) (site_id: 5)>]

The function takes input of an ini_dict which contains lots of plotting options. 
This will generally be derived from a .ini formatted file. However it can also
be set up manually: 

>>> ini_dict = {}

There are various plotting options that are descripted within the documentation
for :func:`adaq_plotting.plot_sitescube_maps`, but here we will just pick a few
common ones:

Set the colour bar to be vertical:

>>> ini_dict['cbar_orientation'] = 'vertical'

Give an extent list to define the domain shown on the plot. If this is not set
it will zoom such that some sites appear right on the edges of the plot.
Note this is set up as [xmin,xmax,ymin,ymax].

>>> ini_dict['extent_list'] = [-12, 4.0, 48., 60.]

Set up a list of contour levels:

>>> import numpy as np
>>> ini_dict['levels_list'] = np.linspace(0., 70, 8)

Increase default marker size (from 20) as there are only a few sites to show:

>>> ini_dict['marker_size'] = 40

And finally ensure there is a directory for the plot output:

>>> ini_dict['plot_dir'] = adaqcode.CODE_DIR + "/adaqdocs/figures/user_guide/"

Now call the routine:

>>> adaqcode.adaq_plotting.plot_sitescube_maps(ini_dict, sitescubelist=sitescubelist,
... short_name='O3') # doctest: +ELLIPSIS
Plotting sitescube maps
Saved figure .../user_guide/Fieldplot_sites_Obs_O3_201403260000.png

.. image:: ../figures/user_guide/Fieldplot_sites_Obs_O3_201403260000.png
   :scale: 75%      

Note by using the keyword 'classifications_dict', different site types can be
plotted using different markers, for example rural sites with a circle and 
urban sites with a square, which are then shown in a legend. An example plot
of this is given below, but for instructions on using this see the documentation
for :func:`adaq_plotting.plot_sitescube_maps`. 

.. image:: ../figures/adaq_plotting/gridded_fields/Fieldplot_sites_aqum_casestudy_O3_201404020000.png
   :scale: 75%      


Plotting Sites Data Exercises
-----------------------------

  1. Produce a time series plot of Volcanic ash air concentration at Bristol.

    a. First load in the data in the name time series file *Time_series_grid1.txt*
       in the *name* directory in the sample data directory
    b. Then use extract to select the data from Bristol
    c. Set up a time series plot
    d. Add the Bristol cube as a line
    e. Plot the time series plot
    f. Now create a new time series plot object and plot the same data again
       but change the colour of the line and the title
       (Hint: to change the title use *tsp.title=*)

  2. Now use the example data to produce a quantile-quantile plot
     (see :class:`statistics_plotting.QQPlot`)

    a. First you will need to load in the full 5 day set of site example data as
       shown at the top of this page.
    b. Select a species to plot (this could be *O3* or *PM2p5*) and extract
       these from both the model and observations site cube. (Hint: the example
       data contains a model data list so you will need to choose which item in
       the list to extract the data from or you can extract both separately)
    c. Set up a quantile-quantile plot object (qq = adaqcode.statistics_plotting.QQPlot())
    d. Add the observations cube to the x-axis (qq.xcube = my_obs_cube)
    e. Add the model cube(s) using the add_line function
    f. Plot the quantile-quantile plot
    
  3. Use the standard short set of example data to produce site data maps for
     both observations for PM2.5.
     
     a. First load in the short set of example data (hint no keywords are 
        required for this) - you should get a 1 day example
     b. Set up a new ini_dict and just give a plot directory
     c. Produce site maps for all times. Hint: use the keyword 'od' rather 
        than using 'sitescubelist'. Also need to give short_name='PM2p5'
     d. Repeat the above, but using an appropriate domain 'extent_list' and
        'levels_list' and a different colour map ('cmap') and marker_size. You 
        to make a prettier looking plot. Note, setting 'filesuffix' will stop
        these plots overwritting the first ones.

Answers to :ref:`plottingsites_answers_ref`
