Plotting gridded data
=====================

In this section, we will explain how to plot gridded fields.
Note for simple plotting it may be easier to use one of the scripts which
are explained in :ref:`adaqscripts`, for example
:mod:`name_field_plot`.py or :mod:`aq_plot`.py.

We begin by ensuring the adaqcode is imported. This is important both for automated
testing of this user guide page, but also if you want to follow this page yourself.

>>> import sys
>>> adaq_path = '../../' #Full address of checked out code
>>> sys.path.append(adaq_path)
>>> import adaqcode

Now load some example gridded data. This could be done using
:func:`adaq_functions.get_exampledata`:

>>> ini_dict, sites_data, od, md_list = adaqcode.adaq_functions.get_exampledata(
... gridded_cube_list=True, sites_cube_list=False) # doctest: +ELLIPSIS
Reading inifile .../adaqcode/example_data_1days.ini

.. topic:: Alternatively (eg in real code!)

   Or alternatively, just read in an ini_dict and use this to load md_list,
   using :func:`inifile.get_inidict` and :func:`adaq_functions.get_models`
   (introduced in :ref:`user_guide_inifiles` and :ref:`loading_using_inifiles`
   respectively).

   >>> inifilename = adaqcode.CODE_DIR + 'adaqcode/example_data_1days.ini'
   >>> ini_dict = adaqcode.inifile.get_inidict(inifilename)  # doctest: +ELLIPSIS
   Reading inifile .../adaqcode/example_data_1days.ini
   >>> md_list = adaqcode.adaq_functions.get_models(ini_dict) # doctest: +ELLIPSIS
   Getting model data for  aqum_oper  at  ...
   Getting model data for  aqum_casestudy  ...

   Using this method you may also want to convert to ug/m3 - a shortcut
   function for this is :func:`adaq_functions.unit_conversion`
   which also expects observations, but this could just be set to None:

   >>> od, md_list = adaqcode.adaq_functions.unit_conversion(None, md_list, chem_units='ug/m3', aerosol_units='ug/m3')

Now check the data:

>>> for md in md_list:
...     print(md) # doctest: +ELLIPSIS
<class '...'> - Subclass of ADAQData - Contains:
sites_cube_list:
< No cubes >
gridded_cube_list:
0: mass_concentration_of_pm2p5_dry_aerosol_in_air / (ug/m3) (time: 25; grid_latitude: 182; grid_longitude: 146)
1: mass_concentration_of_ozone_in_air / (ug/m3) (time: 25; grid_latitude: 182; grid_longitude: 146)
trajectory_cube_list:
< No cubes >
<class '...'> - Subclass of ADAQData - Contains:
sites_cube_list:
< No cubes >
gridded_cube_list:
0: mass_concentration_of_pm2p5_dry_aerosol_in_air / (ug/m3) (time: 25; grid_latitude: 182; grid_longitude: 146)
1: mass_concentration_of_ozone_in_air / (ug/m3) (time: 25; grid_latitude: 182; grid_longitude: 146)
trajectory_cube_list:
< No cubes >

md_list could contain any number of models and the routines described below
will still work, plotting fields individually for all the times in
the data. However, for speed, we will work with the first model in
the list. We will also shorten the total fields used so only a single
time is plotted.

>>> md_list = [md_list[0]]
>>> for i, cube in enumerate(md_list[0].gridded_cube_list):
...     md_list[0].gridded_cube_list[i] = cube[0]
>>> for md in md_list:
...    print(md.gridded_cube_list)
0: mass_concentration_of_pm2p5_dry_aerosol_in_air / (ug/m3) (grid_latitude: 182; grid_longitude: 146)
1: mass_concentration_of_ozone_in_air / (ug/m3) (grid_latitude: 182; grid_longitude: 146)

Plotting gridded fields automatically
-------------------------------------

To do a simple set of plotting, then we can use
:func:`adaq_plotting.plot_md_gridded_fields`. This uses various default option
settings depending on the type of plotting required. It also takes an ini_dict
which contains lots of plotting options. This will generally be derived from
a .ini formatted file. However for the simplest case, it
only requires a plot_dir to state where the plots will be saved:

>>> ini_dict = {}
>>> ini_dict['plot_dir'] = adaqcode.CODE_DIR + "/adaqdocs/figures/user_guide/"

We also need to give the routine a short_name of the field to plot -
here we choose 'PM2p5'.

For example for a default AQ-type plot, which uses the native coordinate
system of the cubes and a linear colour scale.

>>> adaqcode.adaq_plotting.plot_md_gridded_fields(ini_dict, md_list, 'PM2p5',
... defaults='AQ', filesuffix='_basicAQ') # doctest: +ELLIPSIS
Plotting gridded fields
Saved figure  .../Fieldplot_aqum_oper_PM2p5_201404020000_basicAQ.png

.. image:: ../figures/user_guide/Fieldplot_aqum_oper_PM2p5_201404020000_basicAQ.png
   :scale: 75%

We can also plot the same data using NAME-based defaults, which instead use a
log-based colour scale and use the PlateCarree projection. To make room for the
longer numbers on the colorbar, we can also set the figure size using 'figsize'.

>>> ini_dict['figsize'] = [9.0, 6.0]
>>> adaqcode.adaq_plotting.plot_md_gridded_fields(ini_dict, md_list, 'PM2p5',
... defaults='NAME', filesuffix='_basicNAME') # doctest: +ELLIPSIS
Plotting gridded fields
Saved figure  .../Fieldplot_aqum_oper_PM2p5_201404020000_basicNAME.png

.. image:: ../figures/user_guide/Fieldplot_aqum_oper_PM2p5_201404020000_basicNAME.png
   :scale: 75%

Now try also adding some extra options in, such as a list of contour levels etc.
The contour levels are set using 'levels_list' and the colours are set by setting
the cmap to the daily air quality index colours (generated by :func:`aq_indices.daqi_cmap`).
(Note the syntax ini_dict.get(a,b) here gets the value a from the dictionary, but if
it is not in the dictionary, sets it to b by default).
We also modify the colorbar so that it appears to the right of the plot instead
of underneath by setting cbar_orientation to vertical.
A full list of the options available which can be used in the ini_dict can be
seen in the parameters for :func:`adaq_plotting.plot_md_gridded_fields`.
For an example of some of these in use see in adaqcode/inifile.ini.

>>> import aq_indices
>>> import numpy as np
>>> ini_dict['figsize'] = [6.0, 6.0]
>>> ini_dict['cbar_orientation'] = 'vertical'
>>> ini_dict['levels_list'] = np.linspace(0, 80., 11)
>>> ini_dict['cmap'] = ini_dict.get('cmap', adaqcode.aq_indices.daqi_cmap())
>>> adaqcode.adaq_plotting.plot_md_gridded_fields(ini_dict, md_list, 'PM2p5',
... defaults='AQ', filesuffix='_fullAQ') # doctest: +ELLIPSIS
Plotting gridded fields
Saved figure  .../Fieldplot_aqum_oper_PM2p5_201404020000_fullAQ.png

.. image:: ../figures/user_guide/Fieldplot_aqum_oper_PM2p5_201404020000_fullAQ.png
   :scale: 75%


Plotting gridded fields manually
--------------------------------

Fields can also be plotted manually for further control. Here we will try
working directly with NAME data, but the same can be applied for any cubes.
(Note for more information on values for field_attributes,
see :ref:`field_attributes_ref`)

>>> sample_data_path = adaqcode.SAMPLE_DATADIR+'name/'
>>> name = adaqcode.name_data.NAMEData()
>>> name.readdata(sample_data_path + 'Fields_grid3_201511291800.txt',
... field_attributes = {'Quantity': 'Total deposition'})
[<iris 'Cube' of CAESIUM-137_TOTAL_DEPOSITION / (Bq/m2) (latitude: 205; longitude: 245)>]
>>> cube = name.gridded_cube_list[0]
>>> print(cube)
CAESIUM-137_TOTAL_DEPOSITION / (Bq/m2) (latitude: 205; longitude: 245)
     Dimension coordinates:
          latitude                              x               -
          longitude                             -               x
     Scalar coordinates:
          source_latitude: 54.0303 degrees
          source_longitude: -2.9175 degrees
          time: 2015-11-29 18:00:00, bound=(2015-11-27 18:00:00, 2015-11-29 18:00:00)
          z: Boundary layer
     Attributes:
          End of release: 0000UTC 28/11/2015
          Forecast duration: 48 hours
          Met data: NWP Flow.UKV_PT1_flow; NWP Flow.UKV_PT2_flow
          NAME Version: NAME III (version 6.3.1)
          Quantity: Total deposition
          Release height: 0.000 to 50.000m agl
          Release location: 2.9175W   54.0303N
          Release rate: 4.6305558E-05Bq/s
          Run time: 0947UTC 01/12/2015
          Species: CAESIUM-137
          Species Category: RADIONUCLIDE
          Start of release: 1800UTC 27/11/2015
          Time Av or Int: 048 hr time integrated
          Title: Heysham_27112015_1800Z
          label: NAME
          short_name: CAESIUM-137_TOTAL_DEPOSITION
     Cell methods:
          sum: time

Field Layers
^^^^^^^^^^^^

A :class:`field_layer.FieldLayer` object contains a cube and information
about how it should be plotted.

Initialise the field layer by giving it a cube:

>>> fl = adaqcode.field_layer.FieldLayer(cube)

The user can then set information about the layer styling, (for a full
list of options see :func:`field_layer.FieldLayer.set_layerstyle`) eg:

>>> fl.set_layerstyle(nlevels=10,
... plottype='pcolormesh',
... mask=True,
... autozoom=True,
... colorscale='log',
... cmap='YlGnBu')
>>> import numpy as np
>>> np.set_printoptions(formatter={'float':lambda x: '{:5.2e}'.format(x)})
>>> print(fl) # doctest: +ELLIPSIS
<class '...field_layer.FieldLayer'>
autozoom: True
cbar: None
cbar_label: default
cbar_num_fmt: None
cbar_orientation: None
cmap: YlGnBu (<matplotlib.colors.LinearSegmentedColormap object at ...>)
colors: None
colorscale: log
cube: CAESIUM-137_TOTAL_DEPOSITION / (Bq/m2) (latitude: 205; longitude: 245)
     Dimension coordinates:
          latitude                              x               -
          longitude                             -               x
     Scalar coordinates:
          source_latitude: 54.0303 degrees
          source_longitude: -2.9175 degrees
          time: 2015-11-29 18:00:00, bound=(2015-11-27 18:00:00, 2015-11-29 18:00:00)
          z: Boundary layer
     Attributes:
          End of release: 0000UTC 28/11/2015
          Forecast duration: 48 hours
          Met data: NWP Flow.UKV_PT1_flow; NWP Flow.UKV_PT2_flow
          NAME Version: NAME III (version 6.3.1)
          Quantity: Total deposition
          Release height: 0.000 to 50.000m agl
          Release location: 2.9175W   54.0303N
          Release rate: 4.6305558E-05Bq/s
          Run time: 0947UTC 01/12/2015
          Species: CAESIUM-137
          Species Category: RADIONUCLIDE
          Start of release: 1800UTC 27/11/2015
          Time Av or Int: 048 hr time integrated
          Title: Heysham_27112015_1800Z
          label: NAME
          short_name: CAESIUM-137_TOTAL_DEPOSITION
     Cell methods:
          sum: time
cube_crs: <cartopy._crs.Geodetic object at ...>
cube_extent: [-4.10..., 7.95...]
label: None
levels: [3.16e-14 1.00e-13 3.16e-13 1.00e-12 3.16e-12 1.00e-11 3.16e-11 1.00e-10
 3.16e-10 1.00e-09]
marker: o
markersize: 20
mask: True
nlevels: 10
norm: <matplotlib.colors.BoundaryNorm object at ...>
plottype: pcolormesh
step: 0.5
<BLANKLINE>

>>> np.set_printoptions()

The layer can be 'sliced' for use in plotting
(see :func:`field_layer.FieldLayer.layer_slice`).
Each slice returns a 2D segment of the cube, while retaining information
about how to plot the layer (eg levels and colours).

These slices can be used as an iterator in a loop for plotting different images,
particularly important if working with a 3D cube eg:

>>> for fl_slice in fl.layer_slice(['longitude','latitude']):
...     print(fl_slice.cube.summary(shorten=True))
CAESIUM-137_TOTAL_DEPOSITION / (Bq/m2) (longitude: 245; latitude: 205)


Field Plot objects
^^^^^^^^^^^^^^^^^^

A :class:`field_plot.FieldPlot` object contains some 2D slices of field layers
to plot and information about the plot.

>>> from adaqcode import inifile
>>> ini_dict = inifile.get_inidict(defaultfilename='adaqcode/inifile.ini')  # doctest: +ELLIPSIS
Reading inifile .../adaqcode/inifile.ini

>>> fp = adaqcode.field_plot.FieldPlot(ini_dict)

It requires at least one 2D cube slice which is added as a layer to the plot:

>>> fp.add_layer(fl)

Mapping can then also be added (for a full
list of options see :func:`field_plot.FieldPlot.setup_mapping`):

>>> fp.setup_mapping(
... extent=[-12, 4.0, 48., 60.],
... mapping='coastlines',
... gridlines=True)

The plot can then be produced, shown to screen (plt.show), saved if required,
or else the matplotlib figure object is returned for further use or adjustment:

>>> fig = fp.plot()
>>> #fp.save_fig(plotdir="outputdir/")
>>> import matplotlib.pyplot as plt
>>> #plt.show()
>>> plt.close() #Close this figure as we don't need it any longer.

.. image:: ../figures/Fieldplot_NAME_CAESIUM-137_TOTAL_DEPOSITION_Boundarylayer_201511291800.png

.. Note to developers. If you change any of the questions below then you will
   also need to change the answers in answers.rst

.. _plottinggridded_questions_ref:

Plotting Gridded Data Exercises
-------------------------------

  1. Plotting gridded fields automatically. Load in some NAME data using the inifile
     *name_field_plot.ini* in the adaqscripts directory and plot the field with the
     short name *VOLCANIC_ASH_AIR_CONCENTRATION* using the *NAME* defaults. You may
     also find it helpful to take a look at the inifile with a text editor or print
     out the contents of inifile so you can see what variables are set.

     a. Load in the NAME data using the inifile *name_field_plot.ini*. (Hint: You did
        this in question 4 in :ref:`loadingdata_questions_ref`)
     b. Change the location of the plot_dir to somewhere in your home space. Here the
        user_guide directory is only used for the purpose of the doctests.
     c. Plot the data setting the short name to *VOLCANIC_ASH_AIR_CONCENTRATION* and
        using the *NAME* defaults for plotting.

  2. If you look at the inifile *name_field_plot.ini* you will see that only the extent
     and the levels have been set. Try replotting the *VOLCANIC_ASH_AIR_CONCENTRATION*
     with a different colourmap (*cmap*) and changing the extent to focus on Iceland.

     a. Change the colourmap - some alternative colourmaps can be seen here:
        http://matplotlib.org/1.3.1/examples/color/colormaps_reference.html
     b. Change the extent_list. Note that the numbers for the extent list can be provided
        as strings or floats. They are converted into floats (if necessary) in
        *plot_md_gridded_fields*
     c. Re-run *plot_md_gridded_fields*. Don't forget to add a *filesuffix* so you
        don't overwrite your first plot.

Answers to :ref:`plottinggridded_answers_ref`

