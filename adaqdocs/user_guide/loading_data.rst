Loading data
============

In this section, we will explain how to load some data into ADAQData objects.

We begin by ensuring the adaqcode is imported. This is important both for
automated testing of this user guide page, but also if you want to follow
this page yourself.

>>> import sys
>>> adaq_path = '../../' #Full address of checked out code
>>> sys.path.append(adaq_path)
>>> import adaqcode

Loading AQUM pp data
--------------------

Here we give an example of loading AQUM pp data into the ADAQData objects.
This will go through the explict code required to achieve this. There are
shortcut routines to allow this to be done using inifiles which are described
further down.

Set up required variables
^^^^^^^^^^^^^^^^^^^^^^^^^

Firstly set up your data directory containing pp files.
We will limit the species we load using short_name, and the dates we
load using the datetime module.
Here we will make use of data from the SAMPLE_DATADIR directory:

>>> directory = adaqcode.SAMPLE_DATADIR+'aqum_output/casestudy/'
>>> print('Directory:'); print(directory) #  doctest: +ELLIPSIS
Directory:
.../aqum_output/casestudy/

Set up some variables to limit the data by:

>>> import datetime
>>> short_name_list = ['O3','NO2','PM2p5']
>>> start_datetime = datetime.datetime(2014,3,28,0)
>>> end_datetime = datetime.datetime(2014,3,29,0)

Get filenames
^^^^^^^^^^^^^

Generate a list of filenames manually, or by using a directory listing:

>>> import os
>>> filenames = [directory + filename for filename in os.listdir(directory)]
>>> print(sorted(filenames)) #  doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
['.../aqum_output/casestudy/aqum_20140324_T+012_QM18.pp',
'.../aqum_output/casestudy/aqum_20140324_T+024_QM18.pp',
'.../aqum_output/casestudy/aqum_20140325_T+000_QM18.pp',
'.../aqum_output/casestudy/aqum_20140325_T+012_QM18.pp',
'.../aqum_output/casestudy/aqum_20140325_T+024_QM18.pp',
...
'.../aqum_output/casestudy/aqum_20140403_T+024_QM18.pp',
'.../aqum_output/casestudy/aqum_20140404_T+000_QM18.pp',
'.../aqum_output/casestudy/aqum_20140404_T+012_QM18.pp']

Or for a more limited set of files (very AQUM specific), based on
start/end dates, forecast day etc. In this example we will request output from
the first full day of the forecast (day 1).

>>> ppfiles = adaqcode.pp_data.AQUMppFiles()
>>> filenames = ppfiles.get_filenames(directory,
...    start_datetime = start_datetime, end_datetime = end_datetime,
...    forecast_day = 1)
>>> print(filenames) #  doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
['.../aqum_output/casestudy/aqum_20140326_T+024_QM18.pp',
'.../aqum_output/casestudy/aqum_20140327_T+000_QM18.pp',
'.../aqum_output/casestudy/aqum_20140327_T+012_QM18.pp',
'.../aqum_output/casestudy/aqum_20140327_T+024_QM18.pp',
'.../aqum_output/casestudy/aqum_20140328_T+000_QM18.pp']


Read gridded pp data
^^^^^^^^^^^^^^^^^^^^

Initialise a :class:`pp_data.PPData` class - set this to the variable md
(model data). This class is a subclass of ADAQData so it has all the
functionality of the ADAQData objects, but also knows how to read in pp data.

>>> md = adaqcode.pp_data.PPData()
>>> print(md) #  doctest: +ELLIPSIS
<class '...pp_data.PPData'> - Subclass of ADAQData - Contains:
sites_cube_list:
< No cubes >
gridded_cube_list:
< No cubes >
trajectory_cube_list:
< No cubes >

Read data in, giving this model data a label, eg 'Casestudy'. We use the
variables we set up earlier containing short names, dates and filenames.
We also request data just from day 1, which here is classed as the first full
day of data covering 01Z - 24Z.

>>> md.readdata(filenames = filenames,
... short_name_list = short_name_list,
... start_datetime = start_datetime,
... end_datetime = end_datetime,
... forecast_day = 1,
... label='Casestudy')
[<iris 'Cube' of mass_concentration_of_pm2p5_dry_aerosol_in_air / (ug m-3) (time: 25; grid_latitude: 182; grid_longitude: 146)>,
<iris 'Cube' of mass_fraction_of_nitrogen_dioxide_in_air / (kg kg-1) (time: 25; grid_latitude: 182; grid_longitude: 146)>,
<iris 'Cube' of mass_fraction_of_ozone_in_air / (kg kg-1) (time: 25; grid_latitude: 182; grid_longitude: 146)>]

>>> print(md) #  doctest: +ELLIPSIS
<class '...pp_data.PPData'> - Subclass of ADAQData - Contains:
sites_cube_list:
< No cubes >
gridded_cube_list:
0: mass_concentration_of_pm2p5_dry_aerosol_in_air / (ug m-3) (time: 25; grid_latitude: 182; grid_longitude: 146)
1: mass_fraction_of_nitrogen_dioxide_in_air / (kg kg-1) (time: 25; grid_latitude: 182; grid_longitude: 146)
2: mass_fraction_of_ozone_in_air / (kg kg-1) (time: 25; grid_latitude: 182; grid_longitude: 146)
trajectory_cube_list:
< No cubes >

Extracting a sites cube list
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To extract a sites cube from this data, you need to use an ndarray. This is
provided from the :func:`sites_info.SitesInfo.read_from_file` method
of a :class:`sites_info.SitesInfo` object.
This ndarray gives information about the site locations:

>>> si = adaqcode.sites_info.SitesInfo()
>>> sitesfile = adaqcode.SAMPLE_DATADIR+'AURN_obs/aq_sites_GEMSRAQ_v4b_dev.txt'
>>> sites_data = si.read_from_file(sitesfile, allsites=True)
Number of sites:  5
>>> fmtstr = adaqcode.sites_info.format_string(sites_data)
>>> for site in sites_data:
...     print(fmtstr.format(*site))
GB0001, ABD,  57.158, -2.094,  20,URBAN_BACKGROUND,        Aberdeen
GB0003,ACTH,  55.883, -3.347, 260,           RURAL,Auchencorth_Moss
GB0002,  AH,  52.504, -3.034, 370,           RURAL,      Aston_Hill
GB0045, HAR,  51.571, -1.327, 137,           RURAL,         Harwell
GB0128,  YW,  50.597, -3.716, 119,           RURAL,     Yarner_Wood

This can then be used to extract the sites cubes:

>>> md_scl = md.extract_sites(sites_data)
>>> print(md) #  doctest: +ELLIPSIS
<class '...pp_data.PPData'> - Subclass of ADAQData - Contains:
sites_cube_list:
0: mass_concentration_of_pm2p5_dry_aerosol_in_air / (ug m-3) (site_id: 5; time: 25)
1: mass_fraction_of_nitrogen_dioxide_in_air / (kg kg-1) (site_id: 5; time: 25)
2: mass_fraction_of_ozone_in_air / (kg kg-1) (site_id: 5; time: 25)
gridded_cube_list:
0: mass_concentration_of_pm2p5_dry_aerosol_in_air / (ug m-3) (time: 25; grid_latitude: 182; grid_longitude: 146)
1: mass_fraction_of_nitrogen_dioxide_in_air / (kg kg-1) (time: 25; grid_latitude: 182; grid_longitude: 146)
2: mass_fraction_of_ozone_in_air / (kg kg-1) (time: 25; grid_latitude: 182; grid_longitude: 146)
trajectory_cube_list:
< No cubes >

To get a cube at a single site from this, perhaps for a single species,
eg Harwell NO2:

>>> har_no2 = md.extract(site_name='Harwell', short_name='NO2', singlecube=True)
>>> print(har_no2) # doctest: +ELLIPSIS
mass_fraction_of_nitrogen_dioxide_in_air / (kg kg-1) (time: 25)
     Dimension coordinates:
          time                                            x
     Auxiliary coordinates:
          forecast_period                                 x
     Scalar coordinates:
          abbrev: HAR
          forecast_day: 1.0 Days
          grid_latitude: -0.922... degrees
          grid_longitude: 0.729... degrees
          latitude: 51.571... degrees
          level_height: 20.000... m, bound=(0.0, 49.998...) m
          longitude: -1.326... degrees
          model_level_number: 1
          sigma: 0.997..., bound=(1.0, 0.994...)
          site_altitude: 137 m
          site_id: 35867333.141...
          site_name: Harwell
          site_type: RURAL
     Attributes:
          STASH: m01s34i004
          label: Casestudy
          short_name: NO2
          source: Data from Met Office Unified Model
     Cell methods:
          mean: time (1 hour)

Note this extraction of sites has done a simple phenomenon conversion to mass
concentration (which therefore changes the units from kg/kg to ug/m3).
This only occurs for PPData to be consistent with air quality
observations.

Or to get data from the gridded cube instead:

>>> no2 = md.extract(short_name='NO2', singlecube=True, gridded=True)


Phenomenon / Unit Conversion
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To do a simple unit conversion to ug/m3, making the assumption of
standard temperature and pressure (this is usually assumed for
surface fields) use
:func:`cube_chemistry.cube_to_mass_conc_in_air`.

>>> no2_ugm3 = adaqcode.cube_chemistry.cube_to_mass_conc_in_air(no2,at_stp=True)
>>> print('{:.3e}'.format(no2.data.max()))
8.395e-08
>>> print('{:.3f}'.format(no2_ugm3.data.max()))
101.134

Note if you need to convert eg SO2 which is expressed as S, then you also
need to use :func:`cube_chemistry.cube_elemental_to_molecular_mass` first.

For doing all gridded and site data, for observations and model data, then
:func:`adaq_functions.unit_conversion` can be useful.

For converting 3D data, e.g. to examine vertical profiles, it is invalid to
use stp and therefore 3D pressure and temperature cubes must also be supplied.
An example is provided in the documentation for
:func:`adaq_functions.unit_conversion`.
This function can also be used for converting surface fields, either with or
without making the assumption of stp.

Loading NAME data
-----------------

Set up required variables
^^^^^^^^^^^^^^^^^^^^^^^^^

Firstly set up your data directory containing name files.
Here we will use data from the SAMPLE_DATADIR directory:

>>> directory = adaqcode.SAMPLE_DATADIR+'name/'
>>> print('Directory:'); print(directory) #  doctest: +ELLIPSIS
Directory:
.../name/

Set up a field attributes dictionary - this will allow you to limit the
data you extract by setting particular variables from the NAME file.
This is a dictionary of the attributes in a cube once loaded in.
For more information on this, see :ref:`field_attributes_ref`

>>> field_attributes = {'Species':'CAESIUM-137'}

Initalise the NAMEData class (:class:`name_data.NAMEData`):

>>> name = adaqcode.name_data.NAMEData()


Get filenames
^^^^^^^^^^^^^

Filenames can be defined simply by providing a list of filenames:

>>> filenames = [directory + 'Fields_grid1*.txt']
>>> print(filenames) #  doctest: +ELLIPSIS
['.../name/Fields_grid1*.txt']

Or by limiting according to required start and end dates and using the
:func:`name_data.NAMEData.get_filenames` function, which stores the filenames
within the NAMEData object for later use.

>>> filenames = name.get_filenames(directory,  file_pattern = 'Fields_grid1',
... start_datetime=datetime.datetime(2011, 0o3, 13),
... end_datetime=datetime.datetime(2011, 0o3, 18))
>>> print(filenames) #  doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
['.../name/Fields_grid1_C1_T2_201103130000.txt',
'.../name/Fields_grid1_C1_T3_201103140000.txt',
'.../name/Fields_grid1_C1_T4_201103150000.txt',
'.../name/Fields_grid1_C1_T5_201103160000.txt',
'.../name/Fields_grid1_C1_T6_201103170000.txt',
'.../name/Fields_grid1_C1_T7_201103180000.txt']


Reading gridded NAME data
^^^^^^^^^^^^^^^^^^^^^^^^^

Gridded model data can be read in using :func:`name_data.NAMEData.readdata`:

>>> gcl = name.readdata(filenames, field_attributes=field_attributes)
>>> print(gcl)
0: CAESIUM-137_AIR_CONCENTRATION / (Bq / m^3) (time: 6; latitude: 90; longitude: 180)
>>> print(gcl[0])
CAESIUM-137_AIR_CONCENTRATION / (Bq / m^3) (time: 6; latitude: 90; longitude: 180)
     Dimension coordinates:
          time                                  x            -              -
          latitude                              -            x              -
          longitude                             -            -              x
     Scalar coordinates:
          height: 1000.0 m
     Attributes:
          End of release: 20/04/2011 00:00 UTC
          Ensemble Av: No ensemble averaging
          Horizontal Av or Int: No horizontal averaging
          Met data: NWP Flow.Global_PT1_flow; NWP Flow.Global_PT2_flow; NWP Flow.Global_PT3_flow;...
          NAME Version: NAME III (version 6.3)
          Name: Unnamed Field Req 5
          Quantity: Air Concentration
          Release height: Multiple Sources
          Release location: Multiple Sources
          Run duration: 91day 9hr 0min
          Run name: Stohl2degFukushimaTake2
          Run time: 10/02/2014 11:37:15.277 UTC
          Source strength: Multiple Sources
          Sources: All sources
          Species: CAESIUM-137
          Species Category: RADIONUCLIDE
          Start of release: 10/03/2011 12:00 UTC
          Time Av or Int: 1day 0hr 0min average
          Vertical Av or Int: No vertical averaging
          label: NAME
          short_name: CAESIUM-137_AIR_CONCENTRATION
     Cell methods:
          mean: time

Reading NAME time-series data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Name time-series data can be read in a similar way, but will be read into
a sites cube list instead:

>>> filenames = [directory + 'Time_series_grid1.txt']
>>> name_ts = adaqcode.name_data.NAMEData()
>>> scl = name_ts.readdata(filenames)
>>> print(name_ts) #  doctest: +ELLIPSIS
<class '...name_data.NAMEData'> - Subclass of ADAQData - Contains:
sites_cube_list:
0: VOLCANIC_ASH_AIR_CONCENTRATION / (g/m3) (site_id: 2; time: 132)
gridded_cube_list:
< No cubes >
trajectory_cube_list:
< No cubes >
>>> print(name_ts.sites_cube_list[0])
VOLCANIC_ASH_AIR_CONCENTRATION / (g/m3) (site_id: 2; time: 132)
     Dimension coordinates:
          site_id                               x        -
          time                                  -        x
     Auxiliary coordinates:
          latitude                              x        -
          longitude                             x        -
          site_abbrev                           x        -
          site_name                             x        -
     Scalar coordinates:
          site_type: Unknown
          source_latitude: 63.63 degrees
          source_longitude: -19.05 degrees
          z: Boundary layer
     Attributes:
          End of release: 0830UTC 28/03/2022
          Forecast duration: 87648 hours
          Met data: NWP Flow.Global_PT1_flow; NWP Flow.Global_PT2_flow
          NAME Version: NAME III (version 6.0 for VAAC)
          Quantity: Air Concentration
          Release height: 8376.000m asl +/- 6864.000m
          Release location: 19.0500W   63.6300N
          Release rate: 2.5890222E+08g/s
          Run time: 1358UTC 30/04/2012
          Species: VOLCANIC_ASH
          Species Category: VOLCANIC
          Start of release: 0830UTC 28/03/2012
          Title: KATLA_20120401_12Z
          label: NAME
          short_name: VOLCANIC_ASH_AIR_CONCENTRATION

Loading NAME data into an Ensemble
----------------------------------
It is also possible to load NAME data into an ensemble. I.e. if we load
NAME data from runs which use different met data we can merge them into
a single cube with a realization coordinate. It is then possible to use
the functionality in Iris to generate probabilities of exceedance and
percentiles of concentration from this single cube.

The data is loaded using a similar method to the loading of gridded
NAME data but now we need to specify an ensemble dictionary to link
each the directory containing each ensemble member with an ensemble
label (called and ensemble_member coordinate) and a numeric identifier
which is used to create a dimension coordinate in the iris cube. It is
necessary to specify the numeric identifier because python dictionaries
do not guarantee order so if the numeric identifier was omitted the
ensemble members may not be numbered in the order they were listed.

First we define the filenames:

>>> name_ens = adaqcode.name_data.NAMEData()
>>> directory = adaqcode.SAMPLE_DATADIR+'name_ensemble/member[012]/'
>>> filenames = [directory + 'Fields_grid4*.txt']

Then we define the ensemble dictionary. In this case we have a ten
member ensemble and each ensemble member is stored in a subdirectory
with the name member* where * is a number from zero to nine.

>>> ensemble_dict={adaqcode.SAMPLE_DATADIR+'name_ensemble/member0/': ['member0', 0],
... adaqcode.SAMPLE_DATADIR+'name_ensemble/member1/': ['member1', 1],
... adaqcode.SAMPLE_DATADIR+'name_ensemble/member2/': ['member2', 2]}

Now load the data:

>>> gcl = name_ens.readdata(filenames, ensemble_dict=ensemble_dict,
... field_attributes = {'Field Name':'Unnamed Field Req 8'})

We can look at the two new coordinates; the dimension coordinate
'realization' and the auxilary coordinate 'ensemble_member':

>>> print(gcl[0].coord('realization'))
DimCoord(array([0, 1, 2]), standard_name='realization', units=Unit('1'))

>>> print(gcl[0].coord('ensemble_member').points)
['member0' 'member1' 'member2']

Loading NAME trajectory data
----------------------------
Firstly set up your data directory containing name trajectory files.
Here we will use data from the SAMPLE_DATADIR directory:

>>> directory = adaqcode.SAMPLE_DATADIR+'name_trajectory/'
>>> print('Directory:'); print(directory) #  doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
Directory:
.../name_trajectory/

Initialise the TrajData class (trajectory_data.TrajData):

>>> traj = adaqcode.trajectory_data.TrajData()

Filenames can be defined by simply providing a list of filenames:

>>> filenames = [directory + 'Data_Traj_C1_*.txt']
>>> print(filenames)#  doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
['.../name_trajectory/Data_Traj_C1_*.txt']

Trajectory data can be read in using traj_data.TrajData.readdata():

>>> tcl = traj.readdata(filenames)
>>> print(tcl)
0: U Turb / (m/s)                      (time: 193)
1: U Turb / (m/s)                      (time: 193)
>>> print(tcl[0])
U Turb / (m/s)                      (time: 193)
     Dimension coordinates:
          time                           x
     Auxiliary coordinates:
          Travel Time                    x
          altitude                       x
          latitude                       x
          longitude                      x
     Scalar coordinates:
          PP Index: 2.0
          Release Time: 2016-12-02 08:00:00
          Source: Source_2
     Attributes:
          Met data: NWP Flow.Global_PT1_flow; NWP Flow.Global_PT2_flow; NWP Flow.Global_PT3_flow;...
          NAME Version: NAME III (version 6.5.1)
          Run name: London 02122016 0800
          Run time: 08/12/2016 08:27:38.939 UTC
          Trajectory direction: Forward trajectories
          label: TRAJ
          short_name: U Turb


NAME trajectory files contain a large number of columns of data
such as met data, plume rise data, and turbulence data. The trajectory
class will try to read in the 'U Turb variable because it always exist.
However, it is possible to specify alternate variables to read in by simply
specifying a different short_name when calling readdata. Lets try reading
in the temperature variable:

>>> tcl = traj.readdata(filenames, short_name='Temperature')
>>> print(tcl)
0: Temperature / (K)                   (time: 193)
1: Temperature / (K)                   (time: 193)

Loading generic gridded data
----------------------------

.. todo: Add something here about setting directly from iris cubes.

In theory any data that has been read into iris cubes can also be read into
an ADAQData object. However there are are few additional rules that you should
be aware of. Let us try loading some UK NWP data:

>>> import iris
>>> sample_data_path = adaqcode.SAMPLE_DATADIR + 'ukhires/'
>>> cubes = iris.load(sample_data_path + 'uk_hires.pp', 'air_potential_temperature')
>>> print(cubes)
0: air_potential_temperature / (K)     (time: 3; model_level_number: 7; grid_latitude: 204; grid_longitude: 187)

Now set up your ADAQData object:

>>> iris_data = adaqcode.adaq_data.ADAQData()

The cubes loaded in are gridded data (they have an X and Y coordinate), so
put them into the gridded_cube_list:

>>> iris_data.gridded_cube_list = cubes

You should now always check that this is a valid gridded_cube_list. To do
this, use the function :func:`adaq_data.ADAQData.check_gridded_cube`:

>>> iris_data.check_gridded_cube()
Traceback (most recent call last):
  File "<pyshell#10>", line 1, in <module>
    iris_data.check_gridded_cube()
  File "adaqcode/adaq_data.py", line 375, in check_gridded_cube
    self.__check_attribute(gridded_cube, attrib_name, verbose)
  File "adaqcode/adaq_data.py", line 300, in __check_attribute
    raise ValueError, attrib_name + ' is not in attributes'
ValueError: short_name is not in attributes

You can see that this fails with the error "short_name is not in attributes".
This is because one of the key features of ADAQData is that within the
attributes of a cube is the "short_name". For air_potential_temperature this
could be "theta". To check exactly what attributes are required, have a look at
GRIDDED_CUBE_REQUIRED_ATTRIBUTES. You can also check which dimension coordinates
are required to be a gridded cube:

>>> print(iris_data.GRIDDED_CUBE_REQUIRED_ATTRIBUTES)
['short_name', 'label']
>>> print(iris_data.GRIDDED_CUBE_REQUIRED_COORD_AXES)
['T', 'X', 'Y']

We therefore need to add both the short_name and a label to our gridded cube:

>>> for icube in range(len(iris_data.gridded_cube_list)):
...     iris_data.gridded_cube_list[icube].attributes['short_name'] = 'theta'
...     iris_data.gridded_cube_list[icube].attributes['label'] = 'IrisData'
>>> print(iris_data.gridded_cube_list[0])
air_potential_temperature / (K)     (time: 3; model_level_number: 7; grid_latitude: 204; grid_longitude: 187)
     Dimension coordinates:
          time                           x                      -                 -                    -
          model_level_number             -                      x                 -                    -
          grid_latitude                  -                      -                 x                    -
          grid_longitude                 -                      -                 -                    x
     Auxiliary coordinates:
          forecast_period                x                      -                 -                    -
          level_height                   -                      x                 -                    -
          sigma                          -                      x                 -                    -
          surface_altitude               -                      -                 x                    x
     Derived coordinates:
          altitude                       -                      x                 x                    x
     Scalar coordinates:
          forecast_reference_time: 2009-11-19 04:00:00
     Attributes:
          STASH: m01s00i004
          label: IrisData
          short_name: theta
          source: Data from Met Office Unified Model
          um_version: 7.3

Now check again:

>>> iris_data.check_gridded_cube()
short_name - OK
label - OK
T axis - OK - coord= time
X axis - OK - coord= grid_longitude
Y axis - OK - coord= grid_latitude

Sites cubes could also be done in a similar way, but these need to have a
site_id coordinate so are more tricky and therefore best avoided.

.. _loading_using_inifiles:

Loading data using Inifiles
---------------------------

Data can also be loaded in through the use of inifiles.
This can be particularly useful to allow scripts to be reused with multiple
different input options.

As an example, read in an example file in .ini format:

>>> inifilename = adaqcode.CODE_DIR + 'adaqcode/example_data_1days.ini'

The contents of this are:

.. literalinclude:: ../../adaqcode/example_data_1days.ini

Now read this in using :func:`inifile.get_inidict`:

>>> ini_dict = adaqcode.inifile.get_inidict(inifilename) # doctest: +ELLIPSIS
Reading inifile .../adaqcode/example_data_1days.ini

Also as we would like to extract some sites cubes, get some site information.
We have given a sites file as part of the ini_dict which we can use:

>>> print(ini_dict['sites_file']) # doctest: +ELLIPSIS
/.../aq_sites_GEMSRAQ_v4b_dev.txt

We can retrieve this sites information, simply by passing in the ini_dict
into :func:`sites_info.get_siteinfo`

>>> sites_data = adaqcode.sites_info.get_siteinfo(ini_dict)
Number of sites:  5

Now load in some observations into sites cubes. For this we use the format and
directory of observations which are defined in the ini_dict:

>>> print(ini_dict['obs_fmt'])
aurn
>>> print(ini_dict['obs_dir']) # doctest: +ELLIPSIS
/.../AURN_obs

Again though we only need to pass the ini_dict into a routine, along with
the sites_data that we loaded. For this we use the routine
:func:`adaq_functions.get_obs`:

>>> od = adaqcode.adaq_functions.get_obs(ini_dict, sites_data) # doctest: +ELLIPSIS
Creating observation data at  ...
Reading obs data files
Found obs for  .../AURN_obs/ABD_20140101_20140818.txt
Found obs for  .../AURN_obs/ACTH_20140101_20140818.txt
Found obs for  .../AURN_obs/AH_20140101_20140818.txt
Found obs for  .../AURN_obs/HAR_20140101_20140818.txt
Found obs for  .../AURN_obs/YW_20140101_20140818.txt
Creating obs cubes

This returns an ADAQData object, with the sites cube list filled with
observation data:

>>> print(od) #  doctest: +ELLIPSIS
<class '...aurn_data.AURNData'> - Subclass of ADAQData - Contains:
sites_cube_list:
0: mass_concentration_of_ozone_in_air / (ug/m3) (site_id: 5; time: 25)
1: mass_concentration_of_pm2p5_ambient_aerosol_in_air / (ug/m3) (site_id: 5; time: 25)
gridded_cube_list:
< No cubes >
trajectory_cube_list:
< No cubes >

Now get some model data. For this we use other variables given in the inifile.
Often we are trying to compare one model to another, or multiple runs of
the same model, so we will work with a list of model data objects.
The inifile therefore contains a list for each
keyword, with each value in the list corresponding to a model run. For example:

>>> print(ini_dict['models_fmt_list'])
['pp', 'pp']
>>> print(ini_dict['models_dir_list']) # doctest: +ELLIPSIS
['.../aqum_output/oper', '.../casestudy']
>>> print(ini_dict['models_list'])
['aqum_oper', 'aqum_casestudy']

So with this output, the first model is in pp format, and will be given the
label 'aqum_oper'. The second model is also in 'pp' format, and will be given
the label 'aqum_casestudy'. Note the models_list was not given in the inifile
(it is commented out), but is set using the final subdirectory name.
The model data can then be read in easily using
:func:`adaq_functions.get_models`:

>>> md_list = adaqcode.adaq_functions.get_models(ini_dict, sites_data) # doctest: +ELLIPSIS
Getting model data for  aqum_oper  at ...
Getting model data for  aqum_casestudy  at  ...

This returns a list of ADAQData objects, each with gridded model data, but also
(as sites_data was not None), the sites_cube_list also filled as this has been
extracted from the gridded_cube_list:

>>> for md in md_list:
...     print(md) #  doctest: +ELLIPSIS
<class '...pp_data.PPData'> - Subclass of ADAQData - Contains:
sites_cube_list:
0: mass_concentration_of_pm2p5_dry_aerosol_in_air / (ug m-3) (site_id: 5; time: 25)
1: mass_fraction_of_ozone_in_air / (kg kg-1) (site_id: 5; time: 25)
gridded_cube_list:
0: mass_concentration_of_pm2p5_dry_aerosol_in_air / (ug m-3) (time: 25; grid_latitude: 182; grid_longitude: 146)
1: mass_fraction_of_ozone_in_air / (kg kg-1) (time: 25; grid_latitude: 182; grid_longitude: 146)
trajectory_cube_list:
< No cubes >
<class '...pp_data.PPData'> - Subclass of ADAQData - Contains:
sites_cube_list:
0: mass_concentration_of_pm2p5_dry_aerosol_in_air / (ug m-3) (site_id: 5; time: 25)
1: mass_fraction_of_ozone_in_air / (kg kg-1) (site_id: 5; time: 25)
gridded_cube_list:
0: mass_concentration_of_pm2p5_dry_aerosol_in_air / (ug m-3) (time: 25; grid_latitude: 182; grid_longitude: 146)
1: mass_fraction_of_ozone_in_air / (kg kg-1) (time: 25; grid_latitude: 182; grid_longitude: 146)
trajectory_cube_list:
< No cubes >

Loading some example data
-------------------------

To enable quick access to some data for testing, a routine
:func:`adaq_functions.get_exampledata` is available. This returns an ini_dict,
sites_data, observation data and a list of model data objects:

>>> ini_dict, sites_data, od, md_list = adaqcode.adaq_functions.get_exampledata() # doctest: +ELLIPSIS
Reading inifile .../adaqcode/example_data_1days.ini
Number of sites:  5

>>> print(od) #  doctest: +ELLIPSIS
<class '...adaq_data.ADAQData'> - Subclass of ADAQData - Contains:
sites_cube_list:
0: mass_concentration_of_ozone_in_air / (ug/m3) (site_id: 5; time: 25)
1: mass_concentration_of_pm2p5_ambient_aerosol_in_air / (ug/m3) (site_id: 5; time: 25)
gridded_cube_list:
< No cubes >
trajectory_cube_list:
< No cubes >

>>> for md in md_list:
...    print(md) #  doctest: +ELLIPSIS
<class '...adaq_data.ADAQData'> - Subclass of ADAQData - Contains:
sites_cube_list:
0: mass_concentration_of_ozone_in_air / (ug/m3) (site_id: 5; time: 25)
1: mass_concentration_of_pm2p5_dry_aerosol_in_air / (ug/m3) (site_id: 5; time: 25)
gridded_cube_list:
< No cubes >
trajectory_cube_list:
< No cubes >
<class '...adaq_data.ADAQData'> - Subclass of ADAQData - Contains:
sites_cube_list:
0: mass_concentration_of_ozone_in_air / (ug/m3) (site_id: 5; time: 25)
1: mass_concentration_of_pm2p5_dry_aerosol_in_air / (ug/m3) (site_id: 5; time: 25)
gridded_cube_list:
< No cubes >
trajectory_cube_list:
< No cubes >

Note by default this only sets up the sites_cube_list in the model data,
to get the gridded data as well, set the keyword gridded_cube_list=True.
To get a slightly longer period of example data (5 days) with more species,
set the keyword exampletype='full', eg:

>>> ini_dict, sites_data, od, md_list = adaqcode.adaq_functions.get_exampledata(gridded_cube_list=True,
... sites_cube_list=False, exampletype='full') # doctest: +ELLIPSIS
Reading inifile .../adaqcode/example_data_5days.ini
>>> for md in md_list:
...    print(md) #  doctest: +ELLIPSIS
<class '...adaq_data.ADAQData'> - Subclass of ADAQData - Contains:
sites_cube_list:
< No cubes >
gridded_cube_list:
0: mass_concentration_of_nitrogen_dioxide_in_air / (ug/m3) (time: 121; grid_latitude: 182; grid_longitude: 146)
1: mass_concentration_of_ozone_in_air / (ug/m3) (time: 121; grid_latitude: 182; grid_longitude: 146)
2: mass_concentration_of_pm10_dry_aerosol_in_air / (ug/m3) (time: 121; grid_latitude: 182; grid_longitude: 146)
3: mass_concentration_of_pm2p5_dry_aerosol_in_air / (ug/m3) (time: 121; grid_latitude: 182; grid_longitude: 146)
4: mass_concentration_of_sulfur_dioxide_in_air / (ug/m3) (time: 121; grid_latitude: 182; grid_longitude: 146)
trajectory_cube_list:
< No cubes >
<class '...adaq_data.ADAQData'> - Subclass of ADAQData - Contains:
sites_cube_list:
< No cubes >
gridded_cube_list:
0: mass_concentration_of_nitrogen_dioxide_in_air / (ug/m3) (time: 121; grid_latitude: 182; grid_longitude: 146)
1: mass_concentration_of_ozone_in_air / (ug/m3) (time: 121; grid_latitude: 182; grid_longitude: 146)
2: mass_concentration_of_pm10_dry_aerosol_in_air / (ug/m3) (time: 121; grid_latitude: 182; grid_longitude: 146)
3: mass_concentration_of_pm2p5_dry_aerosol_in_air / (ug/m3) (time: 121; grid_latitude: 182; grid_longitude: 146)
4: mass_concentration_of_sulfur_dioxide_in_air / (ug/m3) (time: 121; grid_latitude: 182; grid_longitude: 146)
trajectory_cube_list:
< No cubes >

.. Note to developers. If you change any of the questions below then you will
   also need to change the answers in answers.rst

.. _loadingdata_questions_ref:

Loading Data Exercises
----------------------

  1. Load in data with the short name *PM10*, a validity date between
     0Z 01/04/2014 and 0Z 02/04/2014 from the first full day of forecast. This
     data can be found in the *aqum_output/oper* directory in the sample data
     directory and should be given the label 'Oper' (short for operational).

     a. Create a short_name, start_datetime and end_datetime variable to
        limit the data loading with.
     b. Setup a AQUMppfiles object and use it obtain the filenames for day
        1 forecasts.
     c. Set up a PPData object and pass all the variables (e.g. short_name,
        filenames) to it to load the data for day 1 forecasts.

  2. Load in the XENON-133 data from the Fields_grid1*.txt files stored
     in the *name* directory in the sample data directory. Then using the sites
     data stored in the *ctbto_sites.txt* file in the *ctbto_obs* directory in
     the sample data directory extract the Xenon-133 data at Sacramento and
     Takasaki.

     a. Setup the directory and some field_attributes to restrict loading
     b. Setup a NAMEData object and use it to read some files in
     c. Setup a SitesInfo object and pass the sites file location to it to read
     d. Use the sites_data to extract the sites_cubes from the NAMEData object

  3. Use the information in the inifile example_data_5days.ini to load some
     observations data.

     a. Read the information in adaqcode/example_data_5days.ini into an ini_dict
     b. Use the ini_dict to retrieve the sites information
     c. Load in the observations
     d. Print the returned ADAQ object

  4. Use the information in the inifile name_field_plot.ini to load some NAME
     model data.

     a. Read the information in adaqscripts/name_field_plot.ini into an ini_dict
     b. Now use the ini_dict to load in some NAME model data
     c. Print out the contents of the first item in the model list

Answers to :ref:`loadingdata_answers_ref`
