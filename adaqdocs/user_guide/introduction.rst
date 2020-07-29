Introduction
============

Welcome to ADAQPython!

In this user guide we will guide you through various aspects of ADAQPython
so by the end you should feel confident enough to have a go using this code
in ADAQPython for working with your own data to easily produce high
quality plots and to get statistics for model verification etc.

We will introduce the various component modules of the main code base to
enable you to use the main code directly, with just a little experience
of using python.
However it is also worth noting that if you just want to do some simple
plotting, some of the scripts (see :ref:`adaqscripts`) may be helpful
(for these you don't need any python knowledge!)

Please note that the ADAQPython code should work with the latest version
of iris and older versions are not guaranteed to work. Please therefore
make sure you are using the most recent version. In the Met Office this
can be ensured by doing:

.. code-block:: ksh

     module load scitools

In order to use the plotting code you will need to have access to map
data. This is setup through the use of cartopy config (see
:ref:`cartopy_setup_ref` for more information)


User guide syntax
-----------------

As you go through this user guide there is some syntax that is worth
being aware of:

>>> print("hello")
hello

Anything after ">>>" is python code that can be run in your favourite
python editor (eg idle, emacs) or if python has been lauched on the
command line.
On the following line is the output you would expect to find from running
this command (in this case hello)

You may also see places where "..." is used, eg:

>>> my_list = ["a", "b",
... "c", "d"]

The three dots here are just to say that the line has been continued -
it could equally have been written:

>>> my_list = ["a", "b", "c", "d"]

Alternatively, the three dots in the case of:

>>> for a in range(3):
...    print(a)
0
1
2

are used to illustrate the loop has been continued. Here we ensure that we
are always confirming to the coding standards of 4 spaces for indentation
(see :ref:`coding_standards` for futher information).

Another place where you may see three dots is in the output lines:

>>> print(list(range(100))) # doctest: +ELLIPSIS
[0, 1, 2, ..., 97, 98, 99]

Here it represents that some of the output is missing. This is often
used for user-specific locations, or particuarly long, irrelevant output.

Setting up your code to find ADAQ Python Code
---------------------------------------------

If this is the first time you have worked with ADAQ Python, then you will
need to initially set up a link to the code repository. This can be done
by opening $HOME/.fcm in your favourite text editor (creating this file if it
doesn't already exist) and then copy the following line into it:

   *set::url::adaq_python svn://fcm8/ADAQ_PythonCode_svn/ADAQ_PythonCode*

ADAQPython is not by default on your systems python path
(unlike numpy for example where you can directly import numpy).
You therefore need to tell any code firstly where the ADAQ Python code is.
Within the Met Office, there is a copy of the trunk
kept in ~apdg/PythonCode/ADAQ_Python/trunk. This should be updated whenever
any code changes are made in the trunk. However this could mean that the
code is changed while you are working with it. It is therefore recommended
that you firstly checkout a copy of the trunk for your own use (just be
careful not to make any changes which you then commit directly back to the
trunk!).
To check out a version of the trunk:

.. code-block:: ksh

   cd $DATADIR/adaq_python #Or another suitable directory
   fcm checkout fcm:adaq_python_tr

At this stage if you only need to run the ADAQ python scripts (rather than
write any python code), then this should be sufficient - refer to
:ref:`adaqscripts` for further instructions.

For using the ADAQ Python code in your own code, the location of where you
have checked out the code (for example in the above example this would be
/data/users/username/adaq_python/trunk) should then be added to your python
path at the top of any code, which then allows adaqcode to be imported:

>>> import sys
>>> adaq_path = '../../' #Full address of checked out code
>>> #Eg adaq_path = 'my/fcm/location/trunk'
>>> sys.path.append(adaq_path)
>>> import adaqcode

It is worth noting that there are several important variables that are
set up in :mod:`config` (see adaqcode/config.py) which gives locations where
certain directories can be found locally at your site (eg at the Met Office).
For the following couple of print statements, try printing these to screen for
yourself to see where they are installed for you.

Some of the variables set up are:

>>> print(adaqcode.CODE_DIR) # doctest: +ELLIPSIS
/.../adaqcode/../

This CODE_DIR points at the top level directory of the python - this directory
then contains the sub-directories adaqcode, adaqdocs, adaqscripts and adaqsandpit.
(See :ref:`directory_structure` for more information about what is contained in here.)
Note this location is specific to the version of the branch, or of the trunk, that
you have checked out from fcm.

>>> print(adaqcode.SAMPLE_DATADIR) # doctest: +ELLIPSIS
/.../python_sample_data/

The SAMPLE_DATADIR contains various sub-directories which have some sample data
files contained within them, for example for NAME or AQUM pp files. Try printing
out this variable and then having a look at the contents of this directory.


ADAQData objects
----------------

ADAQData objects are at the very centre of the adaqcode library.
These are built on top of the Iris cube and are extended to give additional
functionality to be useful for ADAQ users.


.. topic:: Aside - more information about Classes and Subclasses

   The ADAQData objects are based on a class, ADAQData. The process of
   extending a class is known in Python as creating subclasses. These
   subclasses have all the same features as the parent subclass (ADAQData),
   but also have additional functionality specific to themselves.
   Within ADAQPython, these subclasses extend the functionality to allow
   different data formats to be read in. For example one such subclass is
   NAMEData (:class:`name_data.NAMEData`), which allows NAME format files
   to be read in. Another subclass is PPData (:class:`pp_data.PPData`) to
   allow PP format files to be read in.

The ADAQData objects are extended (using subclasses) to allow the reading of
different format data, such as pp or NAME. Once these data have been read in
they can all be treated the same.

Try setting up a basic ADAQData object (see :class:`adaq_data.ADAQData`) to
hold some model data:

>>> md = adaqcode.adaq_data.ADAQData()
>>> print(md) # doctest: +ELLIPSIS
<class '...adaq_data.ADAQData'> - Subclass of ADAQData - Contains:
sites_cube_list:
< No cubes >
gridded_cube_list:
< No cubes >
trajectory_cube_list:
< No cubes >

This contains three main attributes: a sites_cube_list, a gridded_cube_list
and a trajectory cube list. Initially these are all set to empty iris cubelists
(iris.cube.CubeList()), which have a length of zero.

Gridded Cube List
^^^^^^^^^^^^^^^^^

A gridded cube list is an iris cube list containing cubes of gridded data.
Now load some example data into the gridded_cube_list to help explain this
component better:

>>> md_gcl = md.load_gridded(adaqcode.SAMPLE_DATADIR+'gridded_cube_list/aqum_oper_1days.nc')

If we again print the ADAQData object, we can see that the gridded_cube_list
has now been filled:

>>> print(md) # doctest: +ELLIPSIS
<class '...adaq_data.ADAQData'> - Subclass of ADAQData - Contains:
sites_cube_list:
< No cubes >
gridded_cube_list:
0: mass_concentration_of_ozone_in_air / (ug/m3) (time: 25; grid_latitude: 182; grid_longitude: 146)
1: mass_concentration_of_pm2p5_dry_aerosol_in_air / (ug/m3) (time: 25; grid_latitude: 182; grid_longitude: 146)
trajectory_cube_list:
< No cubes >

The gridded cube list can be accessed:

>>> gcl = md.gridded_cube_list
>>> print(gcl)
0: mass_concentration_of_ozone_in_air / (ug/m3) (time: 25; grid_latitude: 182; grid_longitude: 146)
1: mass_concentration_of_pm2p5_dry_aerosol_in_air / (ug/m3) (time: 25; grid_latitude: 182; grid_longitude: 146)

This is an iris cubelist. Each cube in the list contains gridded data.
It must contain data on an 'X' (eg longitude or grid_longitude),
'Y' (eg latitude or grid_latitude) and 'T' (eg time) axis.
Try printing the first cube in this list:

>>> print(gcl[0]) # doctest: +ELLIPSIS
mass_concentration_of_ozone_in_air / (ug/m3) (time: 25; grid_latitude: 182; grid_longitude: 146)
     Dimension coordinates:
          time                                    x                  -                    -
          grid_latitude                           -                  x                    -
          grid_longitude                          -                  -                    x
     Auxiliary coordinates:
          forecast_period                         x                  -                    -
          surface_altitude                        -                  x                    x
     Derived coordinates:
          altitude                                -                  x                    x
     Scalar coordinates:
          atmosphere_hybrid_height_coordinate: 20.000... m, bound=(0.0, 49.998...) m
          forecast_day: 1.0 Days
          model_level_number: 1
          sigma: 0.997..., bound=(1.0, 0.994...)
     Attributes:
          Conventions: CF-1.5
          STASH: m01s34i001
          label: aqum_oper
          short_name: O3
          source: Data from Met Office Unified Model
     Cell methods:
          mean: time (1 hour)

This has all the required coordinate axes. It also has two other requirements
which are in the attributes. The first of these is the 'short_name':

>>> print(gcl[0].attributes['short_name'])
O3

This is a shortened name given to represent the cube phenomenon/species, which
can be used for labelling plots, or generating filenames which are not too long.

The second requirement is the 'label':

>>> print(gcl[0].attributes['label'])
aqum_oper

This is a label given to distinguish where the data has come from, for example
a particular model, or it could be 'Control' and in another object 'Casestudy'.
These are used again for plotting and generating filenames as a means of
distinguishing between different model runs or observations.

Sites Cube List
^^^^^^^^^^^^^^^

A sites cube list is an iris cube list whose cubes contain time-series at
multiple sites. Now load some example data into the sites_cube_list to help
explain this component better:

>>> md_scl = md.load_ts(adaqcode.SAMPLE_DATADIR+'sites_cube_list/aqum_oper_1days.nc')

If we again print the ADAQData object, we can see that the sites_cube_list has
now also been filled:

>>> print(md) # doctest: +ELLIPSIS
<class '...adaq_data.ADAQData'> - Subclass of ADAQData - Contains:
sites_cube_list:
0: mass_concentration_of_ozone_in_air / (ug/m3) (site_id: 5; time: 25)
1: mass_concentration_of_pm2p5_dry_aerosol_in_air / (ug/m3) (site_id: 5; time: 25)
gridded_cube_list:
0: mass_concentration_of_ozone_in_air / (ug/m3) (time: 25; grid_latitude: 182; grid_longitude: 146)
1: mass_concentration_of_pm2p5_dry_aerosol_in_air / (ug/m3) (time: 25; grid_latitude: 182; grid_longitude: 146)
trajectory_cube_list:
< No cubes >

We can now access the this sites_cube_list and try printing the first cube in
the list:

>>> scl = md.sites_cube_list
>>> sites_cube = scl[0]
>>> print(sites_cube)  # doctest: +ELLIPSIS
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
          Conventions: CF-1.5
          STASH: m01s34i001
          label: aqum_oper
          short_name: O3
          source: Data from Met Office Unified Model
     Cell methods:
          mean: time (1 hour)

This cube contains two required coordinates: time, plus site_id. These are
iris cube coordinates.

>>> print(type(sites_cube.coord('site_id')))
<class 'iris.coords.DimCoord'>

The site_id coordinate is an identifier to distinguish uniquely between different
sites. It is a number derived from the sites latitude and longitude. However most of
the time you will not need these numbers, they are purely used for organising and
ordering the cube. We can examine the values for this dimension by looking at
the points attribute of the coordinate:

>>> import numpy as np
>>> np.set_printoptions(formatter={'float':lambda x: '{:5.2f}'.format(x)})
>>> print(sites_cube.coord('site_id').points)
[35628361.14 35665278.15 35696583.14 35790611.15 35867333.14]

On the same axis as this site_id are all the other coordinates related to the site. For
example 'latitude', 'longitude' and 'site_name':

>>> print(sites_cube.coord('latitude').points)
[50.60 55.88 52.50 57.16 51.57]
>>> np.set_printoptions()
>>> print(sites_cube.coord('site_name').points)
['Yarner_Wood' 'Auchencorth_Moss' 'Aston_Hill' 'Aberdeen' 'Harwell']

Note this cube also has grid_latitude and grid_longitude.
This comes from the model data grid that this cube was originally derived from,
so for data extracted from AQUM which is on a rotated grid. The coordinates
grid_latitude and grid_longitude therefore remain to enable direct comparison
back to the model.

A particular site can be extracted, eg Harwell, using
:func:`adaq_data.ADAQData.extract`:

>>> scl_harwell = md.extract(site_name='Harwell')
>>> print(scl_harwell.sites_cube_list)
0: mass_concentration_of_ozone_in_air / (ug/m3) (time: 25)
1: mass_concentration_of_pm2p5_dry_aerosol_in_air / (ug/m3) (time: 25)

We can set singlecube=True to ensure a cube is returned as opposed to a cubelist -
if just one cube can not be returned then an error is raised.
So an example of extracting a site and species at the same time into a cube.

>>> har_o3 = md.extract(site_name='Harwell', short_name='O3', singlecube=True)
>>> print(har_o3)  # doctest: +ELLIPSIS
mass_concentration_of_ozone_in_air / (ug/m3) (time: 25)
     Dimension coordinates:
          time                                    x
     Auxiliary coordinates:
          forecast_period                         x
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
          surface_altitude: 104.795... m
     Attributes:
          Conventions: CF-1.5
          STASH: m01s34i001
          label: aqum_oper
          short_name: O3
          source: Data from Met Office Unified Model
     Cell methods:
          mean: time (1 hour)

This now only has a time-dimension so can be useful for plotting time-series
etc.

Trajectory Cube List
--------------------
A trajectory cube list is an iris cube list whose cubes are single trajectories
containing one dimension coordinate (time) and at least three auxiliary
coordinates along the same dimension (an X, Y and Z coordinate). Now load
some example data into the trajectory cube list to help explain this
component better:

>>> md_tcl = md.load_trajectory(adaqcode.SAMPLE_DATADIR+'name_trajectory/trajectory_cube_list.nc')

If we print the ADAQData object, we can see that the trajectory_cube_list has
now been filled with two trajectories and thus two cubes (both with the variable 'pressure').
Note that plotting only one of these will produce a single line.

>>> print(md) # doctest: +ELLIPSIS
<class '...adaq_data.ADAQData'> - Subclass of ADAQData - Contains:
sites_cube_list:
0: mass_concentration_of_ozone_in_air / (ug/m3) (site_id: 5; time: 25)
1: mass_concentration_of_pm2p5_dry_aerosol_in_air / (ug/m3) (site_id: 5; time: 25)
gridded_cube_list:
0: mass_concentration_of_ozone_in_air / (ug/m3) (time: 25; grid_latitude: 182; grid_longitude: 146)
1: mass_concentration_of_pm2p5_dry_aerosol_in_air / (ug/m3) (time: 25; grid_latitude: 182; grid_longitude: 146)
trajectory_cube_list:
0: Pressure / (1)                      (time: 193)
1: Pressure / (1)                      (time: 193)

We can now access this trajectory_cube_list and try printing the first cube
in the list:

>>> tcl = md.trajectory_cube_list
>>> trajectory_cube = tcl[0]
>>> print(trajectory_cube)
Pressure / (1)                      (time: 193)
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
          source: Source_2
     Attributes:
          Conventions: CF-1.5
          label: TRAJ
          short_name: Pressure

The cube contains the required dimension coordinate:time.

>>> print(type(trajectory_cube.coord('time')))
<class 'iris.coords.DimCoord'>

We can also see that longitude, latitude and altitude are auxiliary
coordinates on the same dimension. e.g.:

>>> print(type(trajectory_cube.coord('longitude')))
<class 'iris.coords.AuxCoord'>

The cube also contains another auxiliary coordinate: Travel Time, and a number
of scalar coordinates.

Site Information Data
---------------------

Site information can be read in from ascii text files, for example:

>>> sitesfile = adaqcode.SAMPLE_DATADIR+'AURN_obs/aq_sites_GEMSRAQ_v4b_dev.txt'

The exact format of this file can be seen:

>>> with open(sitesfile,"r") as fin:
... 	for line in fin:
...         print(line) # doctest: +NORMALIZE_WHITESPACE
GEMS_code  Abbrev    Lat        Lon      Altit        Type      NAME
GB0001  ABD   57.15750122  -2.093888760   20  URBAN_BACKGROUND  Aberdeen
GB0003  ACTH  55.88333511  -3.347222328  260  RURAL             Auchencorth_Moss
GB0002  AH    52.50361252  -3.034166574  370  RURAL             Aston_Hill
GB0045  HAR   51.57110977  -1.326666594  137  RURAL             Harwell
GB0128  YW    50.59749985  -3.716388941  119  RURAL             Yarner_Wood

Note there is a line of headers, followed by rows of site information -
one row per site. The key headers that are required are 'latitude' (or 'lat'),
'longitude' (or 'lon'), and 'site_name' (or 'name') - note all headers are
converted to lowercase when read in.

This can be read in to adaqcode using :class:`sites_info.SitesInfo`:

>>> si = adaqcode.sites_info.SitesInfo()

If reading from this file, without any observations, you need to set allsites=True to
ensure that all sites are read in:

>>> sites_data = si.read_from_file(sitesfile, allsites=True)
Number of sites:  5

This returns a numpy ndarray. This could be printed directly, however
to print this data nicely, use the function :func:`sites_info.format_string`:

>>> fmtstr = adaqcode.sites_info.format_string(sites_data)
>>> print(fmtstr)
{:>6},{:>4},{:8.3f},{:7.3f},{:4d},{:>16},{:>16}

This then allows you to loop over every site and print it nicely:

>>> for site in sites_data:
...    print(fmtstr.format(*site))
GB0001, ABD,  57.158, -2.094,  20,URBAN_BACKGROUND,        Aberdeen
GB0003,ACTH,  55.883, -3.347, 260,           RURAL,Auchencorth_Moss
GB0002,  AH,  52.504, -3.034, 370,           RURAL,      Aston_Hill
GB0045, HAR,  51.571, -1.327, 137,           RURAL,         Harwell
GB0128,  YW,  50.597, -3.716, 119,           RURAL,     Yarner_Wood

The numpy ndarray, sites_data can be treated a bit like a dictionary in some
respects. For example you can get the keys (which were the headers):

>>> print(sites_data.dtype.names)
('gems_code', 'abbrev', 'latitude', 'longitude', 'site_altitude', 'site_type', 'site_name')

Which can then be used to get columns of data:

>>> print(sites_data['site_name'])
['Aberdeen' 'Auchencorth_Moss' 'Aston_Hill' 'Harwell' 'Yarner_Wood']

Or single sites can be accessed from the rows, by indexing the data:

>>> print(fmtstr.format(*sites_data[0]))
GB0001, ABD,  57.158, -2.094,  20,URBAN_BACKGROUND,        Aberdeen

Note if you have a directory of observations and want to limit the returned
list to only sites where observations are available, then can use:

>>> sites_data = si.read_from_file(sitesfile, allsites=False,
... obsdir=adaqcode.SAMPLE_DATADIR+'AURN_obs/')
Number of sites:  5

Note a simpler routine to return the sites_data array, using .ini files
which is described below (:ref:`user_guide_inifiles`),
:func:`sites_info.get_siteinfo`.

Extracting a sites_cube_list from a gridded_cube_list
-----------------------------------------------------

Model data is generally produced as gridded data which can be read in easily
to a gridded_cube_list. Observations however are often produced at specific
sites so are read into a sites_cube_list. For comparision between the two,
model data should be extracted from the gridded data to site specific data.
This is done using bilinear interpolation which then gives us a sites_cube_list.
This interpolation is only done in the X and Y coordinates (not on T, Z or any
other coordinate).
This is easily done given a gridded_cube_list in an ADAQData object and site
information in the form of an ndarray from a :class:`sites_info.SitesInfo`
object.

Firstly get our ADAQData object and our sites data:

>>> md = adaqcode.adaq_data.ADAQData()
>>> md_gcl = md.load_gridded(adaqcode.SAMPLE_DATADIR+'gridded_cube_list/aqum_oper_1days.nc')
>>> print(md) # doctest: +ELLIPSIS
<class '...adaq_data.ADAQData'> - Subclass of ADAQData - Contains:
sites_cube_list:
< No cubes >
gridded_cube_list:
0: mass_concentration_of_ozone_in_air / (ug/m3) (time: 25; grid_latitude: 182; grid_longitude: 146)
1: mass_concentration_of_pm2p5_dry_aerosol_in_air / (ug/m3) (time: 25; grid_latitude: 182; grid_longitude: 146)
trajectory_cube_list:
< No cubes >

>>> si = adaqcode.sites_info.SitesInfo()
>>> sites_data = si.read_from_file(sitesfile, allsites=True)
Number of sites:  5

Now extract these sites. This particular extraction routine needs to know the
names of the X coordinate and the Y coordinate:

>>> md_scl = md.extract_scl_from_gridded(sites_data, "grid_longitude", "grid_latitude")
>>> print(md) # doctest: +ELLIPSIS
<class '...adaq_data.ADAQData'> - Subclass of ADAQData - Contains:
sites_cube_list:
0: mass_concentration_of_ozone_in_air / (ug/m3) (site_id: 5; time: 25)
1: mass_concentration_of_pm2p5_dry_aerosol_in_air / (ug/m3) (site_id: 5; time: 25)
gridded_cube_list:
0: mass_concentration_of_ozone_in_air / (ug/m3) (time: 25; grid_latitude: 182; grid_longitude: 146)
1: mass_concentration_of_pm2p5_dry_aerosol_in_air / (ug/m3) (time: 25; grid_latitude: 182; grid_longitude: 146)
trajectory_cube_list:
< No cubes >

Subclasses of ADAQData which have gridded data, also know what their X and Y
coordinates are, eg PPData, so can use their simpler routine which will be
called extract_sites (eg :func:`pp_data.PPData.extract_sites` or
:func:`name_data.NAMEData.extract_sites`)

>>> md = adaqcode.pp_data.PPData()
>>> md_gcl = md.load_gridded(adaqcode.SAMPLE_DATADIR+'gridded_cube_list/aqum_oper_1days.nc')
>>> md_scl = md.extract_sites(sites_data)
>>> print(md) # doctest: +ELLIPSIS
<class '...pp_data.PPData'> - Subclass of ADAQData - Contains:
sites_cube_list:
0: mass_concentration_of_ozone_in_air / (ug/m3) (site_id: 5; time: 25)
1: mass_concentration_of_pm2p5_dry_aerosol_in_air / (ug/m3) (site_id: 5; time: 25)
gridded_cube_list:
0: mass_concentration_of_ozone_in_air / (ug/m3) (time: 25; grid_latitude: 182; grid_longitude: 146)
1: mass_concentration_of_pm2p5_dry_aerosol_in_air / (ug/m3) (time: 25; grid_latitude: 182; grid_longitude: 146)
trajectory_cube_list:
< No cubes >

.. _user_guide_inifiles:

Inifiles
--------

Inifiles are designed to allow a user to set up values such
as dates, filenames, directories etc for use within a script, without having
to change any code each time you want to do something slightly different.

Inifiles are simple ascii text files, which look like:

.. literalinclude:: ../../adaqcode/inifile.ini

These have keys given on the left hand side of an equals sign,
and values given on the right. Any lines beginning with # are comments.

These can be read in by setting up an :class:`inifile.Inifile` class
and passing in the filename of your inifile:

>>> inifilename = adaqcode.CODE_DIR + 'adaqcode/inifile.ini'
>>> ini_dict = adaqcode.inifile.Inifile(inifilename) # doctest: +ELLIPSIS
Reading inifile .../inifile.ini

This reads the file and puts all the keys/values into a dictionary:

>>> print('Keys:'); print(sorted(ini_dict.keys())) # doctest: +NORMALIZE_WHITESPACE
Keys:
['aerosol_units', 'calc_stats', 'calc_stats_format_list', 'chem_units', 
'contours', 'daqi', 'diurnal', 'end_date', 'end_datetime', 'extent_list', 
'field_attribute_dict', 'forecast_day', 'histograms', 'html_dir', 
'levels_list', 'line_colours_list', 'models_dir_list', 'models_fmt_list', 
'models_list', 'obs_dir', 'obs_dir_list', 'obs_fmt', 'obsdir_fixed_csv', 
'plot_dir', 'pollen_units', 'qqplots', 'range_days', 'rolling_stats_file_list',
'scaling_factors_list', 'short_name_list', 'site_types_list', 'sites_file', 
'soccer_plots', 'start_date', 'start_datetime', 'strict_statistics', 
'timeseries', 'timeseries_multiple_short_names', 
'timeseries_multiple_short_names_dict', 'timeseries_of_stats', 
'timeseries_of_stats_list']

Most values are left as their given type, for example:

>>> print(ini_dict['obs_fmt'])
aurn

But any True, False or None values are converted to python True and False
and None:

>>> print(ini_dict['timeseries'])
True

Also, any keys whose name ends in '_list' are converted to python lists:

>>> print(ini_dict['models_fmt_list'])
['pp', 'maccens']

Note the individual values are read as their correct type as well, for example:

>>> print(ini_dict['extent_list'])
[120, 240, -20, 80]

And any keys whose name ends in '_dict' are converted to python dictionaries:

>>> print(sorted(ini_dict['field_attribute_dict'].keys()))
['Quantity', 'Species']
>>> print(ini_dict['field_attribute_dict']['Species'])
PM10

Note an alternative routine that can be used to do the same thing,
:func:`inifile.get_inidict`:

>>> ini_dict = adaqcode.inifile.get_inidict(inifilename) # doctest: +ELLIPSIS
Reading inifile .../inifile.ini

As mentioned in the section above, a simpler routine exists for returning
the sites_data array, given the ini_dict. This uses the 'sites_file'
given in the inifile, alongwith possibly also 'obs_dir' and 'sites_types_list'
(which limits the returned list based on 'site_type'):

>>> sites_data = adaqcode.sites_info.get_siteinfo(ini_dict)
Number of sites:  5

.. Note to developers. If you change any of the questions below then you will
   also need to change the answers in answers.rst


.. _intro_questions_ref:

Introduction Exercises
----------------------

  1. Print the short name and the label of the 5th gridded cube in the
     sample data file called *aqum_casestudy_5days.nc*

     a. Load the data file *aqum_casestudy_5days.nc* which is in the
        *gridded_cube_list* directory in the sample data directory.
        (Hint: the sample data directory is *adaqcode.SAMPLE_DATADIR*)
     b. Select the 5th gridded cube and print the short name and the label

  2. Now load in the sites_cube_list *aqum_oper_1days.nc* from the
     *sites_cube_list* directory in the sample data directory and print
     the abbreviated names of all the sites in the first sites cube list.
     (Hint: The coordinate of abbreviated names is called *abbrev*).

  3. Read in the sites files *ctbto_sites.txt* in the *ctbto_obs* subfolder in
     the sample data directory and print out the sites data.

  4. Read in the inifile *name_field_plot.ini* which is in the *adaqscripts*
     directory in the ADAQ code. Print the *models_fmt_list* and the
     *field_attribute_dict*. (Hint: the code directory is *adaqcode.CODE_DIR*)

Answers to :ref:`intro_answers_ref`
