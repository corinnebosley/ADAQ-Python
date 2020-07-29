"""
Module to hold ADAQData class.
"""
from __future__ import division
from __future__ import print_function

from six.moves.builtins import str
from six.moves.builtins import zip
from six.moves.builtins import object

import os
import copy
import warnings

import numpy as np
import iris

import cube_statistics

from distutils.version import LooseVersion
if LooseVersion(iris.__version__) <= '1.13.0':
    iris.FUTURE.cell_datetime_objects = True


def generate_siteids(lons, lats):
    """
    Generate list of siteids, by combining longitudes and latitudes.

    :param lons: list of longitude values
    :param lats: list of latitude values

    Assumes -90.<=lats<=+90.

    Returns a site id which looks like lon.lat,
    where lat is in range 0<=lat<180,
    and lon is in the range 0.<=lon<360.

    Rounds both lat and lon to 5 decimal places. So note, that
    any two sites which match both lat and lon to 5 dp will return
    the same site_id.

    For example:

    >>> lon=10.
    >>> lat=50.
    >>> siteids = generate_siteids([lon],[lat])
    >>> print("%18.8f" % siteids[0])
      1000000.14000000

    >>> lons=[-25.123456,0.0,367.12]
    >>> lats=[-10.123456,0.0,89.999]
    >>> siteids=generate_siteids(lons,lats)
    >>> print("%18.8f"*3 % tuple(siteids))
     33487654.07987654        0.09000000   712000.17999900
    """

    siteids = [float(int(np.round((lon % 360.) * 1e5)))
               + float(int(np.round((90. + lat) * 1e5))) / 1.e8
               for lon, lat in zip(lons, lats)]

    return siteids

def site_cubes_comparable(sitecube1, sitecube2):
    """
    Checks that sitecubes are comparable.
    In this case 'comparable' means that:

      *  site_id dimension contains same values for coordinates
      *  site_name dimension contains same names for coordinate points

    >>> import config
    >>> sample_datadir = config.SAMPLE_DATADIR+'sites_cube_list/'
    >>> md = ADAQData()
    >>> md_scl = md.load_ts(sample_datadir+'aqum_oper_1days.nc')
    >>> site_cubes_comparable(md_scl[0],md_scl[1])

    Should complain if they haven't got the same number of sites:

    >>> sitecube1 = md.extract(short_name='PM2p5', singlecube=True)
    >>> sitecube2 = md.extract(short_name='PM2p5',
    ... singlecube=True, site_type='RURAL')
    >>> site_cubes_comparable(sitecube1,sitecube2) # doctest: +ELLIPSIS
    Traceback (most recent call last):
      File ".../doctest.py", line ..., in __run
        compileflags, 1) in test.globs
      File "<doctest __main__.ADAQData.site_cubes_comparable[...]>", \
line 1, in <module>
        md.site_cubes_comparable(sitecube1,sitecube2) # doctest: +ELLIPSIS
      File ".../adaq_data.py", line ..., in site_cubes_comparable
        raise ValueError,"Site lists do not match: different lengths"
    ValueError: Site lists do not match: different lengths

    Or have different site names:

    >>> sitecube3 = md.extract(short_name='O3',singlecube=True)
    >>> sitecube3.coord('site_name').points[0]='Random_Site'
    >>> site_cubes_comparable(sitecube1,sitecube3) # doctest: +ELLIPSIS
    Traceback (most recent call last):
      File ".../doctest.py", line ..., in __run
        compileflags, 1) in test.globs
      File "<doctest __main__.ADAQData.site_cubes_comparable[...]>", \
line 1, in <module>
        md.site_cubes_comparable(sitecube1,sitecube3)
      File ".../adaq_data.py", line ..., in site_cubes_comparable
        + coord_name + ') ' + v1 + ' != ' + v2
    ValueError: Site lists do not match: (site_name) Yarner_Wood != Random_Site
    """

    if not len(sitecube1.coord('site_id').points) \
           == len(sitecube2.coord('site_id').points):

        raise ValueError("Site lists do not match: different lengths")

    for coord_name in ['site_name', 'site_id']:

        for v1, v2 in zip(sitecube1.coord(coord_name).points,
                          sitecube2.coord(coord_name).points):

            try:
                assert v1 == v2
            except:
                raise ValueError("Site lists do not match: (" \
                      + coord_name + ') ' + v1 + ' != ' + v2)

def get_derived_sn_specs(short_name_list, derived_short_name_dict=None):
    """
    Function to identify short names that need to be derived from others.

    :param short_name_list: list of required short names.
    :param derived_short_name_dict:
        (optional) dictionary whose keys are short names, and whose values
        are specifications (see below).

    :returns: **specifications** - a list with the same length as
              `short_name_list`. Each element is a dictionary for the
              corresponding short name, specifying how to calcuate it.
              See :any:`ADAQData.derive` for details on the format.

    In the first instance, this function simply looks up each short name
    in `derived_short_name_dict`.

    >>> spec = get_derived_sn_specs(
    ...     ["NOx"],
    ...     derived_short_name_dict={"NOx": {"components": ["NO", "NO2"]}})[0]
    >>> for key in sorted(spec):
    ...     print(key+":", spec[key]) # doctest: +ELLIPSIS
    components: ['NO', 'NO2']
    short_name: NOx

    This parameter therefore allows the function to be easily extended,
    for example the :any:`ADAQData.get_derived_sn_specs` method.

    Additionally, several simple statistics can be automatically
    recognised. For example, "mean_O3" is interpreted as a request
    to collapse the time dimension of each ozone cube using
    `iris.analysis.MEAN`:

    >>> spec = get_derived_sn_specs(["mean_O3"])[0]
    >>> for key in sorted(spec):
    ...     print(key+":", spec[key]) # doctest: +ELLIPSIS
    components: ['O3']
    function: <function ...>
    kwargs: {'aggregator': <iris.analysis.Aggregator object ...>}
    short_name: mean_O3

    Recognised statistics are currently:

    * mean
    * sum
    * max
    * min

    Optionally, one of the following time periods can also be specified:

    * day (or daily)
    * month (or monthly)
    * year (or yearly)

    These words can appear before or after the short name they depend on,
    and must separated by underscores. If the original name contains
    underscores, a word can also be inserted at any such point.
    When specifying both a statistic and a time period, there is no
    restriction on which order.

    For example, "total_precip_daily_sum", "daily_total_precip_sum",
    and "total_sum_precip_daily" are all valid ways of requesting the sum
    of "total_precip", though some read better than others:

    >>> specs = get_derived_sn_specs(["total_precip_daily_sum",
    ...                               "total_sum_precip_daily"])
    >>> specs[0]["components"] == specs[1]["components"]
    True
    >>> specs[0]["function"] == specs[1]["function"]
    True
    >>> specs[0]["kwargs"] == specs[1]["kwargs"]
    True
    """

    #default to an empty dict
    if derived_short_name_dict is None:
        derived_short_name_dict = {}

    stats = {'mean', 'sum', 'max', 'min'}
    periods = {
        'hour': 'hour', 'hourly': 'hour',
        'day': 'day', 'daily': 'day',
        'month': 'month', 'monthly': 'month',
        'year': 'year', 'yearly': 'year',
        }

    specifications = []

    for short_name in short_name_list:
        #check provided dictionary
        if short_name in derived_short_name_dict:
            spec = derived_short_name_dict[short_name]
            spec["short_name"] = short_name
            specifications.append(spec)
            continue

        #check for a named statistic
        spec = {}
        kwargs = {}
        stat = None
        period = None
        words = []
        for word in short_name.split("_"):
            if word.lower() in stats:
                stat = word.lower()
            elif word.lower() in periods:
                period = periods[word.lower()]
            else:
                words.append(word)

        spec["components"] = ["_".join(words)]

        #if no special words are recognised, assume that the name was
        # already a known short name to be loaded directly
        if stat is None:
            specifications.append({'components': [short_name]})
            continue

        #choose a function and arguments to carry out the calculation
        if period is None:
            stat = "NAN"+stat.upper()
            kwargs["aggregator"] = cube_statistics.CUBE_AGGREGATORS[stat]
            spec["function"] = (lambda cube, aggregator:
                                cube.copy().collapsed("time", aggregator))
        else:
            spec["function"] = cube_statistics.periodic_stat
            kwargs["stat"] = stat
            kwargs["period"] = period

        spec["short_name"] = short_name
        if kwargs:
            spec["kwargs"] = kwargs
        specifications.append(spec)

    return specifications


class ADAQData(object):
    """
    Class to hold cubes from model or observation data.
    Will contain one or more of a sites_cube_list, a gridded_cube_list and
    a trajectory_cube_list

    To initialise an ADAQData object:

    >>> md = ADAQData()
    >>> print(md) # doctest: +ELLIPSIS
    <class '....ADAQData'> ... Contains:
    sites_cube_list:
    < No cubes >
    gridded_cube_list:
    < No cubes >
    trajectory_cube_list:
    < No cubes >

    By default the sites_cube_list, gridded_cube_list and trajectory_cube_list
    are all set to an empty iris cubelist (iris.cube.CubeList) of length zero:

    >>> print(md.sites_cube_list)
    < No cubes >
    >>> print(type(md.gridded_cube_list))
    <class 'iris.cube.CubeList'>
    >>> print(len(md.gridded_cube_list))
    0
    >>> print(md.trajectory_cube_list)
    < No cubes >

    **ADAQData.sites_cube_list**

    This is designed to hold time-series data at specific sites.
    This should be an iris cube list, where each cube in the list has certain
    attributes and coordinates:

    >>> print(ADAQData.SITES_CUBE_REQUIRED_ATTRIBUTES)
    ['short_name', 'label']
    >>> print(ADAQData.SITES_CUBE_REQUIRED_COORDS)
    ['latitude', 'longitude', 'site_name']
    >>> print(ADAQData.SITES_CUBE_REQUIRED_DIMCOORDS)
    ['time', 'site_id']

    Each site should have a name, latitude and longitude associated with it,
    as a minimum. Other coordinates, such as a site abbreviation or altitude
    could also be given, but these are not compulsory.

    Load time-series data into sites_cube_list:

    >>> import config
    >>> sample_datadir = config.SAMPLE_DATADIR+'sites_cube_list/'
    >>> md_scl = md.load_ts(sample_datadir+'aqum_oper_1days.nc')

    View sites_cube_list:

    >>> print(md.sites_cube_list)
    0: mass_concentration_of_ozone_in_air / (ug/m3) \
(site_id: 5; time: 25)
    1: mass_concentration_of_pm2p5_dry_aerosol_in_air / (ug/m3) \
(site_id: 5; time: 25)

    View a single sites_cube within this list:

    >>> print(md.sites_cube_list[0]) # doctest: +ELLIPSIS
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
              level_height: 20.00... m, bound=(0.0, 49.99...) m
              model_level_number: 1
              sigma: 0.99..., bound=(1.0, 0.99...)
         Attributes:
              Conventions: CF-1.5
              STASH: m01s34i001
              label: aqum_oper
              short_name: O3
              source: Data from Met Office Unified Model
         Cell methods:
              mean: time (1 hour)

    Check this is a valid sites_cube_list using check_sites_cube(),
    or pass in a single sites_cube:

    >>> md.check_sites_cube(sites_cube=md.sites_cube_list[0])
    short_name - OK
    label - OK
    latitude  coord - OK
    longitude  coord - OK
    site_name  coord - OK
    time  dim_coord - OK
    site_id  dim_coord - OK

    **ADAQData.gridded_cube_list**

    This is designed to hold gridded fields data.
    This should be an iris cube list, where each cube in the list has certain
    attributes and coordinates:

    >>> print(ADAQData.GRIDDED_CUBE_REQUIRED_ATTRIBUTES)
    ['short_name', 'label']
    >>> print(ADAQData.GRIDDED_CUBE_REQUIRED_COORD_AXES)
    ['T', 'X', 'Y']

    Note this should contain a time coordinate. Other coordinates,
    such as height are allowed, but are not compulsory.

    Load gridded data into a gridded_cube_list:

    >>> sample_datadir = config.SAMPLE_DATADIR+'gridded_cube_list/'
    >>> md_gcl = md.load_gridded(sample_datadir+'aqum_oper_1days.nc')
    >>> print(md_gcl)
    0: mass_concentration_of_ozone_in_air / (ug/m3) \
(time: 25; grid_latitude: 182; grid_longitude: 146)
    1: mass_concentration_of_pm2p5_dry_aerosol_in_air / (ug/m3) \
(time: 25; grid_latitude: 182; grid_longitude: 146)
    >>> md.check_gridded_cube()
    short_name - OK
    label - OK
    T axis - OK - coord= time
    X axis - OK - coord= grid_longitude
    Y axis - OK - coord= grid_latitude
    short_name - OK
    label - OK
    T axis - OK - coord= time
    X axis - OK - coord= grid_longitude
    Y axis - OK - coord= grid_latitude

    View a single gridded cube within the gridded_cube_list:

    >>> print(md_gcl[0]) # doctest: +ELLIPSIS
    mass_concentration_of_ozone_in_air / (ug/m3) \
(time: 25; grid_latitude: 182; grid_longitude: 146)
         Dimension coordinates:
              time                                    \
x                  -                    -
              grid_latitude                           \
-                  x                    -
              grid_longitude                          \
-                  -                    x
         Auxiliary coordinates:
              forecast_period                         \
x                  -                    -
              surface_altitude                        \
-                  x                    x
         Derived coordinates:
              altitude                                \
-                  x                    x
         Scalar coordinates:
              atmosphere_hybrid_height_coordinate: \
20.00... m, bound=(0.0, 49.99...) m
              forecast_day: 1.0 Days
              model_level_number: 1
              sigma: 0.99..., bound=(1.0, 0.99...)
         Attributes:
              Conventions: CF-1.5
              STASH: m01s34i001
              label: aqum_oper
              short_name: O3
              source: Data from Met Office Unified Model
         Cell methods:
              mean: time (1 hour)

    >>> print(md) # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    <class '....ADAQData'> ... Contains:
    sites_cube_list:
    0: mass_concentration_of_ozone_in_air / (ug/m3) \
(site_id: 5; time: 25)
    1: mass_concentration_of_pm2p5_dry_aerosol_in_air / (ug/m3) \
(site_id: 5; time: 25)
    gridded_cube_list:
    0: mass_concentration_of_ozone_in_air / (ug/m3) \
(time: 25; grid_latitude: 182; grid_longitude: 146)
    1: mass_concentration_of_pm2p5_dry_aerosol_in_air / (ug/m3) \
(time: 25; grid_latitude: 182; grid_longitude: 146)
    trajectory_cube_list:
    < No cubes >

    **ADAQData.trajectory_cube_list**

    This is designed to hold trajectory data with each trajectory held as
    a single cube in a cube list.
    This should be an iris cube list, where each cube in the list has certain
    attributes and coordinates:

    >>> print(ADAQData.TRAJECTORY_CUBE_REQUIRED_ATTRIBUTES)
    ['short_name', 'label']
    >>> print(ADAQData.TRAJECTORY_CUBE_REQUIRED_COORD_AXES)
    ['X', 'Y']

    Note this should contain a time coordinate and a height coordinate.
    Other coordinates, such as time of arrival are allowed, but are not
    compulsory.

    Load trajectory data into a trajectory_cube_list:

    >>> sample_datadir = config.SAMPLE_DATADIR+'name_trajectory/'
    >>> md_tcl = md.load_trajectory(sample_datadir+'trajectory_cube_list.nc')
    >>> print(md_tcl)
    0: Pressure / (1)                      (time: 193)
    1: Pressure / (1)                      (time: 193)

    View a single cube within the trajectory cube list:

    >>> traj_cube = md_tcl.extract(iris.Constraint(Source='Source_3'),
    ... strict=True)
    >>> print(traj_cube)
    Pressure / (1)                      (time: 193)
         Dimension coordinates:
              time                           x
         Auxiliary coordinates:
              Travel Time                    x
              altitude                       x
              latitude                       x
              longitude                      x
         Scalar coordinates:
              PP Index: 3.0
              Release Time: 2016-12-02 08:00:00
              Source: Source_3
         Attributes:
              Conventions: CF-1.5
              label: TRAJ
              short_name: Pressure

    """

    SITES_CUBE_REQUIRED_ATTRIBUTES = ['short_name', 'label']
    SITES_CUBE_REQUIRED_COORDS = ['latitude', 'longitude', 'site_name']
    SITES_CUBE_REQUIRED_DIMCOORDS = ['time', 'site_id']
    GRIDDED_CUBE_REQUIRED_ATTRIBUTES = ['short_name', 'label']
    GRIDDED_CUBE_REQUIRED_COORD_AXES = ['T', 'X', 'Y']
    # Note that TRAJECTORY_CUBE_REQUIRED_COORD_AXES should be X Y Z but iris
    # does not currently recognise NAME vertical coordinates as Z axes.
    TRAJECTORY_CUBE_REQUIRED_COORD_AXES = ['X', 'Y']
    TRAJECTORY_CUBE_REQUIRED_DIMCOORDS = ['time']
    TRAJECTORY_CUBE_REQUIRED_ATTRIBUTES = ['short_name', 'label']

    DERIVED_SHORT_NAME_DICT = {}

    def __init__(self, sites_cube_list=None, gridded_cube_list=None,
                 trajectory_cube_list=None):
        """
        Initialise the ADAQ Data class.
        """
        if sites_cube_list is None:
            sites_cube_list = iris.cube.CubeList()
        if gridded_cube_list is None:
            gridded_cube_list = iris.cube.CubeList()
        if trajectory_cube_list is None:
            trajectory_cube_list = iris.cube.CubeList()
        self.sites_cube_list = sites_cube_list
        self.gridded_cube_list = gridded_cube_list
        self.trajectory_cube_list = trajectory_cube_list

    def __str__(self):
        """
        Output from print
        """
        string = str(type(self))
        if string != "<class '__main__.ADAQData'>":
            string += ' - Subclass of ADAQData'
        string += ' - Contains:\nsites_cube_list:'
        string += '\n' + iris.cube.CubeList.__str__(self.sites_cube_list)
        string += '\ngridded_cube_list:'
        string += '\n' + iris.cube.CubeList.__str__(self.gridded_cube_list)
        string += '\ntrajectory_cube_list:'
        string += '\n' + iris.cube.CubeList.__str__(self.trajectory_cube_list)
        return string

    def __check_attribute(self, cube, attrib_name, verbose=True):
        """
        Internal function to check all required attributes exist
        """
        if attrib_name in cube.attributes:
            if verbose:
                print(attrib_name, '- OK')
        else:
            raise ValueError(attrib_name + ' is not in attributes')


    def __check_coord(self, cube, coord_name):
        """
        Internal function to check all required coordinates exist
        """
        if cube.coords(coord_name):
            print(coord_name, ' coord - OK')
        else:
            raise ValueError(coord_name + ' coord does not exist')

    def __check_dim_coord(self, cube, coord_name):
        """
        Internal function to check all required dim coordinates exist
        """
        dim_coord_names = [coord.name() for coord in cube.dim_coords]
        if coord_name in dim_coord_names:
            print(coord_name, ' dim_coord - OK')
        else:
            raise ValueError(coord_name + ' dim_coord does not exist')

    def check_sites_cube(self, sites_cube=None):
        """
        Checks a single sitecube if given (otherwise checks all in
        self.sites_cube_list)
        to determine if it has all required attributes and coordinates.
        Should also be an iris cube, so by default has a phenomenon name
        and units
        Raises an error if any missing.
        """

        if sites_cube is not None:
            cubelist = iris.cube.CubeList([sites_cube])
        else:
            cubelist = self.sites_cube_list

        for sites_cube in cubelist:

            if not isinstance(sites_cube, iris.cube.Cube):
                raise ValueError('Item in list is not an iris cube')

            for attrib_name in self.SITES_CUBE_REQUIRED_ATTRIBUTES:
                self.__check_attribute(sites_cube, attrib_name)

            for coord_name in self.SITES_CUBE_REQUIRED_COORDS:
                self.__check_coord(sites_cube, coord_name)

            for coord_name in self.SITES_CUBE_REQUIRED_DIMCOORDS:
                self.__check_dim_coord(sites_cube, coord_name)

    def check_gridded_cube(self, gridded_cube=None, gridded_cube_list=None,
                           verbose=True):
        """
        Checks a single gridded cube if given (otherwise checks
        all in self.gridded_cube_list)
        to determine if it has all required attributes and coordinates.
        Should also be an iris cube, so by default has a phenomenon name
        and units
        Raises an error if any missing.
        """

        if gridded_cube is not None:
            cubelist = iris.cube.CubeList([gridded_cube])
        elif gridded_cube_list is not None:
            cubelist = gridded_cube_list
        else:
            cubelist = self.gridded_cube_list

        for gridded_cube in cubelist:

            if not isinstance(gridded_cube, iris.cube.Cube):
                raise ValueError('Item in list is not an iris cube')

            for attrib_name in self.GRIDDED_CUBE_REQUIRED_ATTRIBUTES:
                self.__check_attribute(gridded_cube, attrib_name, verbose)


            for axis in self.GRIDDED_CUBE_REQUIRED_COORD_AXES:
                found_axis = False
                for coord in gridded_cube.coords():
                    if iris.util.guess_coord_axis(coord) == axis:
                        if verbose:
                            print(axis, 'axis - OK - coord=', coord.name())
                        found_axis = True
                        if axis == 'X' or axis == 'Y':
                            if not coord in gridded_cube.dim_coords:
                                found_axis = False
                        break

                if not found_axis:
                    raise ValueError(axis + ' axis coord does not exist')


    def check_trajectory_cube(self, trajectory_cube=None,
                              trajectory_cube_list=None, verbose=True):
        """
        Checks a single trajectory cube if given (otherwise checks
        all in self.trajectory_cube_list)
        to determine if it has all required attributes and coordinates.
        Should also be an iris cube, so by default has a phenomenon name
        and units
        Raises an error if any missing.
        """

        if trajectory_cube is not None:
            cubelist = iris.cube.CubeList([trajectory_cube])
        elif trajectory_cube_list is not None:
            cubelist = trajectory_cube_list
        else:
            cubelist = self.trajectory_cube_list

        for trajectory_cube in cubelist:

            if not isinstance(trajectory_cube, iris.cube.Cube):
                raise ValueError('Item in list is not an iris cube')

            for attrib_name in self.TRAJECTORY_CUBE_REQUIRED_ATTRIBUTES:
                self.__check_attribute(trajectory_cube, attrib_name, verbose)


            for axis in self.TRAJECTORY_CUBE_REQUIRED_COORD_AXES:
                found_axis = False
                for coord in trajectory_cube.coords():
                    if iris.util.guess_coord_axis(coord) == axis:
                        if verbose:
                            print(axis, 'axis - OK - coord=', coord.name())
                        found_axis = True
                        break

                if not found_axis:
                    raise ValueError(axis + ' axis coord does not exist')

            for coord_name in self.TRAJECTORY_CUBE_REQUIRED_DIMCOORDS:
                self.__check_dim_coord(trajectory_cube, coord_name)


    def extract(self, singlecube=False, return_cube_list=False,
                short_name=None, label=None, species=None,
                quantity=None, gridded=False, **kwargs):
        """
        Extracts from sites_cube_list (default) using iris.Constraint keywords.

        Returns same class of object as passed in (apart from the
        exceptions given according to the parameters below).

        :param singlecube: Set to True to return a single cube. Note \
                          this raises an error if more than one cube in \
                          cubelist. If there are no cubes to return, then
                          None is returned instead.

        :param return_cube_list: Set to True to return the cube list instead \
                                 of an ADAQData object.

        :param short_name: Set to a short name required to be extracted. \
                           This is treated differently as it is a cube \
                           attribute.

        :param label: Set to a label required to be extracted. \
                      This is treated differently as it is a cube attribute.

        :param species: Set to a species required to extracted. \
                      This is treated differently as it is a cube attribute.

        :param quantity: Set to a quantity required to extracted. \
                      This is treated differently as it is a cube attribute.

        :param gridded: Set to True to modify/return the gridded_cube_list \
                        instead of the default sites_cube_list.



        For more complicated extractions, instead use:
         self.sites_cube_list.extract(iris.Constraint(...))

        >>> import config
        >>> sample_datadir = config.SAMPLE_DATADIR+'sites_cube_list/'
        >>> OD = ADAQData()
        >>> od_slc = OD.load_ts(sample_datadir+'aurn_1days.nc')
        >>> Rural_OD = OD.extract(site_name='Harwell')
        >>> print(Rural_OD.sites_cube_list[0]) # doctest: +ELLIPSIS
        mass_concentration_of_ozone_in_air / (ug/m3) (time: 25)
             Dimension coordinates:
                  time                                    x
             Scalar coordinates:
                  abbrev: HAR
                  latitude: 51.57110977 degrees
                  longitude: -1.326666594 degrees
                  site_altitude: 137 m
                  site_id: 35867333.14...
                  site_name: Harwell
                  site_type: RURAL
             Attributes:
                  Conventions: CF-1.5
                  label: Obs
                  short_name: O3
                  source: AURN
             Cell methods:
                  mean: time (1 hour)

        Example of extracting from the gridded_cube_list:

        >>> sample_datadir =  config.SAMPLE_DATADIR+'gridded_cube_list/'
        >>> od_gcl = OD.load_gridded(sample_datadir+'aqum_oper_1days.nc')
        >>> o3_OD = OD.extract(short_name='O3',gridded=True)
        >>> print(o3_OD.gridded_cube_list)
        0: mass_concentration_of_ozone_in_air / (ug/m3) \
(time: 25; grid_latitude: 182; grid_longitude: 146)
        """

        newself = copy.copy(self)

        if gridded:
            input_cube_list = self.gridded_cube_list
        else:
            if not self.sites_cube_list:
                warnings.warn("sites_cube_list has not been set")
            input_cube_list = self.sites_cube_list

        if short_name is not None:
            input_cube_list = input_cube_list.extract(
                iris.AttributeConstraint(short_name=short_name))
        if label is not None:
            input_cube_list = input_cube_list.extract(
                iris.AttributeConstraint(label=label))

        if species is not None:
            input_cube_list = input_cube_list.extract(
                iris.AttributeConstraint(Species=species))

        if quantity is not None:
            input_cube_list = input_cube_list.extract(
                iris.AttributeConstraint(Quantity=quantity))

        cubelist = input_cube_list.extract(iris.Constraint(**kwargs))

        if singlecube:
            if return_cube_list:
                warnings.warn("singlecube=True and cubelist=True, will"
                              "return a singlecube (may be None).")
            #Check if single cube to return, otherwise raise error
            if not cubelist:
                cube = None
            elif len(cubelist) == 1:
                cube = cubelist[0]
            else:
                raise ValueError("Unable to return single cube - more than \
                                 one element in extracted_cubelist\n \
                                 Try running with singlecube=False")
            return cube

        elif return_cube_list:
            return cubelist

        else:
            if gridded:
                newself.gridded_cube_list = cubelist
            else:
                newself.sites_cube_list = cubelist
            return newself


    def extract_scl_from_gridded(self, sites_data, xcoord_name, ycoord_name):
        """
        Extracts site-specific data on a lat-lon grid, aka sites_cube_list (scl)
        from the gridded cube list, which may be on a different coordinate
        system.

        :param sites_data: site information from :class:`SiteInfo`
        :param xcoord_name: coordinate name for x dimension, eg 'grid_longitude'
        :param ycoord_name: coordinate name for y dimension, eg 'grid_latitude'

        Return sites_cube_list (also set as self.sites_cube_list)

        >>> import config
        >>> sample_datadir = config.SAMPLE_DATADIR+'gridded_cube_list/'
        >>> MD = ADAQData()
        >>> md_gcl = MD.load_gridded(sample_datadir+'aqum_oper_1days.nc')
        >>> print(md_gcl[0]) # doctest: +ELLIPSIS
        mass_concentration_of_ozone_in_air / (ug/m3) \
(time: 25; grid_latitude: 182; grid_longitude: 146)
             Dimension coordinates:
                  time                                    \
x                  -                    -
                  grid_latitude                           \
-                  x                    -
                  grid_longitude                          \
-                  -                    x
             Auxiliary coordinates:
                  forecast_period                         \
x                  -                    -
                  surface_altitude                        \
-                  x                    x
             Derived coordinates:
                  altitude                                \
-                  x                    x
             Scalar coordinates:
                  atmosphere_hybrid_height_coordinate: \
20.00... m, bound=(0.0, 49.99...) m
                  forecast_day: 1.0 Days
                  model_level_number: 1
                  sigma: 0.99..., bound=(1.0, 0.99...)
             Attributes:
                  Conventions: CF-1.5
                  STASH: m01s34i001
                  label: aqum_oper
                  short_name: O3
                  source: Data from Met Office Unified Model
             Cell methods:
                  mean: time (1 hour)

        Read in some sites data:

        >>> import sites_info
        >>> sitesfilename = config.SAMPLE_DATADIR+'AURN_obs/aq_sites_DEFRA.txt'
        >>> obsdir = config.SAMPLE_DATADIR+'AURN_obs/'
        >>> sites = sites_info.SitesInfo()
        >>> sites_data = sites.read_from_file(sitesfilename,
        ...         allsites=False, obsdir=obsdir)
        Number of sites:  15

        Can now extract from the gridded data, site cubes at these sites.
        Note, to do this, need to know the names of the x and y coordinates

        >>> scl = MD.extract_scl_from_gridded(sites_data,
        ... 'grid_longitude','grid_latitude')
        >>> print(MD.sites_cube_list)
        0: mass_concentration_of_ozone_in_air / (ug/m3) \
(site_id: 15; time: 25)
        1: mass_concentration_of_pm2p5_dry_aerosol_in_air / (ug/m3) \
(site_id: 15; time: 25)

        The site information has been added to the cubes:

        >>> print(MD.sites_cube_list[0]) # doctest: +ELLIPSIS
        mass_concentration_of_ozone_in_air / (ug/m3) (site_id: 15; time: 25)
             Dimension coordinates:
                  site_id                                    x         -
                  time                                       -         x
             Auxiliary coordinates:
                  abbrev                                     x         -
                  grid_latitude                              x         -
                  grid_longitude                             x         -
                  latitude                                   x         -
                  longitude                                  x         -
                  site_altitude                              x         -
                  site_name                                  x         -
                  site_type                                  x         -
                  surface_altitude                           x         -
                  forecast_period                            -         x
             Scalar coordinates:
                  atmosphere_hybrid_height_coordinate: 20.00... m, \
bound=(0.0, 49.99...) m
                  forecast_day: 1.0 Days
                  model_level_number: 1
                  sigma: 0.99..., bound=(1.0, 0.99...)
             Attributes:
                  Conventions: CF-1.5
                  STASH: m01s34i001
                  label: aqum_oper
                  short_name: O3
                  source: Data from Met Office Unified Model
             Cell methods:
                  mean: time (1 hour)

        Site names can be extracted (similarly other coordinate information):

        >>> print(MD.sites_cube_list[0].coord('site_name').points.astype(str))
        ... # doctest: +NORMALIZE_WHITESPACE
         ['Lullington_Heath' 'Wicken_Fen' 'Rochester' 'Sibton' 'Lough_Navar'
          'Strath_Vaich' 'Yarner_Wood' 'Eskdalemuir' 'Bush_Estate' 'Aston_Hill'
          'Glazebury' 'Ladybower' 'Harwell' 'Bottesford' 'High_Muffles']
        """

        #Set up regular lat-lon coord system
        lat_lon_coord_sys = iris.coord_systems.GeogCS(
            semi_major_axis=iris.fileformats.pp.EARTH_RADIUS)

        #Get gridded coord sys
        if not self.gridded_cube_list:
            raise ValueError("gridded_cube_list has no cubes")
        xcoord = self.gridded_cube_list[0].coord(xcoord_name)
        gridded_coord_sys = xcoord.coord_system.as_cartopy_crs()

        #Convert obs lat-lon onto model grid
        result = gridded_coord_sys.transform_points(
            lat_lon_coord_sys.as_cartopy_crs(),
            np.asarray(sites_data['longitude']),
            np.asarray(sites_data['latitude']))
        lon_modelgrid = result[:, 0]
        lat_modelgrid = result[:, 1]

        #Check if any other coordinates exist on x/y axes
        # - these will need to be moved to site_id axes
        xcoord_dim = self.gridded_cube_list[0].coord_dims(xcoord_name)[0]
        ycoord_dim = self.gridded_cube_list[0].coord_dims(ycoord_name)[0]
        xdim_coords = self.gridded_cube_list[0].coords(
            contains_dimension=xcoord_dim)
        ydim_coords = self.gridded_cube_list[0].coords(
            contains_dimension=ycoord_dim)
        # Generate list of x/y-axis coordinate names
        xydim_coords = []
        for coord in xdim_coords + ydim_coords:
            xydim_coords.append(coord[0].name())
        xydim_coords = list(set(xydim_coords)) #Remove duplicates

        #Get list of site ids, derived from site lon and lat
        site_ids = generate_siteids(sites_data['longitude'],
                                    sites_data['latitude'])

        #Can now loop over gridded cubes
        site_cubes = iris.cube.CubeList()
        for cube in self.gridded_cube_list:
            for isite, site_data in enumerate(sites_data):

                #Interpolate to point on gridded cube
                point = [(xcoord_name, lon_modelgrid[isite]),
                         (ycoord_name, lat_modelgrid[isite])]

                #Interpolate to point
                #Any points outside of gridded cube are set to nan
                sitecube = cube.interpolate(point,
                                            iris.analysis.Linear
                                            (extrapolation_mode='nan'))


                # Add site id to cube as dim_coord
                sitecube.add_aux_coord(iris.coords.DimCoord(
                    site_ids[isite], long_name='site_id'))

                #Add new axis - this is always on dimension 0
                sitecube = iris.util.new_axis(sitecube, 'site_id')

                #Remove any aux factory coordinates
                # (these can't easily be moved to site_id dimension
                #  - may want to revisit this if they are needed at somepoint)
                for aux_factory in sitecube.aux_factories:
                    sitecube.remove_aux_factory(aux_factory)

                #Move any other coordinates from x & y coord
                # to auxillary coords on same data dim.
                for coord_name in xydim_coords:
                    #Check coord hasn't been removed as aux_factory
                    try:
                        coord = sitecube.coord(coord_name)
                        sitecube.remove_coord(coord)
                        if coord_name != 'latitude' \
                           and coord_name != 'longitude':
                            #Don't add latitude or longitude back in as these
                            # are added as part of site location
                            sitecube.add_aux_coord(coord, 0)
                    except:
                        pass


                #   Add site coordinates to cube
                sitecube.add_aux_coord(iris.coords.AuxCoord(
                    np.array(site_data['site_name'],
                             dtype=sites_data['site_name'].dtype),
                    long_name='site_name'), 0)
                if 'abbrev' in site_data.dtype.names:
                    sitecube.add_aux_coord(iris.coords.AuxCoord(
                        np.array(site_data['abbrev'],
                                 dtype=sites_data['abbrev'].dtype),
                        long_name='abbrev'), 0)
                if 'site_type' in site_data.dtype.names:
                    sitecube.add_aux_coord(iris.coords.AuxCoord(
                        np.array(site_data['site_type'],
                                 dtype=sites_data['site_type'].dtype),
                        long_name='site_type'), 0)
                if 'site_altitude' in site_data.dtype.names:
                    sitecube.add_aux_coord(iris.coords.AuxCoord(
                        site_data['site_altitude'],
                        long_name='site_altitude',
                        units='m'), 0)

                #Add site lat and lon
                lon_coord = iris.coords.AuxCoord(site_data['longitude'],
                                                 standard_name='longitude',
                                                 units='degrees',
                                                 coord_system=lat_lon_coord_sys)
                sitecube.add_aux_coord(lon_coord, 0)
                lat_coord = iris.coords.AuxCoord(site_data['latitude'],
                                                 standard_name='latitude',
                                                 units='degrees',
                                                 coord_system=lat_lon_coord_sys)
                sitecube.add_aux_coord(lat_coord, 0)

                site_cubes.append(sitecube)

        self.sites_cube_list = site_cubes.concatenate()

        return self.sites_cube_list

    def save_ts(self, filename='sites_cube_list.nc', sites_cube_list=None):
        """
        Saves sites_cube_list to netcdf file
        If sites_cube_list passed in, then this is saved instead.
        """
        print('Saving to', filename)
        if sites_cube_list:
            iris.save(sites_cube_list, filename)
        else:
            iris.save(self.sites_cube_list, filename)


    def load_ts(self, filename='sites_cube_list.nc'):
        """
        Restore sites_cube_list from netcdf file
        Returns sites_cube_list directly, but also sets
        sites_cube_list in ADAQdata object
        Cubes are returned in alphabetical order according
        to cube name.
        """
        #print('Restoring from', filename)
        self.sites_cube_list = iris.load(filename)
        #Sort cubes alphabetically to simplify tests
        #(Note iris load does not guarantee to always load
        # cubes in the same order)
        self.sites_cube_list = iris.cube.CubeList(sorted(list(
            self.sites_cube_list), key=lambda cube: cube.name()))

        return self.sites_cube_list


    def save_gridded(self, filename='gridded_cube_list.nc',
                     gridded_cube_list=None):
        """
        Saves gridded_cube_list to netcdf file
        If gridded_cube_list is passed in then this is saved instead
        """
        print('Saving to', filename)
        if gridded_cube_list is None:
            gridded_cube_list = self.gridded_cube_list
        iris.save(gridded_cube_list, filename)

    def __gridded_callback(self, cube, field, filename):
        """
        Callback function to perform a check on the gridded cube.
        Checks the gridded cube to determine if it has all
        required attributes and coordinates.
        Cubes are returned in alphabetical order according
        to cube name.        """
        try:
            self.check_gridded_cube(cube, verbose=False)
        except:
            raise iris.exceptions.IgnoreCubeException

    def load_gridded(self, filename='gridded_cube_list.nc'):
        """
        Restores gridded_cube_list
        Returns gridded_cube_list directly,
        but also sets gridded_cube_list in ADAQdata object
        """
        #print('Restoring from', filename)
        self.gridded_cube_list = iris.load(filename, callback=self.__gridded_callback)
        #Sort cubes alphabetically to simplify tests
        #(Note iris load does not guarantee to always load
        # cubes in the same order)
        self.gridded_cube_list = iris.cube.CubeList(sorted(list(
            self.gridded_cube_list), key=lambda cube: cube.name()))

        return self.gridded_cube_list

    def save_trajectory(self, filename='trajectory_cube_list.nc',
                        trajectory_cube_list=None):
        """
        Saves trajectory_cube_list to netcdf file
        If trajectory_cube_list is passed in then this is saved instead
        The trajectory NetCDF file is not pretty but it loads back OK
        using Iris.
        """
        print('Saving to', filename)
        if trajectory_cube_list is None:
            trajectory_cube_list = self.trajectory_cube_list
        iris.save(trajectory_cube_list, filename)

    def load_trajectory(self, filename='trajectory_cube_list.nc'):
        """
        Restores trajectory_cube_list
        Returns trajectory_cube_list directly,
        but also sets trajectory_cube_list in ADAQdata object
        Cubes are returned in order according
        to its PP Index coordinate.
        """
        #print('Restoring from', filename)
        self.trajectory_cube_list = iris.load(filename)
        #Sort cubes alphabetically to simplify tests
        #(Note iris load does not guarantee to always load
        # cubes in the same order)
        self.trajectory_cube_list = iris.cube.CubeList(sorted(
            list(self.trajectory_cube_list),
            key=lambda cube: cube.coord('PP Index').points[0]))

        return self.trajectory_cube_list

    def save_all(self, path='./', prefix='', suffix=''):
        """
        Save each of sites_cube_list, gridded_cube_list,
        and trajectory_cube_list, if present.
        The cube lists will be saved as "sites_cube_list.nc", etc,
        with an optional prefix and suffix.

        :param path: directory to save all cubelists to
        :param prefix: string to prepend to the default file names
        :param suffix: string to append to each name, before the file extension

        Tested as part of :any:`aq_plot`.
        """
        if self.sites_cube_list:
            self.save_ts(
                os.path.join(path, prefix+'sites_cube_list'+suffix+'.nc'))
        if self.gridded_cube_list:
            self.save_gridded(
                os.path.join(path, prefix+'gridded_cube_list'+suffix+'.nc'))
        if self.trajectory_cube_list:
            self.save_trajectory(
                os.path.join(path, prefix+'trajectory_cube_list'+suffix+'.nc'))


    def derive(self, *specifications):
        """
        Method to compute new cubes from existing ones,
        and immediately add them to this :any:`ADAQData` object.
        Both gridded cubes and site cubes will be modified, if non-empty.

        :param \\*specifications: all positional arguments should be
            dictionaries describing how to derive a particular cube,
            such as those returned by :any:`get_derived_sn_specs`.
            They may contain any of the following keys:

            * **components** - a list of short names that it depends on.
              If not specified, this quantity is expected to be directly
              available to load with no further processing.
            * **function** - a function that returns a cube of the required
              quantity. It should accept either a cubelist or a cube as its
              first argument, depending on whether there is more than one
              component or not.
            * **kwargs** - any other arguments to pass to the function.
            * **short_name** - expected short name of the new cube.
            * **standard_name** - expected standard name of the new cube.

            For convenience, if a string is encountered instead of a
            dictionary, it will automatically be passed to
            :any:`get_derived_sn_specs`.

        >>> import adaq_functions
        >>> ini_dict, sites_data, od, md_list = adaq_functions.get_exampledata(
        ...   sites_cube_list=True, gridded_cube_list=True)
        ...   # doctest: +ELLIPSIS
        Reading inifile .../example_data_1days.ini
        Number of sites:  5
        >>> md = md_list[0]
        >>> print(md) # doctest: +ELLIPSIS
        <class '....ADAQData'> - Subclass of ADAQData - Contains:
        sites_cube_list:
        0: mass_concentration_of_ozone_in_air / (ug/m3) \
(site_id: 5; time: 25)
        1: mass_concentration_of_pm2p5_dry_aerosol_in_air / (ug/m3) \
(site_id: 5; time: 25)
        gridded_cube_list:
        0: mass_concentration_of_ozone_in_air / (ug/m3) \
(time: 25; grid_latitude: 182; grid_longitude: 146)
        1: mass_concentration_of_pm2p5_dry_aerosol_in_air / (ug/m3) \
(time: 25; grid_latitude: 182; grid_longitude: 146)
        trajectory_cube_list:
        < No cubes >

        Now derive some new cubes.

        >>> md.derive("PM2p5_daily_max", "mean_O3")
        >>> print(md) # doctest: +ELLIPSIS
        <class '....ADAQData'> - Subclass of ADAQData - Contains:
        sites_cube_list:
        ...
        2: mass_concentration_of_pm2p5_dry_aerosol_in_air / (ug/m3) \
(site_id: 5; time: 2)
        3: mass_concentration_of_ozone_in_air / (ug/m3) \
(site_id: 5)
        gridded_cube_list:
        ...
        2: mass_concentration_of_pm2p5_dry_aerosol_in_air / (ug/m3) \
(time: 2; grid_latitude: 182; grid_longitude: 146)
        3: mass_concentration_of_ozone_in_air / (ug/m3) \
(grid_latitude: 182; grid_longitude: 146)
        trajectory_cube_list:
        < No cubes >
        """

        for spec in specifications:
            if isinstance(spec, str):
                spec = self.get_derived_sn_specs([spec])[0]

            if "components" not in spec:
                continue

            constraints = [iris.AttributeConstraint(short_name=sn)
                           for sn in spec["components"]]
            func = spec.get("function")
            kwargs = spec.get("kwargs", {})

            #loop over each relevant cubelist
            for cubelist in (self.sites_cube_list,
                             self.gridded_cube_list):
                if not cubelist:
                    continue

                #extract a cubelist containing the required cubes
                #note that if there is one constraint, then iris will
                # return a single cube (instead of a cubelist containing
                # a single cube).
                try:
                    components = cubelist.extract_strict(constraints)
                except iris.exceptions.ConstraintMismatchError:
                    warnings.warn("could not find all components {} required "
                                  "for derivation".format(spec['components']))
                    continue

                #call the derivation function
                #if not specified, sum the cubes instead
                if func:
                    output = func(components, **kwargs)
                elif len(constraints) > 1:
                    output = components[0].copy()
                    for cube in components[1:]:
                        attrs = dict(output.attributes)
                        output += cube

                        #iris removes attributes after adding, so put back
                        # the ones that both cubes agreed on
                        for key, val in cube.attributes.items():
                            if key in attrs and attrs[key] == val:
                                output.attributes[key] = val
                else:
                    continue

                if (output is components
                        or (output is None
                            and 1 < len(constraints) < len(components))):
                    #these conditions both suggest that function worked
                    # by appending to the input list
                    #we should therefore pass any extra cubes on to the
                    # main list
                    cubelist.extend(components[len(constraints):])

                elif isinstance(output, (list, tuple)):
                    #the output is a new list - this includes cubelists
                    #we assume that they are all new cubes, but it
                    # may be worth checking that this isn't just a copy
                    # of the input
                    cubelist.extend(output)

                elif isinstance(output, iris.cube.Cube):
                    #set attributes that the function may not have done
                    if "standard_name" in spec:
                        output.rename(spec["standard_name"])
                    if "short_name" in spec:
                        output.attributes["short_name"] = spec["short_name"]

                    cubelist.append(output)

    def get_derived_sn_specs(self, short_name_list):
        """
        Wrapper method to call :func:`get_derived_sn_specs` with
        :data:`ADAQData.DERIVED_SHORT_NAME_DICT`.
        This attribute should be overwritten by subclasses, to define
        which derivations can be handled.
        """
        return get_derived_sn_specs(
            short_name_list,
            derived_short_name_dict=self.DERIVED_SHORT_NAME_DICT)



if __name__ == '__main__':

    import doctest
    doctest.testmod()
