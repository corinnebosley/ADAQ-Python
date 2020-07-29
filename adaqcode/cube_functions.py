"""
Functions which work on cubes (and don't belong in any other cube... modules)
"""

import numpy as np
import iris
from iris.analysis.cartography import rotate_pole, unrotate_pole
import iris.analysis.trajectory as trajectory


def get_latlon_cube(lonmin, lonmax, latmin, latmax, dlon, dlat,
                    exact_maxlatlon=False):
    '''
    Get a cube with a regular lat-lon grid.

    :param lonmin: float, mininum longitude to be included in cube
    :param lonmax: float, maximum longitude to be included in cube
    :param latmin: float, mininum latitude to be included in cube
    :param latmax: float, maximum latitude to be included in cube
    :param dlon: float, required distance between longitude values.
    :param dlat: float, required distance between latitude values.
                 Note if this does not give lonmax as the maximum
                 value, then dlon is adjusted to ensure latmax is
                 the final value.
    :param exact_maxlatlon: logical. Controls behaviour in the case where
                            lonmax != lonmin + N * dlon or
                            latmax != latmin + N * dlat and N is an integer.
                            If True, then ensures that lonmax and lonmin
                            are definately included in the output cube and
                            dlon/dlat are adjusted (by increasing in size).
                            If False, then ensures that dlon and dlat
                            are fixed and lonmax and latmax are adjusted
                            (increased to ensure requested lonmax and latmax
                            are still included in cube).

    :returns: cube, with regular lat-lon coordinates and data all set to 1.

    Example usage:

    >>> llcube = get_latlon_cube(lonmin=1., lonmax=5, latmin=1., latmax=4.,
    ...     dlon=1., dlat=1.)
    >>> print(llcube)
    Regular lat-lon grid / (unknown)    (longitude: 6; latitude: 5)
         Dimension coordinates:
              longitude                           x            -
              latitude                            -            x
    >>> print(llcube.coord('longitude'))
    DimCoord(array([1., 2., 3., 4., 5., 6.]), standard_name='longitude', \
units=Unit('degrees'), coord_system=GeogCS(6371229.0))
    >>> print(llcube.coord('latitude'))
    DimCoord(array([1., 2., 3., 4., 5.]), standard_name='latitude', \
units=Unit('degrees'), coord_system=GeogCS(6371229.0))

    Example with lonmax(=5.5) not matching lonmin+N*dlon:

    >>> llcube = get_latlon_cube(lonmin=1., lonmax=5.5, latmin=1., latmax=4.,
    ...    dlon=1., dlat=1., exact_maxlatlon=False)
    >>> print(llcube.coord('longitude').points)
    [1. 2. 3. 4. 5. 6.]
    >>> llcube = get_latlon_cube(lonmin=1., lonmax=5.5, latmin=1., latmax=4.,
    ...    dlon=1., dlat=1., exact_maxlatlon=True)
    >>> print(llcube.coord('longitude').points)
    [1.    2.125 3.25  4.375 5.5  ]
    '''

    if lonmin < 0 or lonmax < 0:
        lonmin += 360.
        lonmax += 360.

    lat_lon_coord_system = iris.coord_systems.GeogCS(
        semi_major_axis=iris.fileformats.pp.EARTH_RADIUS)

    #Calculate longitude and latitude points to be used
    if exact_maxlatlon:
        #Set nlon such that when linspace is used,
        #requested lonmax is the maximum longitude
        nlon = 1 + int(np.round((lonmax - lonmin) / dlon))
    else:
        #Adjust lonmax to ensure dlon is kept as requested
        nlon = 2 + int(np.round((lonmax - lonmin) / dlon))
        lonmax = lonmin + (nlon-1)*dlon
    lons = np.linspace(lonmin, lonmax, num=nlon)

    if exact_maxlatlon:
        #Set nlon such that when linspace is used,
        #requested lonmax is the maximum longitude
        nlat = 1 + int(np.round((latmax - latmin) / dlat))
    else:
        #Adjust latmax to ensure dlat is kept as requested
        nlat = 2 + int(np.round((latmax - latmin) / dlat))
        latmax = latmin + (nlat-1)*dlat
    lats = np.linspace(latmin, latmax, num=nlat)

    #Set up coordinates
    lon_coord = iris.coords.DimCoord(lons,
                                     standard_name='longitude',
                                     units='degrees',
                                     coord_system=lat_lon_coord_system)
    lat_coord = iris.coords.DimCoord(lats,
                                     standard_name='latitude',
                                     units='degrees',
                                     coord_system=lat_lon_coord_system)

    #Put into a cube
    data = np.ones((len(lons), len(lats)))
    cube = iris.cube.Cube(data,
                          long_name='Regular lat-lon grid',
                          dim_coords_and_dims=[(lon_coord, 0),
                                               (lat_coord, 1)])
    return cube

def guess_coord_names(cube, axes):
    """
    Guess the name of the coordinate corresponding to the required axes

    :param cube: iris Cube
    :param axes: List of axes, eg 'X','Y','Z','T'

    :returns: List of coordinate names corresponding to these axes.
              If an axes not found, then value in list is None.
              Will try to return dimension coordinates if possible.

    >>> import config
    >>> sample_datadir = config.SAMPLE_DATADIR+'gridded_cube_list/'
    >>> cube = iris.load_cube(sample_datadir+'aqum_oper_1days.nc',
    ... 'surface_altitude')
    >>> print(cube.summary(True)) # doctest: +NORMALIZE_WHITESPACE
    surface_altitude / (m)  (grid_latitude: 182; grid_longitude: 146)
    >>> coord_names = guess_coord_names(cube,['X','Y','Z'])
    >>> coord_names == ['grid_longitude', 'grid_latitude', None]
    True
    """

    coord_names = [None]*len(axes)
    for coord in cube.coords():
        axis = iris.util.guess_coord_axis(coord)
        for i, ax in enumerate(axes):
            if axis == ax:
                if coord_names[i] is None:
                    coord_names[i] = coord.name()

    return coord_names


def extract_section(cube, waypoints, n_sample_points=30):
    """
    Extract a section from the cube given lat-lon waypoints. This works
    on any cube with x and y dimension coordinates. If the cube also contains
    Z or T etc dimension coordinates then these are retained in the extracted
    section.

    :param cube: Iris cube
    :param waypoints: List of dictionaries, whose keys are 'latitude' and
                      'longitude' and whose values are the lat/lon points
                      that should be included in the section.
    :param n_sample_points: Number of sample points to include in section

    :returns: An iris cube which no longer has X and Y as dimension coordinates
              but instead has an 'i_sample_point' coordinate.
              If latitude and longitude were not in the original cube, they
              are added to to the new section cube.

    Load sample data - a gridded surface-only, time-varying cube:

    >>> import config
    >>> sample_datadir = config.SAMPLE_DATADIR+'gridded_cube_list/'
    >>> cube = iris.load_cube(sample_datadir+'aqum_oper_1days.nc',
    ... 'mass_concentration_of_ozone_in_air')
    >>> print(cube) # doctest: +NORMALIZE_WHITESPACE
    mass_concentration_of_ozone_in_air / (ug/m3) (time: 25; grid_latitude: 182; \
grid_longitude: 146)
         Dimension coordinates:
              time                 x                  -                    -
              grid_latitude        -                  x                    -
              grid_longitude       -                  -                    x
         Auxiliary coordinates:
              forecast_period      x                  -                    -
              surface_altitude     -                  x                    x
         Derived coordinates:
              altitude             -                  x                    x
         Scalar coordinates:
              atmosphere_hybrid_height_coordinate: 20.000338 m, \
bound=(0.0, 49.998882) m
              forecast_day: 1.0 Days
              model_level_number: 1
              sigma: 0.9977165, bound=(1.0, 0.99429625)
         Attributes:
              Conventions: CF-1.5
              STASH: m01s34i001
              label: aqum_oper
              short_name: O3
              source: Data from Met Office Unified Model
         Cell methods:
              mean: time (1 hour)

    Set up waypoints list:

    >>> waypoints = [
    ... {'latitude': 50.8, 'longitude': -1.8},
    ... {'latitude': 51.2, 'longitude':	-1.2},
    ... {'latitude': 51.4, 'longitude':	-0.9}]

    Extract section along these waypoints:

    >>> section = extract_section(cube, waypoints)
    >>> print(section)
    mass_concentration_of_ozone_in_air / (ug/m3) (time: 25; i_sample_point: 30)
         Dimension coordinates:
              time                                    x                   -
              i_sample_point                          -                   x
         Auxiliary coordinates:
              forecast_period                         x                   -
              grid_latitude                           -                   x
              grid_longitude                          -                   x
              latitude                                -                   x
              longitude                               -                   x
              surface_altitude                        -                   x
         Derived coordinates:
              altitude                                -                   x
         Scalar coordinates:
              atmosphere_hybrid_height_coordinate: 20.000338 m, \
bound=(0.0, 49.998882) m
              forecast_day: 1.0 Days
              model_level_number: 1
              sigma: 0.9977165, bound=(1.0, 0.99429625)
         Attributes:
              Conventions: CF-1.5
              STASH: m01s34i001
              label: aqum_oper
              short_name: O3
              source: Data from Met Office Unified Model
         Cell methods:
              mean: time (1 hour)

    """

    xcoord, ycoord = guess_coord_names(cube, ['X', 'Y'])

    if None in [xcoord, ycoord]:
        raise ValueError('Can not calculate section for this cube')

    lat_lon_coord_system = iris.coord_systems.GeogCS(
        semi_major_axis=iris.fileformats.pp.EARTH_RADIUS)
    lons = np.array([waypoint['longitude'] for waypoint in waypoints])
    lats = np.array([waypoint['latitude'] for waypoint in waypoints])

    #Calculate waypoints in cube coordinates
    if xcoord == 'grid_longitude' and ycoord == 'grid_latitude':
        #Rotated pole
        pole_lon = cube.coord(xcoord).coord_system.grid_north_pole_longitude
        pole_lat = cube.coord(ycoord).coord_system.grid_north_pole_latitude

        #Perform rotation
        rot_lons, rot_lats = rotate_pole(lons, lats, pole_lon, pole_lat)

        #Put back into waypoints format
        cube_waypoints = [{xcoord: lon, ycoord: lat}
                          for lat, lon in zip(rot_lats, rot_lons)]
    elif xcoord == 'projection_x_coordinate' and \
         ycoord == 'projection_y_coordinate':
        #Other coordinate system (note this may work for x/ycoords other than
        #those considered here
        ll_crs = lat_lon_coord_system.as_cartopy_crs()
        cube_crs = cube.coord(xcoord).coord_system.as_cartopy_crs()
        #Convert to lat/lon points
        cube_lonlats = ll_crs.transform_points(cube_crs, lons, lats)
        cube_lons = cube_lonlats[:, 0]
        cube_lats = cube_lonlats[:, 1]
        #Put back into waypoints format
        cube_waypoints = [{xcoord: lon, ycoord: lat}
                          for lat, lon in zip(cube_lats, cube_lons)]

    elif xcoord == 'longitude' and ycoord == 'latitude':
        cube_waypoints = waypoints
    else:
        raise ValueError('Unable to convert cube x/y points to lat/lon')

    #Specify the trajectory (straight-line) by giving the
    # start and end of the line and specifying how many points there
    # should be on the line (in this case 30)
    sample_points = trajectory.Trajectory(cube_waypoints,
                                          sample_count=n_sample_points)

    #Extract the data from the cube at the sample points
    #(Note this only works with iris vn2.2+)
    section = sample_points.interpolate(cube)

    #Rename new index coordinate
    section.coord('index').rename('i_sample_point')

    #Also add latitude and longitude coords to section
    #if not in original cube
    if xcoord != 'longitude' and ycoord != 'latitude':

        cube_lons_samplepts = np.array(
            [d[xcoord] for d in sample_points.sampled_points])
        cube_lats_samplepts = np.array(
            [d[ycoord] for d in sample_points.sampled_points])


        if xcoord == 'grid_longitude' and ycoord == 'grid_latitude':
            section_lons, section_lats = unrotate_pole(cube_lons_samplepts,
                                                       cube_lats_samplepts,
                                                       pole_lon, pole_lat)
        elif xcoord == 'projection_x_coordinate' and \
             ycoord == 'projection_y_coordinate':
            #Other coordinate system (note this may work for x/ycoords other than
            #those considered here
            lonlats = cube_crs.transform_points(ll_crs,
                                                cube_lons_samplepts,
                                                cube_lats_samplepts)
            section_lons = lonlats[:, 0]
            section_lats = lonlats[:, 1]


        #Set up lat/lon coords and add as aux coords
        lon_coord = iris.coords.AuxCoord(section_lons,
                                         standard_name='longitude',
                                         units='degrees',
                                         coord_system=lat_lon_coord_system)

        lat_coord = iris.coords.AuxCoord(section_lats,
                                         standard_name='latitude',
                                         units='degrees',
                                         coord_system=lat_lon_coord_system)
        section.add_aux_coord(lon_coord, section.coord_dims(xcoord))
        section.add_aux_coord(lat_coord, section.coord_dims(ycoord))



    return section




if __name__ == '__main__':

    import doctest
    doctest.testmod()
