"""
Code required for checking that the domain defined by the UM namelist items
completely contains the lat-lon box specified. This therefore includes
code for creating a cube from given UM grid specifications, plus
code for plotting regular lat-lon boxes.
"""
from __future__ import division
from __future__ import print_function
from six.moves.builtins import str

import numpy as np
import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

import os
import sys
#Find adaq directory -
# should be in directory above the one containing this file
adaq_path = os.path.dirname(os.path.realpath(__file__))+'/../'
sys.path.append(adaq_path)
from adaqcode.cube_functions import get_latlon_cube

def make_a_cube(data, name, units, latarr, lonarr, polelata, polelona):
    '''
    Make an iris cube, given lat and lon points, pole coordinates,
    data array, cube name and units.
    '''

    # set up longitude and latitude coords usng the arrays
    # passed in
    lon_coord = iris.coords.DimCoord(lonarr,
        standard_name='grid_longitude',
        units='degrees',
        coord_system=
            iris.coord_systems.RotatedGeogCS(polelata, polelona))
    lat_coord = iris.coords.DimCoord(latarr,
        standard_name='grid_latitude',
        units='degrees',
        coord_system=
            iris.coord_systems.RotatedGeogCS(polelata, polelona)
        )

    cube = iris.cube.Cube(data, long_name = name,
                                units = units,
                                dim_coords_and_dims=[
                                    (lon_coord, 0),
                                    (lat_coord, 1)])
    
    return cube
    

def get_cubes(delta_lat=0.110000, delta_lon=0.110000,
            frstlata=-8.4100, frstlona=352.9500,
            polelata=37.5000, polelona=177.5000,
            global_row_length=146, global_rows=182,
            rimweightsa=[1.0,1.0,1.0,1.0,1.0,0.750,0.500,0.250,0.0],
            l_endgame=True):
    '''
    Get a theta, u and v iris cubes as defined by the UM namelist items.
    Do not include rim or blending zone.
    '''

    # first calculate rim/blending zone region size
    # based on number of values greater than zero
    rimsize = len([w for w in rimweightsa if w>0.0])

    # number of valid x and y points are reduced by 2 * rim size
    nx = global_row_length - 2*rimsize
    ny = global_rows - 2*rimsize

    # get revised start lon and lat
    # increase start lat/lon by rim width
    startlon = frstlona + rimsize*delta_lon
    stoplon = startlon + (nx-1)*delta_lon
    startlat = frstlata + rimsize*delta_lat
    stoplat = startlat + (ny-1)*delta_lat

    # use these to make theta lons and lats
    theta_lon = np.linspace(startlon, stoplon, num=nx)
    theta_lat = np.linspace(startlat, stoplat, num=ny)
    
    # make dummy data to allow plotting.
    # Use random data so individual grid boxes can be seen easier.
    theta_data = np.random.random((len(theta_lon), len(theta_lat)))

    # now make theta cube
    theta_cube = make_a_cube(theta_data, 'theta_grid',
                             'K', theta_lat, theta_lon,
                             polelata, polelona)

    if l_endgame:
        # u longitudes same number as theta 
        # but 1/2 gridbox to west
        u_lon = theta_lon.copy() - delta_lon/2.0

        # u latuitudes same number and value as theta
        u_lat = theta_lat

        # data - can reuse theta random data as same no points
        u_data = theta_data

        u_cube = make_a_cube(u_data, 'u_grid',
                             'm s-1', u_lat, u_lon,
                             polelata, polelona)

        # v longitudes same as theta
        v_lon = theta_lon

        # v latitudes one extra with first lat 1/2 grid 
        # box to south
        v_lat = np.resize(theta_lat, (len(theta_lat)+1)) - delta_lat/2.0
        # need to set the last latitude from last but one
        v_lat[-1] = v_lat[-2] + delta_lat

        v_data = np.random.random((len(v_lon), len(v_lat)))
       
        v_cube = make_a_cube(v_data, 'v_grid',
                             'm s-1', v_lat, v_lon,
                             polelata, polelona)

        # b grid longitudes are the same as the u longitudes
        b_lon = u_lon
        # b grid latitudes are the same as the v latitudes
        b_lat = v_lat
        b_data = np.random.random((len(b_lon), len(b_lat)))
        b_cube = make_a_cube(b_data, 'b_grid',
                             'm s-1', b_lat, b_lon,
                             polelata, polelona)
       

    else:
        raise ValueError("New dynamics not yet supported")

    return theta_cube, u_cube, v_cube, b_cube

def plot_cube(cube):
    '''
    Plot cube with random data and grid points marked on a map
    '''
    qplt.points(cube, s=0.2, c='k')
    qplt.pcolormesh(cube, cmap='rainbow', alpha=0.2)
    # Add coastlines to the map created by pcolormesh.
    plt.gca().coastlines()

def plot_latlon_box(box_lonmin=-18.0, box_lonmax=16.9,
                    box_latmin=43.3, box_latmax=63.2):
    '''
    Overplot a lat-lon box on a map
    '''

    plt.plot([box_lonmin, box_lonmin], [box_latmin, box_latmax],
             linewidth=2, transform=ccrs.PlateCarree(), color='r')
    plt.plot([box_lonmax, box_lonmax], [box_latmin, box_latmax],
             linewidth=2, transform=ccrs.PlateCarree(), color='r')
    plt.plot([box_lonmin, box_lonmax], [box_latmin, box_latmin],
             linewidth=2, transform=ccrs.PlateCarree(), color='r')
    plt.plot([box_lonmin, box_lonmax], [box_latmax, box_latmax],
             linewidth=2, transform=ccrs.PlateCarree(), color='r')


def get_latlon_cube_fromfile(
                        filename='/data/nwp1/apdl/HRES_ENS_2015100600+090.nc'):
    '''
    Load a cube from file instead for the lat-lon box domain
    '''

    macc_ens = iris.load('/data/nwp1/apdl/HRES_ENS_2015100600+090.nc',
                         iris.Constraint(level=0))[0]
    box_cube = macc_ens[0, :, :] #Remove time dimension
    box_cube.transpose() #Get first dim as lon, second dim as lat.
    lat_lon_coord_system = iris.coord_systems.GeogCS(semi_major_axis =
                         iris.fileformats.pp.EARTH_RADIUS)
    box_cube.coord('latitude').coord_system = lat_lon_coord_system
    box_cube.coord('longitude').coord_system = lat_lon_coord_system
    box_cube.data[:] = 1.    

    return box_cube

def interp_cube_box(model_cube, box_cube):
    '''
    Interpolate the model cube onto the required lat-lon box cube.
    Count the number of nan points - if zero, then can set success=True.
    If any nan points, these are instead set to -2. for clarity in plotting.
    '''
    
    interp_cube = model_cube.regrid(box_cube,
                                iris.analysis.Linear(extrapolation_mode='nan'))
    nfailed = np.count_nonzero(np.isnan(interp_cube.data))
    
    print('Number of failed points: ', nfailed, \
          'from max of ', interp_cube.data.size, 'points')
    
    i = np.where(np.isnan(interp_cube.data))
    interp_cube.data[i] = -2.

    return interp_cube, nfailed  

def plot_interp(domain_cube, box_lonmin, box_lonmax,
                    box_latmin, box_latmax, box_cube):
    '''Plot a cube of data on a proposed new domain,
    test if it can be interpolated
    onto the required lat-lon box and plot the results
    of the interpolation
    '''

    plt.figure(figsize=(12, 8))
    plt.subplot(121) 
    #Plot this cube, along with required lat-lon domain
    plot_cube(domain_cube)
    plot_latlon_box(box_lonmin=box_lonmin, box_lonmax=box_lonmax,
                    box_latmin=box_latmin, box_latmax=box_latmax)
    plt.title('Model cube with required box', fontsize=10)
    #plt.show()


    #Now try interpolating the model (theta) cube onto the required domain
    interp_cube, nfailed = interp_cube_box(domain_cube, box_cube)

    #Plot results
    plt.subplot(122) 
    plot_cube(interp_cube)
    plot_latlon_box(box_lonmin=box_lonmin, box_lonmax=box_lonmax,
                    box_latmin=box_latmin, box_latmax=box_latmax)
    plt.title('Interpolated model grid on required grid '\
              +'\n Failed points='+str(nfailed),
              fontsize=10)
    plt.suptitle(domain_cube.name(), fontsize=15)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.05)
    plt.show()

    return nfailed
    

def check_domain(delta_lat=0.110000, delta_lon=0.110000,
                 frstlata=-8.4100, frstlona=352.9500,
                 polelata=37.5000, polelona=177.5000,
                 global_row_length=146, global_rows=182,
                 rimweightsa=[1.0,1.0,1.0,1.0,1.0,0.750,0.500,0.250,0.0],
                 l_endgame=True,
                 box_lonmin=-18.0, box_lonmax=16.9,
                 box_latmin=43.3, box_latmax=63.2,
                 box_dlon=0.11, box_dlat=0.11):
    '''
    Check that the domain defined by the UM namelist items
    completely contains the lat-lon box specified.
    Do not include rim or blending zone in valid
    boxes
    '''

    # get cubes to be tested
    model_cubes = get_cubes(delta_lat=delta_lat,
                          delta_lon=delta_lon,
                          frstlata=frstlata, frstlona=frstlona,
                          polelata=polelata, polelona=polelona,
                          global_row_length=global_row_length,
                          global_rows=global_rows,
                          rimweightsa=rimweightsa,
                          l_endgame=l_endgame)

    #Get required lat-lon domain in a cube
    box_cube = get_latlon_cube(lonmin=box_lonmin, lonmax=box_lonmax,
                               latmin=box_latmin, latmax=box_latmax,
                               dlon=box_dlon, dlat=box_dlat,
                               exact_maxlatlon=True)
    #box_cube = get_latlon_cube_fromfile()

    # loop over cubes to be tested
    total_failed = 0
    for cube in model_cubes:
        nfailed = plot_interp(cube, box_lonmin, box_lonmax,
                              box_latmin, box_latmax, box_cube)
        total_failed += nfailed
        #break

    if total_failed == 0:
        success = True
        print('Model cubes all successfully within required domain')
    else:
        success = False
        print('Model cubes not all within required domain')
    
    return model_cubes, box_cube, success

if __name__ == '__main__' :

    #model_cubes, box_cube, success = check_domain()

    #CAMS - minimum enclosing grid + DA considerations
    # DA - delta lat = delta lon
    #    - n gridpts in both directions must be a multiple of 2, 3 and/or 5 ONLY
    #    - eg 576 = 3*3*2*2*2*2*2*2 = 3^2 * 2^6
    #    -    432 = 2^4 * 3^3
    model_cubes, box_cube, success = check_domain(#Required lat-lon box:
                                                  box_lonmin=-25.0,
                                                  box_lonmax=45.0,
                                                  box_latmin=30.0,
                                                  box_latmax=70.0,
                                                  box_dlon=0.1,
                                                  box_dlat=0.1,
                                                  #Model domain:
                                                  delta_lat=0.11,
                                                  delta_lon=0.11,
                                                  frstlata=-12.43,
                                                  frstlona=328.44,
                                                  global_row_length=576,
                                                  global_rows=432,
                                                  polelata=50.0,
                                                  polelona=190.0)
