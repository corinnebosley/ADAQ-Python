#!/usr/bin/env python

'''
Test script to plot the extent of overlap between multiple site back run plumes.

Designed to investigate where/by how much back run plumes overlap if, for example, a contaminant is detected at several locations.
If for each location a back run is performed, the extent of overlap between the plumes could indicate the most likely source location.
This example handles 10-day back run output. Output files for each day (24 hr period) for each site, over different periods.

Calculate one plume per site (combine time-separated fields files over total run duration) and give a uniform value of 1, e.g. plume exists here,
then sum all site plumes. An area where two plumes overlap will have a value of 2, and where three overlap will equal 3 etc.

Distribution / density of sites should be taken into account. 
The greater the density of sites, the more plume overlap - subjective selection may be required to obtain an even distribution.

'''
from __future__ import print_function

from six.moves.builtins import str
from six.moves.builtins import range
import os.path
import iris
import iris.coords
import numpy as np
import iris.plot as iplt
import iris.quickplot as qplt
import datetime as dt
from datetime import date
import matplotlib
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import matplotlib.colors as mcols
import glob

def get_source_location(run):
    '''Each back run should have only one source. '''
    sl_file = run + '/sourcelocations.txt'
    with open(sl_file) as f:
        line1 = f.readline()
        line2 = f.readline()
        line3 = f.readline()
        row = line3.strip().split(',') 
        sl_x = float(row[2])
        sl_y = float(row[3])

    return sl_x, sl_y

def set_latlon_constr(lat_min, lat_max, lon_min, lon_max):
    ''' Set constraints for lat and lon - e.g. for Europe. Make sure all source locations are contained within this '''
    constraint = iris.Constraint(coord_values={'latitude':lambda cell:lat_min<cell<lat_max, 
                                               'longitude':lambda cell:lon_min<cell<lon_max})
    return constraint



if __name__ == '__main__':

    lat_min=30
    lat_max=85
    lon_min=-20
    lon_max=55

    workdir = '/data/users/apdg/python_sample_data/name_plume_overlap/'

    rundirs = sorted(glob.glob(workdir+'BACK_RUN_*'))

    if len(rundirs) == 0:
            raise AssertionError('No back run directories found')

    constraint = set_latlon_constr(lat_min, lat_max, lon_min, lon_max)

    source_locations_x = []
    source_locations_y = []
    
    # set counter for number of sites
    N = 0
    
    # choose height level index if applicable - if not, leave as 0
    height_idx = 0 

    # for each site (i.e. each back run) plot total plume extent (over all time steps)
    # then sum all the site cubes to make one cube - observe overlap

    for run in rundirs:
    
        # get source location info and add to list to plot marker
        sl_x, sl_y = get_source_location(run)
        source_locations_x.append(sl_x)
        source_locations_y.append(sl_y)
            
        # read in all NAME output fields files for that run and create cube list
        # should only be a fields grid 1 , need to specify which fields grid required
        output_fields = glob.glob(run+'/Fields_grid1*.txt')
			
        if len(output_fields) == 0:
            raise AssertionError('No Fields_grid1 files found. If you wish to use Fields_grid2 files for example, please modify code to the correct Fields_grid file. Must specify only one type of Fields_grid file.')
        
        # load a cube of all the fields files for that run (e.g. for all height levels and timesteps)
        cube = iris.load_cube(output_fields)
        
        # apply lat-lon constraint to plot only the region we are interested in
        cube = iris.Constraint.extract(constraint,cube)
        #print cube
        for coord in cube.coords():
            if coord.name() == 'time':
                ntimes = coord.shape[0]
            elif coord.name() == 'latitude':
                continue
            elif coord.name() == 'longitude':
                continue
            elif coord.name() == 'height':
                cube = cube[height_idx] # can only have one height level - also assuming height dim comes before time, lat, lon dims!
            else:
                cube = cube[0] # assuming any other dim will come ahead of the others
                print('WARNING: too many dimensions. A dimension has been removed to proceed with plotting. Results may be incorrect.')
                
        print(cube)
        # sum through each timestep to get one plume over entire back run duration - we will overwrite the values after
        for t in range(ntimes):
            print(t)
            cube_t = cube[t]
            ndims = cube_t.ndim
            if ndims > 2:
                raise AssertionError('Too many cube dimensions to plot!')
            if t == 0:
                cube_tsum = cube_t
            else:
                cube_tsum = iris.analysis.maths.add(cube_tsum,cube_t, dim=None, ignore=True, in_place=True)  # sum next cube in list and current cube

        # create a copy of the time-summed cube before changing the values
        cube_copy = cube_tsum.copy()
        
        # overwrite values with 1 or 0 depending on whether the value at that lat-lon point is greater than zero (i.e. if plume exists)
        for i in range(cube_tsum.coord('latitude').shape[0]):
            for j in range(cube_tsum.coord('longitude').shape[0]):
                if cube_tsum.data[i][j] > 0.0:
                    cube_copy.data[i][j] = 1.0
                else:
                    cube_copy.data[i][j] = 0.0  
         
        N+=1

        # now sum the binary value site cubes to create 'hotspot' cube to show extent of overlap between site plumes
        if N == 1:    
            sites_cube = cube_copy
        else:
            sites_cube = iris.analysis.maths.add(sites_cube,cube_copy, dim=None, ignore=True, in_place=True)

        # plot 
        plt.figure(figsize=(10,5))
        ax = iplt.plt.axes(projection=iplt.ccrs.PlateCarree())
        plot = iplt.pcolormesh(sites_cube,cmap='YlGnBu', vmin=0, vmax=len(rundirs)+2)
        ax.coastlines(resolution='50m', color='black', linewidth=1)
        ax.add_feature(iplt.cartopy.feature.BORDERS)
        ax.gridlines()
        ax.set_xlim(lon_min,lon_max)
        ax.set_ylim(lat_min,lat_max)
        sl_plot = plt.scatter(source_locations_x, source_locations_y, c='r', marker=(5, 1),s=50,linewidth=0.5)
        plt.gca().coastlines()
        cbar = plt.colorbar(plot, ticks = np.arange(len(rundirs)+2))
        cbar.set_label('Plume overlap (number of sites)')

        plt.title('Plume overlap between '+str(N)+' sites')
        plt.savefig('Plot_plume_overlap_between_'+str(N)+'_sites_for_sandpit.png')
        plt.close()	



