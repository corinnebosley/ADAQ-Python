#!/usr/bin/env python

"""
Test script to plot overlap of daily plume extent combining multiple site plumes.

Designed to plot plume extent for each day (from multiple sites at any stage of their back run), then overlay daily plume extents for entire period to indicate hotspot region of most likely source location in terms of duration of plume presence.

"""
from __future__ import print_function

from six.moves.builtins import str
from six.moves.builtins import range
import iris
import iris.coords
import numpy as np
import iris.plot as iplt
import iris.quickplot as qplt
import matplotlib
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import matplotlib.colors as mcols
import glob
import os
from datetime import date
import pandas as pd


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
    
    # set period we want to consider (whole period for which we have results)

    start = date(2017,1,6)  
    end = date(2017,2,7)  
    daterange = pd.date_range(start,end) # set hourly interval if wanting more frequently than every 24h

    workdir = '/data/users/apdg/python_sample_data/name_plume_overlap/'

    rundirs = sorted(glob.glob(workdir+'BACK_RUN_*'))
    
    if len(rundirs) == 0:
            raise IOError('No back run directories found')
    
    source_locations_x = []
    source_locations_y = []
    
    for run in rundirs:
        sl_x, sl_y = get_source_location(run)
        source_locations_x.append(sl_x)
        source_locations_y.append(sl_y)

    constraint = set_latlon_constr(lat_min, lat_max, lon_min, lon_max)
    
    height_idx = 0

    overlap_cube = None
    
    # for each timestep (day)
    
    for DATE in daterange:
        DATE_str = str(DATE)
        fname_date = DATE_str[:4]+DATE_str[5:7]+DATE_str[8:10]  # convert date to format of output file name
        fsearch = glob.glob(workdir+'BACK_RUN_*/Fields_grid1_*'+fname_date+'*.txt') # return all fields files for the date
	
        # if there are one or more fields files returned, add together and plot
        if len(fsearch) == 0:
            continue
        print('{} files found for {}'.format(len(fsearch), DATE)) 
        # load all site plumes for that timestep
        timestep_cl = iris.load(fsearch, constraint)

        timestep_cube = None
        
        # sum the cubes from the different sites on that day
        for cube in timestep_cl:

            for coord in cube.coords():
                if coord.name() == 'height':
                    cube = cube[height_idx]
            if cube.ndim > 2:
                raise AssertionError('Cube has too many dimensions!')
                
            if timestep_cube == None:
                timestep_cube = cube
            else:
                timestep_cube = timestep_cube + cube
        
        # set all values to 1 (normalise)
        
        binary_cube = timestep_cube.copy()
        
        for i in range(timestep_cube.coord('latitude').shape[0]):
            for j in range(timestep_cube.coord('longitude').shape[0]):
                if timestep_cube.data[i][j] > 0.0:
                    binary_cube.data[i][j] = 1.0
                else:
                    binary_cube.data[i][j] = 0.0
        
        # gather all the binary cubes at each timestep
        if overlap_cube == None:
            overlap_cube = binary_cube
            n=1
        else:
            overlap_cube = overlap_cube + binary_cube
            n+=1
    # n is number of overlaps e.g. timesteps for which there is any fields file
    # n is used to define colorbar lims - if too small, plot looks bad
    if n < 5:
        n = 5

    # plot end result of all combined
    plt.figure(figsize=(10,5))
    ax = iplt.plt.axes(projection=iplt.ccrs.PlateCarree())
    plot = iplt.pcolormesh(overlap_cube,cmap='YlGnBu', vmin=0, vmax=n)
    ax.coastlines(resolution='50m', color='black', linewidth=1)
    ax.add_feature(iplt.cartopy.feature.BORDERS)
    ax.gridlines()
    ax.set_xlim(lon_min,lon_max)
    ax.set_ylim(lat_min,lat_max)
    plt.gca().coastlines()
    cbar = plt.colorbar(plot)
    cbar.set_label('Plume overlap (days)')
    sl_plot = plt.scatter(source_locations_x, source_locations_y, c='r', marker=(5, 1),s=50,linewidth=0.5)

    plt.title('Overlap of daily plume extent '+str(daterange[0])[:10]+' to '+str(daterange[-1])[:10])
    plt.savefig('Plot_all_sites_daily_plume_overlap_for_sandpit.png')
    plt.close()

