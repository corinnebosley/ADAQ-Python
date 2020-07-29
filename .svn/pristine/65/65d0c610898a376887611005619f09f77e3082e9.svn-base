#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-
# 
# Python script to plot max plume heights
#

'''
MaxHeight.py
Author: Frances Beckett
Last Modified: 28 July 2014

Program to plot the maximum height of a plume. 
Loops through mutliple files with different output times.  
Options to vary

  * Height Levels
  *  Number of Height Levels
  *  Extent of axes
  *  Lat and Long Grid lines
  *  Ticks on ColorBar
  *  Contour levels
  *  Threshold Air Concentration
'''
from __future__ import print_function


import iris
import iris.plot as iplt
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import glob
import datetime
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

UTC_format = '%H%M%Z %d/%m/%Y'

def MaxHeight(workdir):

    times = ['201309160100', '201309160900'] 
 

    for time in times:
        filenames = glob.glob(workdir+'*'+time+'.txt')
        print(filenames)
        filename = filenames[0] 

        # Identify cube attributes
        cubes = iris.load(filename)
        print(cubes[2].attributes)

        names=['Level1', 'Level2', 'Level3', 'Level4', 'Level5', 'Level6', 'Level7', 'Level8', 'Level9', 'Level10']
        
        colorscale = ('#b4dcff','#04fdff','#00ff00','#fdff00','#ffbd02','#ff6a00','#fe0000', '#0000A0', '#800080', '#006400')


        # Using contourf to provide my colorbar info, then clearing the figure
        Z = [[0,0],[0,0]]
        levels = (0,200,400,600,800,1000,1200,1400,1600,1800,2000)
        CS3 = plt.contourf(Z, levels, colors=colorscale)
        plt.clf()


        # Set up axes
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.set_extent([-24, -14, 60, 67])



        for name in names:
            print(name)
            attConstraint = iris.AttributeConstraint(Name=name)
            cube=iris.load_cube(filename,attConstraint) 


            # Mask data less than threshold which represents 1 real particle in /m3
            conc = cube 
            conc.data = np.ma.masked_less(conc.data, 1e-12)

              

            # Identify time
            phenom_time = conc.coord('time')
            phenom_time_date1 = phenom_time.units.num2date(phenom_time.bounds[0][0]).strftime(UTC_format)
            phenom_time_date2 = phenom_time.units.num2date(phenom_time.bounds[0][1]).strftime(UTC_format)



            # Plot Data
            if name == 'Level1':
                colorscale = '#b4dcff'
            if name == 'Level2':
                colorscale = '#04fdff'
            if name == 'Level3':
                colorscale = '#00ff00'
            if name == 'Level4':
                colorscale = '#fdff00'
            if name == 'Level5':
                colorscale = '#ffbd02'
            if name == 'Level6':
                colorscale = '#ff6a00'
            if name == 'Level7':
                colorscale = '#fe0000'
            if name == 'Level8':
                colorscale = '#0000A0'
            if name == 'Level9':
                colorscale = '#800080'
            if name == 'Level10':
                colorscale = '#006400'
            cf1=iplt.contourf(conc,colors=colorscale) 


    
        
        # Add county outlines
        countries = cfeature.NaturalEarthFeature(
          category='cultural',
          name='admin_0_countries',
          scale='10m',
          facecolor='none')
        ax.add_feature(countries, edgecolor='black',zorder=2)

        # Set-up the gridlines
        gl = ax.gridlines(draw_labels=True, 
                          linewidth=0.8, 
                          alpha=0.9)    
        gl.xlabels_top = False
        gl.ylabels_right = False
        gl.xlocator = mticker.FixedLocator([-24,-23,-22,-21,-20,-19,-18,-17,-16,-15,-14]) # lat long lines
        gl.ylocator = mticker.FixedLocator([60,62,64,66,68]) #lat and long lines
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER


        #Set-up the Colour Bar        
        cb = plt.colorbar(CS3)
        cbartext='Height m asl'     
        cb.set_label(cbartext)
        cb.ax.set_xticklabels([200,400,600,800,1000,1200,1400,1600,1800,2000], rotation = 'vertical')
        plt.title(phenom_time_date2, fontsize=12)


        plt.show()

