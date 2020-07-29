#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-
#
'''
PlumeandMet.py
Author: Frances Beckett
Last Modified: 29 July 2014

Program for time integrated plume position represented as contours. 

Overlaying field grid met data. 

Loops through mutliple files with different output times. 
 
Options to vary:

   * Contour levels
   * Extent of axes
   * Lat and Long Grid lines
   * Ticks on ColorBar
'''

from six.moves.builtins import str
import iris
import iris.plot as iplt
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import glob
import datetime
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

UTC_format = '%H%M%Z %d/%m/%Y'


def PlumeandMet(workdir, metDir):

    times = ['201303060600'] 

    for time in times:
        filenames = glob.glob(workdir+'*'+time+'.txt')
        filename = filenames[0] 
        attConstraint = iris.AttributeConstraint(Name='TotalAC')
        cube=iris.load_cube(filename,attConstraint)
        conc = cube 

        filename=metDir+'*'+time+'.txt' 
        precip = iris.load_cube(filename)

        colorscale = ('#ffffff', '#b4dcff','#04fdff','#00ff00','#fdff00','#ffbd02','#ff6a00','#fe0000', '#0000FF', '#800080', '#008000')

        # Set up axes
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.set_extent([-23, -13, 63, 67])

        # Set up country outlines
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
        gl.xlocator = mticker.FixedLocator([-23,-22,-21,-20,-19,-18,-17,-16,-15,-14]) 
        gl.ylocator = mticker.FixedLocator([63,64,65,66,67]) 
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER

        # Plot
        cf1 = iplt.contourf(precip,
                    levels=[0.0, 0.01, 0.1, 1.0, 10],colors=colorscale)

        cf = iplt.contour(conc,
                  levels=[1e-9, 3.16e-8, 1e-8, 3.16e-7, 1e-7, 3.16e-6, 1e-6, 3.16e-5, 1e-5]) 

        cb = plt.colorbar(cf1, orientation='horizontal',shrink=0.9)
        cb.set_label(str(precip.units))
        plt.title('Precipitation (coloured) and Air Concentration (contoured) \n' +time, fontsize=12)

        plt.show()
