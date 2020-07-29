#!/usr/bin/env python
#
#  Tephi like plot of met data produced by VAAC for various plume rise schemes
#    by Ayoe Buus Hansen, Met Office, ADAQ
#     
#   Last updated: October, 2015
#______________________________________________________________________________

'''
tephibarbs.py
Author: Ayoe Buus Hansen
Last Modified: December, 2015

Tephi like plot of met data produced by VAAC for various plume rise schemes
Produces two plots to capture both the lower and the upper part of the atmosphere.
Reads from Met file:

  * height, wind speed, wind direction, pressure, temperature (C), and humidity
  
Options to vary
  * min and max pressure
  * min and max theta
'''
from __future__ import division

#Import python libraries
from six.moves.builtins import zip
import six

import numpy as np
import os.path
#import tephi
import datetime, os, glob, linecache, re, argparse, csv
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.dates as mdates
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.backends.backend_pdf import PdfPages

if __name__ == '__main__':

    if six.PY3:
        raise ValueError('Not implemented for Python 3: tephi module not available')
    else:
        import tephi
        
    pp = PdfPages('TephiAtVent.pdf')
    #______________________________________________________________________________
    # SET DATA PATH
    #______________________________________________________________________________
    met_filename = '/net/home/h04/ahansen/MyResearch/VAACsprint/MetData.txt'

    height       = np.loadtxt(met_filename, unpack = True, dtype = float, delimiter = ',', skiprows = 38, usecols=(3,))
    windspeed    = np.loadtxt(met_filename, unpack = True, dtype = float, delimiter = ',', skiprows = 38, usecols=(4,))
    direction    = np.loadtxt(met_filename, unpack = True, dtype = float, delimiter = ',', skiprows = 38, usecols=(5,))
    pressure     = np.loadtxt(met_filename, unpack = True, dtype = float, delimiter = ',', skiprows = 38, usecols=(6,))
    Ctemperature = np.loadtxt(met_filename, unpack = True, dtype = float, delimiter = ',', skiprows = 38, usecols=(7,))
    humidity     = np.loadtxt(met_filename, unpack = True, dtype = float, delimiter = ',', skiprows = 38, usecols=(8,))
    hPressure    = np.loadtxt(met_filename, unpack = True, dtype = float, delimiter = ',', skiprows = 38, usecols=(6,))
    hPressure    = hPressure/100.

    HeightTemp = list(zip(hPressure, Ctemperature))
    column_titles = [('pressure', 'temperature')]
    tephi.MIN_PRESSURE = 70 
    tephi.MAX_PRESSURE = 1000
    tephi.ISOBAR_SPEC = [(50, 0.50), (100, 1.5), (200, None)]
    # tephi.ISOBAR_SPEC = [(25, 0.45), (50, 0.50), (100, 1.5), (200, None)]
    tephi.MIN_THETA = 0 
    tephi.MAX_THETA = 700
    tephi.MIN_WET_ADIABAT = 0
    tephi.MAX_WET_ADIABAT = 200
    tephi.WET_ADIABAT_SPEC = [(5, None)]
    tephi.MIXING_RATIO_SPEC = [(1, 0.05), (2, 0.18), (4, 0.3), (8, 1.5)]
    #tephi.MIXING_RATIO_FIXED = [0.001, 0.01, 0.05, 0.1, 0.5, 1.0, 5.0, 10.0, 20.0, 50.0, 75.0,100.0,150.0]
    tpg = tephi.Tephigram(anchor=[(1100, 0), (70, 0)], isotherm_locator=tephi.Locator(25), dry_adiabat_locator=tephi.Locator(25))
    profile = tpg.plot(HeightTemp, label='Temperature', color='blue', linewidth=2, linestyle='-', marker='')
    Wind = list(zip(windspeed, direction, hPressure))
    profile.barbs(Wind[::15], length=8, pivot='middle', color='green', linewidth=3, gutter=0.15)
    # profile.barbs(Wind[], length=8, pivot='middle', color='red', linewidth=3, gutter=0.15)
    ## export image
    plot_out1 = ('VAAC_Tephi_Barbs_plot_full.png')
    plt.savefig(pp, format='pdf')
    plt.savefig(plot_out1, bbox_inches='tight', dpi = 150)
    plt.show()
    plt.clf()
    #______________________________________________________________________________
    HeightTemp = list(zip(hPressure, Ctemperature))
    column_titles = [('pressure', 'temperature')]
    tephi.MIN_PRESSURE = 500
    tephi.MAX_PRESSURE = 1000
    tephi.ISOBAR_SPEC = [(50, 0.50), (100, 1.5), (200, None)]
    tephi.MIN_THETA = 0 
    tephi.MAX_THETA = 100
    tephi.MIN_WET_ADIABAT = 0
    tephi.MAX_WET_ADIABAT = 50
    tephi.WET_ADIABAT_SPEC = [(5, None)]
    tpg2 = tephi.Tephigram(anchor=[(1000, 0), (500, 0)], isotherm_locator=tephi.Locator(5), dry_adiabat_locator=tephi.Locator(5))
    profile = tpg2.plot(HeightTemp, label='Temperature', color='blue', linewidth=2, linestyle='-', marker='')
    Wind = list(zip(windspeed, direction, hPressure))
    profile.barbs(Wind[::5], length=8, pivot='middle', color='green', linewidth=3, gutter=0.15)
    ## export image
    plot_out2 = ('VAAC_Tephi_Barbs_plot_focus.png')
    plt.savefig(pp, format='pdf')
    plt.savefig(plot_out2, bbox_inches='tight', dpi = 150)
    #plt.show()
    pp.close()
    #______________________________________________________________________________
    # END OF CODE
    #______________________________________________________________________________
