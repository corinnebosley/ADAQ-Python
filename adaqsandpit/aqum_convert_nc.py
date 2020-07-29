#!/usr/bin/env python
"""
aqum_convert_nc.py -
Simple routine to use adaqcode for converting aqum pp files to netcdf
"""
from __future__ import print_function

from six.moves.builtins import zip
import os
import sys
#Find adaq directory -
# should be in directory above the one containing this file
adaq_path = os.path.dirname(os.path.realpath(__file__))+'/../'
sys.path.append(adaq_path)
import adaqcode as adaq 
import iris
import datetime


if __name__ == "__main__":

   print('Current version of python:',sys.version)
   print('Current version of iris:',iris.__version__)
   print('Beginning at ',datetime.datetime.now())
   Verbose=False #Print extra output

   #--------------------------------------------------
   # Get ini data - the following will work via cron...
   if len(sys.argv) > 1 :
      inifilename = sys.argv[1]
   else:
      #For testing, using plot_test.ini
      try:
         curr_path = os.path.dirname(os.path.realpath(__file__))
      except:
         #For running interactively
         curr_path='./'
      inifilename=curr_path + '/aqum_convert_nc.ini'   

   ini = adaq.Inifile ( inifilename,  setdefaults=True  )
   ini_data = ini.get ()

   #--------------------------------------------------
   #Get end and start dates. End date must be midnight of current date

   end_datetime = ini_data['end_datetime'] 
   start_datetime = ini_data['start_datetime']

   #--------------------------------------------------
   # Get model data
   for model, datadir in zip(ini_data['models_list'], ini_data['models_dir_list']):

      print("Getting model data for ",model," at ",datetime.datetime.now())

      #--------------------------------------------------
      # Get a list of pp filenames between the specified times

      if Verbose : print("Generating filenames at ",datetime.datetime.now())

      ppfiles = adaq.PPFiles ( datadir, start_datetime, end_datetime,
                             forecast_day=ini_data['forecast_day'] )
      filenames = ppfiles.get ()

      #--------------------------------------------------
      # Read model data
      if Verbose :
         print("Creating model field cubes at ",datetime.datetime.now())
      md = adaq.PPData ( filenames=filenames,
                         short_name_list=ini_data['short_name_list'],
                         start_datetime=start_datetime,
                         end_datetime=end_datetime,
                         label=model)
      md.readdata()

      print(md)

      #--------------------------------------------------
      # Convert to ug/m3
      for cube in md.gridded_cube_list:
         md.convert(cube)

      #--------------------------------------------------
      #Split into daily files
      date = start_datetime.date()
      while date <= end_datetime.date():
         print(date)
         datestr = date.strftime("%Y%m%d")
 
         tconstraint = iris.Constraint(time=lambda cell:
                                       cell.point.date() == date )

         gcl = iris.cube.CubeList()
         for cube in md.gridded_cube_list:
            newcube = cube.extract(tconstraint)
            if newcube != None:
               gcl.append(newcube)
         if len(gcl) > 0 :
            # Save to file
            md.save_gridded(gridded_cube_list = gcl,
                            filename = datadir+'/'+model+'_'+datestr+'.nc')


         date += datetime.timedelta(days=1)






















      
