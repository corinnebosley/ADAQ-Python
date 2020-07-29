'''
Routine to read in pp field and plot contour of this data with AURN
observations plotted over the top for the same time and same colour scale.
Also makes use of user-defined DAQI colour scale.

Lucy Davis, Paul Agnew, Bryon Blay, September 2013

.. note::
   from pylab import * has been commented out as it causes the documentation to
   fail.
   Susan, Aug 2015
'''
from __future__ import division
from __future__ import print_function

from six.moves.builtins import str
from six.moves import cPickle as pickle

import iris
import numpy as np
import iris.quickplot as qplt
import iris.plot as iplt
import cartopy.crs as ccrs
import matplotlib.colors as mplcolors
import matplotlib.pyplot as plt
#from pylab import *
import glob

def readAURNffobs(cube_dt,species,sitelist=None,obsdir=None):
    '''Read in observations from reformatted AURN text file'''

    if sitelist == None:
        sitelist=np.genfromtxt('/home/h03/apdg/AQUM/AQcases_fcm/code/aq_sites_all.txt',dtype=None,names=True)
    if obsdir == None:
        obsdir='/data/nwp1/apdl/AQobs/'+str(cube_dt.year)+'_reformatted/'

    lats=[]
    lons=[]
    values=[]
    #Loop through sites
    for iabbrev,abbrev in enumerate(sitelist['Abbrev']):
        files=glob.glob(obsdir+abbrev+'_*')
        #Check if file exists - extract data if so
        if len(files) == 1:
            obsdata=np.genfromtxt(files[0],dtype=None,names=True)
            #Find date that matches cube
            for idata, data in enumerate(obsdata['date']):
                date=str(obsdata['date'][idata])
                time=str(obsdata['time'][idata])
                dt=fixfmt2dt(date,time)
                if dt == cube_dt and np.isfinite(obsdata[species][idata]):
                    #Pull out data value if not nan
                    print(abbrev,sitelist['Lat'][iabbrev],sitelist['Lon'][iabbrev])
                    print(obsdata[species][idata])
                    lats.append(sitelist['Lat'][iabbrev])
                    lons.append(sitelist['Lon'][iabbrev])
                    values.append(obsdata[species][idata])
                    break   

    return lats,lons,values


def unit_conversion_cube(cube,outputunit,speciesname=None,inputunit=None):
    '''Routine to convert units in a cube'''

    if speciesname == None:
        speciesname=cube.standard_name
  
    if inputunit == None:
        inputunit=cube.units    
    
    conv,finalunit,conversionfactors_dict = unit_conversion(inputunit,outputunit,speciesname)

    cube.data=cube.data*conv
    cube.units=finalunit                   
    return cube


def readcube(stash='m01s34i001',ppfile=None,hour=13):
    '''Read in cube of pp data for single stash code
    Note - currently only plots first time in file
    '''
    if ppfile == None:
        ppfile='/project/NAME/apnm/NAMEIIIout/PythonHackathon/DataFiles/ppData/aqum_20060614.pp'
    stashconstraint=iris.AttributeConstraint(STASH=stash)
    cube=iris.load(ppfile,stashconstraint)
    cube=cube[0] #Should only have loaded one cube, so pick cube out from list

    #Get cube time
    cube_time=cube.coord('time')
    cube_dt=cube_time.units.num2date(cube_time.bounds[:,1])

    #Find required hour
    hourfound=0
    for i,dt in enumerate(cube_dt):
        if cube_dt[i].hour == hour :
            cube=cube[i]
            hourfound=1
            cube_dt=cube_dt[i]
            break
    if hourfound == 0:
        print('Required hour not found in file (',hour,')')
        print('Available hours:',[ dt.hour for dt in cube_dt ])
        return
    
    print('Time:',cube_dt)

    #Convert units
    cube=unit_conversion_cube(cube,'1e-6 g / m3',speciesname='mass_fraction_of_ozone_in_air',inputunit='kg / kg')
    print('Max model data in converted units', cube.data.max())

    return cube, cube_dt

def convert_grids(cube,obs_lons,obs_lats):
    '''Convert observation lons and lats onto rotated model grid'''

    coord_sys=cube.coord(axis='x').coord_system.as_cartopy_crs() #get target coord system
    obs_coord_sys=ccrs.Geodetic() #regular lat-lon coord system

    #Convert obs lat-lon onto model grid
    result=coord_sys.transform_points(obs_coord_sys,np.asarray(obs_lons),np.asarray(obs_lats))
    lons_rot=result[:,0]
    lats_rot=result[:,1]

    return lons_rot,lats_rot
    

def DAQIcmap():
    '''Generate colour map based on DAQI colours'''
    colours=[[146,208,80 ],
            [0  ,176,79 ],
            [0  ,116,50 ],
            [250,240,0  ],
            [255,192,0  ],
            [227,108,10 ],
            [255,0  ,0  ],
            [192,0  ,0  ],
            [112,0  ,0  ],
            [111,47 ,160]]
    colours=np.array(colours)/255.
    cmap=mplcolors.ListedColormap(colours, name='from_list')

    return cmap

def plotdata(modelcube,obsvalues,obs_lons_rot,obs_lats_rot,
             levels=[0,33,66,100,120,140,160,180,200,220,250],
             title=None):
    '''Plot contour plot with observations overlayed'''
    #Set up contour levels and colours
    norm = mplcolors.BoundaryNorm(levels, len(levels))
    #cmap=plt.get_cmap('jet',len(levels))
    cmap=DAQIcmap()

    #Plot model data
    cs=iplt.contourf(modelcube,levels=levels,cmap=cmap,norm=norm)
    ax=qplt.plt.gca()

    #Overplot observations
    qplt.plt.scatter(obs_lons_rot,obs_lats_rot,c=obsvalues,cmap=cmap,norm=norm)

    #Make plot look nice
    cb=plt.colorbar() #plt.colorbar(cs) if not overplotting obs
    cb.set_label('$\mu gm^{-3}$')
    
    ax.coastlines('50m')
    ax.gridlines()
    
    ax.set_title(title)
    plt.tight_layout()

    qplt.plt.show()

def dt2fixfmt(dt):
    '''Routine to convert between datetime format and fixed format time-stamp
    Input:
    datetime.datetime(2012, 2, 26, 0, 0)
    Output:
    dates,times
    Where:
    dates='20120225'
    times='24'
    Notes:
    - Also does conversion from 00Z to 24Z previous day
    - If given a list as entries, then calls dts2fixmt looping routine instead
    '''       
    if not isinstance(dt,datetime):
        #Need to call looping routine instead
        dates,times=dts2fixfmt(dt)
        return dates,times
    else:    
        date=dt.date().strftime("%Y%m%d")
        time=dt.time().strftime("%H")
        if time == '00':
            date=(dt.date()-timedelta(days=1)).strftime("%Y%m%d")
            time='24'
        return date,time


def fixfmt2dt(date,time):
    '''Routine to convert fixed-format date-time stamp to datetime format
    Input:
    date='20120225'
    time='24:00'
    Output:
    datetime.datetime(2012, 2, 26, 0, 0)
    Notes:
    - Also does conversion from 24Z to 00Z next day
    - If given a list of entries, then calls aurn2dts instead
    '''       
    if not (isinstance(date,str) and isinstance(time,str)):
        #Need to call looping routine instead
        dts=fixfmts2dt(date,time)
        return dts
    else:    
        datedt=datetime.strptime(date,"%Y%m%d")
        time=time[0:2]
        if time == '24':
            datedt = datedt + timedelta(days=1)
            time='00'
        timedt=timedelta(hours=int(time))
        dt=datedt+timedt
        return dt


    
if __name__ == "__main__":
    
    species='O3'
    ppfile='/project/NAME/apnm/NAMEIIIout/PythonHackathon/DataFiles/ppData/aqum_20060614.pp'
    requiredhour=13
    rereadobs=1 #Read in obs from file, or set to zero to read in from previously generated pickle file

    #Read cube
    stash='m01s34i001'
    cube,cube_dt=readcube(stash=stash,ppfile=ppfile,hour=requiredhour)
    print(cube_dt)
    
    #Read in observations
    #- to save time when developing, pickle output so dont need to read everytime
    if rereadobs == 1:
        lats,lons,values=readAURNffobs(cube_dt,species+'_ugm3')
        with open('obsfile.pkl',"wb") as fp:
            obs=(lats,lons,values)
            pickle.dump(obs,fp)
    else:
        with open('obsfile.pkl',"rb") as fp:
            obs=pickle.load(fp)
        lats,lons,values=obs

    #Convert obs lon and lats onto rotated grid
    lons_rot,lats_rot = convert_grids(cube,lons,lats)

    #Determine plot title to include date-time stamp
    date,time=dt2fixfmt(cube_dt)
    title=species+' '+str(date)+' '+str(time)+'Z'

    #Plot model data
    plotdata(cube,values,lons_rot,lats_rot,title=title)


