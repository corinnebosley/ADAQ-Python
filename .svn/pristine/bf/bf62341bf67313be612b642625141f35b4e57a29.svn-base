#!/usr/bin/env python
'''Script to reformat metdb land synops observations from raw metdb output,
   into format to match AURN reformatted obs.

   Input arguments:
      1. metdb output filename
      2. output directory to contain reformatted observations and sitelist.txt

  .. todo ::
      Convert this to a class of its own so there is no longer a need to
      reformat data, but rather this is read in as a ADAQData object.
'''     
from __future__ import division
from __future__ import print_function

from six.moves.builtins import str

import numpy as np
import sys
import os
from datetime import datetime,timedelta

Nan=np.nan

def read_data(filein):
    '''Reads in data from metdb out file
       Returns dictionary containing site names as keys and associated data'''
    convert=lambda x: np.nan if float(x) == -9999999.00 else float(x)
    data=np.genfromtxt(filein,dtype=None,names=True,delimiter=[25,11,11,5,3,3,3,3,15,15,15,15,15],
                       converters={'WindDir':convert,
                                   'WindSpeed':convert,
                                   'Temperature':convert,
                                   'DewpointTemp':convert,
                                   'PMSL':convert}) #Read in from file 

    sitesdict={}

    for line in data:
        #print line['SiteName']
        if line['SiteName'] not in sitesdict: 
            sitesdict[line['SiteName']]=[]
        sitesdict[line['SiteName']].append(line)
    #sitesdict['KESWICK                  '][0]['Lat'] to access the data

    return sitesdict

def generate_abbrevs(sitesdict):
    '''Generates unique site abbreviations for each site returns as a dictionary'''
    siteabbrevs={}
    for site in sitesdict.keys():
        if site != '                         ':
            abbrev=site[0:4]
            abbrev=abbrev.replace(' ','')
            abbrev=abbrev.replace('.','')
            abbrev=abbrev.replace('-','')
            count=0
            while abbrev in siteabbrevs:
                if count == 0:
                    abbrev=abbrev[0:3]
                else:
                    abbrev=abbrev[0:3]+str(count)
                count+=1
            abbrev='m'+abbrev
            siteabbrevs[abbrev]=site
    return siteabbrevs




def generate_sitelist(sitesdict,siteabbrevs,sitelistfile,directory='obs'):
    '''Generate site list with abbreviations,lat,lon and full site name'''
    print('generating site list of ',len(siteabbrevs),'sites')
    if not os.path.exists(directory): os.mkdir(directory)
    with open(directory+'/'+sitelistfile,"w") as fout: 
        fout.write('GEMS_code  Abbrev    Lat        Lon      Altit        Type      NAME\n')

        for abbrev,site in siteabbrevs.items():
            lat=sitesdict[site][0]['Lat']
            lon=sitesdict[site][0]['Lon']
            sitename=site.strip()
            sitename=sitename.replace(' (','_')
            sitename=sitename.replace(', ','_')
            sitename=sitename.replace(' ','_')
            sitename=sitename.replace('(','_')
            sitename=sitename.replace(')','')
            sitename=sitename.replace('/','_')
            sitename=sitename.replace(',','_')
            #print sitename
            outputstr="%-8s%-5s%12.8f %13.9f    0  MET               %-25s\n" % (abbrev,abbrev,lat,lon,sitename)
            fout.write(outputstr)

            
def calc_rh(T,DT):
    '''Calculate relative humidity from Temperature(T), Dewpoint temperature(DT) in kelvin'''
    RH = 100*np.exp(19.8*(DT-273)/DT)/(np.exp(19.8*(T-273)/T))
    return RH


def output_data(sitesdict,siteabbrevs,directory='obs'):
    '''NOTE: Output data in format that can be read in by verification code'''
    if not os.path.exists(directory): os.mkdir(directory)
    print('outputting data') 
    str_format='%-11s%-13s%11.1f%11.1f%11.1f%11.1f%11.1f   \n'
    for abbrev,site in siteabbrevs.items():
        #print abbrev,site
        firstline=sitesdict[site][0]
        firstdate=str(firstline['year'])+str('%02d' % firstline['mm'])+str('%02d' % firstline['dd'])
        lastline=sitesdict[site][-1]
        lastdate=str(lastline['year'])+str('%02d' % lastline['mm'])+str('%02d' % lastline['dd'])
        filename=directory+'/'+abbrev+'_'+firstdate+'_'+lastdate+'.txt'
        with open(filename,"w") as fout:
            fout.write('locationID date     time     T1.5_C   WDir_deg    WSpd_ms  Press_hPa       rh_%\n')    
            for iline,line in enumerate(sitesdict[site]):
                dt=datetime(line['year'],line['mm'],line['dd'],line['hh'])
                if iline != 0:
                    #skip any spot observations and just use hourly means
                    if dt == dt_prev or line['mn'] !=0:
                        print('spot obs found but skipping')
                        continue
                    #infill missing datestamps with NANs    
                    if dt != dt_prev+timedelta(hours=1):
                        print('missingtime')
                        dt_iter=dt_prev+timedelta(hours=1)
                        while dt_iter != dt:
                            dt_str_tmp=datetime.strftime(dt_iter,"%Y%m%d %H")
                            dt_hr = int(dt_str_tmp[9:11])
                            if dt_hr == 0:  
                                date_prev=dt_iter-timedelta(hours=1)    
                                dt_str_tmp=datetime.strftime(date_prev,"%Y%m%d")+' 24'
                            else:
                                dt_str_tmp=datetime.strftime(dt_iter,"%Y%m%d %H")
                            output_str=str_format % (abbrev,dt_str_tmp,Nan,Nan,Nan,Nan,Nan)
                            print(dt_str_tmp)
                            fout.write(output_str)
                            dt_iter+=timedelta(hours=1)
                if line['hh'] == 0:
                    date_prev=dt-timedelta(hours=1)
                    dt_str=datetime.strftime(date_prev,"%Y%m%d")+' 24'
                    line['hh'] = 24
                    line['dd'] = int(dt_str[6:8])
                    line['mm'] = int(dt_str[4:6])
                    line['year'] = int(dt_str[0:4])
                else:
                    dt_str=datetime.strftime(dt,"%Y%m%d %H")
                if iline == 0:
                    #Fill in missing data at the beginning of file to start from 1Z
                    if line['hh'] != 1:
                        for hh in np.arange(1,line['hh']):
                            dt_str_tmp=dt_str[:8]+" %02d" % hh
                            output_str=str_format % (abbrev,dt_str_tmp,Nan,Nan,Nan,Nan,Nan)
                            fout.write(output_str)
                rh=calc_rh(float(line['Temperature']),float(line['DewpointTemp']))
                if np.isnan(rh): rh=Nan
                temp=float(line['Temperature'])-273.15
                if np.isnan(temp): temp=Nan
                wind_dir=float(line['WindDir'])
                if np.isnan(wind_dir): wind_dir=Nan
                wind_sp=float(line['WindSpeed'])
                if np.isnan(wind_sp): wind_sp=Nan
                press=float(line['PMSL'])/100.
                if np.isnan(press): press=Nan
               
                output_str=str_format % (abbrev,dt_str,temp,wind_dir,wind_sp,press,rh)
                fout.write(output_str)

                dt_prev=dt
            #fill in missing data at end of file to reach 24Z 
            if line['hh'] != 24:
                for hh in np.arange(line['hh']+1,25):
                    dt_str_tmp=dt_str[:8]+" %02d" % hh
                    output_str=str_format % (abbrev,dt_str_tmp,Nan,Nan,Nan,Nan,Nan)
                    fout.write(output_str)
                            
        #sitesdict['KESWICK                  '][0]['Lat'] to access the data
        #break
        #with open(sitelistfile,"w") as fout:
                    
if __name__ == '__main__':
    '''Reformats metdb obs and also outputs sitelist.txt'''
    
    filein=sys.argv[1] #metdb output filename
    outdir=sys.argv[2] #output directory
    sitesdict=read_data(filein)
    siteabbrevs=generate_abbrevs(sitesdict)
    generate_sitelist(sitesdict,siteabbrevs,sitelistfile='sitelist_met.txt',directory=outdir)
    output_data(sitesdict,siteabbrevs,directory=outdir)

