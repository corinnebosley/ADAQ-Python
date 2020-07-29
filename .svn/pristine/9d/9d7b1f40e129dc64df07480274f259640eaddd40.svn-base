#!/usr/bin/env python
"""
Script to retrieve AURN data from website and reformat to fixed-format files.
Required dates, short_names, sites and directories should be setup in
inifile, for example see adaqscripts/retrieve_aurn.ini.
Fixed format files will be stored in obs_dir/year_fixed_format/abbrev_year.txt.

From the 'adaqscripts' directory, run via:

.. code-block:: ksh

    ./retrieve_aurn.py [inifilename]

where inifilename includes the full path to the ini file.

Alternatively, from within other python code, run using:

.. code-block:: python

    from adaqscripts import retrieve_aurn
    aq_plot.retrieve_reformat_aurn(inifilename='inifilename')

"""
from __future__ import print_function
from six.moves.builtins import str
import os
import sys
import datetime

import iris

#Find adaq directory -
# should be in directory above the one containing this file
adaq_path = os.path.dirname(os.path.realpath(__file__))+'/../'
sys.path.append(adaq_path)
import adaqcode

def retrieve_yearly_web(year, sites_data, yearly_dir):
    """
    Retrieve yearly observation file from AURN website, UK-AIR.
    This retrieves all available species.

    :param year: float, year to retrieve
    :param sites_data: numpy ndarray containing site information data
                       from a :class:`sites_info.SitesInfo` object.
    :param yearly_dir: string, directory to retrieve into

    """

    if not os.path.isdir(yearly_dir):
        print('Creating output yearly directory:', yearly_dir)
        os.makedirs(yearly_dir)

    for abbrev in sites_data['abbrev']:
        print('Retrieving ' + abbrev + '...')
        url = 'http://uk-air.defra.gov.uk/datastore/data_files/site_data/'
        filename = abbrev + '_' + str(year) + '.csv'
        #Delete pre-existing file first
        full_filename = yearly_dir+'/'+filename
        if os.path.isfile(full_filename):
            os.remove(full_filename)
        #Issue wget command to retrieve.
        #Note -q is used to quiet the output - if trying to debug, then
        #remove this option!
        cmd = ('/usr/bin/wget -q --no-check-certificate --directory-prefix='
               + yearly_dir + ' ' + url + '/' + filename)
        return_code = adaqcode.shell_commands.call_shell(cmd)
        if return_code == 0:
            print('File retrieved: ' + full_filename)
        else:
            print('Error retrieving: ' + filename + ' (file may not exist)')



def retrieve_reformat_aurn(inifilename=None, ini_dict=None, sites_data=None,
                           download_csv=True, convert_ff=True):
    """
    Routine to read inifile, download AURN csv files from UK-AIR and
    convert to fixed-format files.

    Csv files are downloaded into obs_dir/year_csv/abbrev_year.txt.
    These will contain all available species for the entire year.

    Fixed format files are written to obs_dir/year_fixed_format/abbrev_year.txt
    These will contain only the required species and the required dates.

    :param inifilename: String giving filename of ini file.
    :param ini_dict: Dictionary of a :class:`inifile` object. If this is given,
                     then it is used instead of reading in from inifilename.
    :param sites_data: numpy ndarray containing site information data
                       from a :class:`sites_info.SitesInfo` object.
                       If this is given, then it is used instead of reading in
                       from a file/cube.
    :param download_csv: Logical. If True, then downloads csv files from UK-AIR
                         Can be set to False if already downloaded into
                         ini_dict['obs_dir'] + '/' + str(year) + '_csv'.
                         (May also be required for doctests if no web access
                         available?)
    :param convert_ff: Logical. If True then converts csv files into
                       fixed-format files. Can be set to False if these
                       files are not required.

    :returns: (ini_dict, sites_data, od)

     * **ini_dict** - Dictionary of an :class:`inifile` object.
     * **sites_data** - numpy ndarray containing site information data
       from a :class:`sites_info.SitesInfo` object.
     * **od** - AURN observations from :class:`aurn_data` for the
       required dates.

    >>> ini_dict, sites_data, od = retrieve_reformat_aurn() # doctest: +ELLIPSIS
    Current version of python: ...
    Current version of iris: ...
    Beginning at ...
    Reading inifile .../retrieve_aurn.ini
    Number of sites:  5
    Year: 2017...
    Retrieving ABD...
    File retrieved: .../2017_csv/ABD_2017.csv
    Retrieving ACTH...
    File retrieved: .../2017_csv/ACTH_2017.csv
    Retrieving AH...
    File retrieved: .../2017_csv/AH_2017.csv
    Retrieving HAR...
    Error retrieving: HAR_2017.csv (file may not exist)
    Retrieving YW...
    File retrieved: .../2017_csv/YW_2017.csv
    Reading obs data files
    Found obs for  .../2017_csv/ABD_2017.csv
    Found obs for  .../2017_csv/ACTH_2017.csv
    Found obs for  .../2017_csv/AH_2017.csv
    No file found for HAR
    Found obs for  .../2017_csv/YW_2017.csv
    Creating obs cubes
    Creating missing data cubes
    Written to .../2017_fixed_format/ABD_2017.txt
    Written to .../2017_fixed_format/ACTH_2017.txt
    Written to .../2017_fixed_format/AH_2017.txt
    HAR - all nan data, no file output
    Written to .../2017_fixed_format/YW_2017.txt
    Finished retrieving AURN files

    Check the output of one of these files:

    >>> outputfile = ini_dict['obs_dir']+'/2017_fixed_format/ACTH_2017.txt'
    >>> with open(outputfile,"r") as fin:
    ...     for iline, line in enumerate(fin):
    ...         print(line, end='')
    ...         if iline > 3: break
    locationID date     time       O3_ug/m3 PM2p5_ug/m3
    ACTH       20170601 00         69.74972     3.70000
    ACTH       20170601 01         69.18427     3.90000
    ACTH       20170601 02         61.46756     6.00000
    ACTH       20170601 03         58.82326     3.80000

    """

    print('Current version of python:', sys.version)
    print('Current version of iris:', iris.__version__)
    print('Beginning at ', datetime.datetime.now())

    #Get ini_dict
    if ini_dict is None:
        #Read in inifile
        ini_dict = adaqcode.inifile.get_inidict(
            inifilename=inifilename,
            defaultfilename='adaqscripts/retrieve_aurn.ini')

    #Get list of sites, by reading in site information data file
    if sites_data is None:
        sites = adaqcode.sites_info.SitesInfo()
        sites_data = sites.read_from_file(ini_dict['sites_file'],
                                          allsites=True)

    #Generate set of required years
    years = set()
    date = ini_dict['start_datetime']
    while date <= ini_dict['end_datetime']:
        years.add(date.year)
        date += datetime.timedelta(days=1)

    yearly_csv_dir_list = [] #list of all years directories
    for year in years:
        print('Year:', year)
        #Retrieve csv files from website
        yearly_csv_dir = ini_dict['obs_dir'] + '/' + str(year) + '_csv'
        yearly_csv_dir_list.append(yearly_csv_dir)
        if download_csv:
            retrieve_yearly_web(year, sites_data, yearly_csv_dir)

    #Convert to fixed-format files:
    #Read data in
    od = adaqcode.aurn_data.AURNData()
    od.readdata(obsdir_list=yearly_csv_dir_list,
                short_name_list=ini_dict['short_name_list'],
                sites=sites_data,
                start_datetime=ini_dict['start_datetime'],
                end_datetime=ini_dict['end_datetime'])
    #Convert to ff file
    if convert_ff:
        od.write_ff_files(ini_dict['obs_dir'])


    print('Finished retrieving AURN files')

    return ini_dict, sites_data, od



if __name__ == '__main__':

    retrieve_reformat_aurn()

    #ini_dict, sites_data, od = retrieve_reformat_aurn(download_csv=False, convert_ff=True)

    #import doctest
    #doctest.testmod()
