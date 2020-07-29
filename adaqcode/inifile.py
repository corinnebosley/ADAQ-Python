"""
inifile.py - contains class for reading inifile and converting to dictionary
"""
#Import __future__ first
from __future__ import print_function
from __future__ import absolute_import

#Standard library imports
import ast
import os
import sys
import warnings
from datetime import datetime, timedelta

#Third party imports
from six.moves.builtins import str

#Local imports
import config
#from . import config

def get_filename(defaultfilename=None):
    """
    Find automatic filename:
      * If given as an argument to calling script, then set to this
      * Otherwise, uses :data:`config.CODE_DIR` + defaultfilename.

    >>> filename = get_filename(defaultfilename='adaqcode/inifile.ini')
    >>> print('The inifile is', filename) # doctest: +ELLIPSIS
    The inifile is .../adaqcode/inifile.ini

    """

    if len(sys.argv) > 1:
        #Get value of first argument given when run as a script
        filename = sys.argv[1]
    elif defaultfilename is not None:
        filename = config.CODE_DIR + defaultfilename

    return filename


class Inifile(dict):

    """
    Wrapper around an 'ini' file.

    This class merely reads a file containing key/value pairs and returns
    a dictionary:

    >>> #from . import config
    >>> import config
    >>> inifilename = config.CODE_DIR + 'adaqcode/inifile.ini'
    >>> ini_dict = Inifile(inifilename) # doctest: +ELLIPSIS
    Reading inifile .../inifile.ini

    Print some example types of data held in the dictionary:

    Numbers:

    >>> print(ini_dict['forecast_day'])
    1
    >>> isinstance(ini_dict['forecast_day'], int)
    True

    String:

    >>> print(ini_dict['obs_fmt'])
    aurn

    Datetime-format:

    >>> print(ini_dict['end_datetime'])
    2014-03-29 00:00:00

    Logical:

    >>> print(ini_dict['calc_stats'])
    True

    List:

    >>> print(ini_dict['line_colours_list'])
    ['k', 'green', 'orange', 'purple']

    Dictionary:

    >>> print(ini_dict['field_attribute_dict'])
    {'Species': 'PM10', 'Quantity': 'Air concentration'}

    .. note::

      * Blank lines in the file are ignored, as are lines beginning with hash.
      * If input file name ends with '.nl', file is read as a fortran namelist.
      * Any lines that end with '\' are treated as continuing on the following
        line.
      * Environment variables are expanded for values.
      * Any variables defined in config.py are converted to the value of the
        variable
        Eg. SAMPLE_DATADIR/oper is converted to
        /data/users/apdg/python_sample_data//oper
      * Any keys which end in "_list" are automatically converted into a list,
        even if they only contain 1 item.
      * If a list contains more than one item, it must be setup in the input
        file as a list,
        eg my_list = ["item1","item2"]
      * Any keys which end in "_dict" are automatically converted into a
        dictionary.
      * All variables are converted to python True/False/None if possible.
      * All strings in input file must be quoted, eg my_string = "hello"

    """

    #List of all known keys.
    #(All keys must be in lowercase, as required by Rose.)
    KEYS = ['range_days', 'start_date', 'end_date', 'forecast_day',
            'short_name_list', 'models_list', 'site_types_list',
            'models_dir_list',
            'models_fmt_list',
            'nimrod_data_type_list',
            'obs_dir', 'obs_dir_list',
            'obs_fmt', 'obsdir_fixed_csv',
            'obsdir_fixed_csv',
            'plot_dir',
            'units',
            'chem_units', 'aerosol_units', 'pollen_units',
            'sites_file', 'regimes_file', 'calc_stats',
            'histograms', 'soccer_plots', 'qqplots', 'timeseries',
            'diurnal', 'timeseries_of_stats', 'timeseries_of_stats_list',
            'timeseries_of_aggregate_stats',
            'timeseries_multiple_short_names',
            'timeseries_multiple_short_names_dict',
            'regime_plots', 'regime_stats_scatter',
            'contours', 'strict_statistics', 'html_dir', 'daqi',
            'line_colours_list', 'calc_stats_format_list',
            'rolling_stats_file_list', 'monthly_stats',
            'daqi_site_maps', 'mobrand', 'figsize', 'units',
            'scaling_factors_list',
            'save_cubes_nc',
            #For Difference plotting
            'plot_dir_case1', 'plot_dir_case2', 'plot_dir_diff',
            'plot_dir_reldiff', 'plot_dir_montage',
            'abs_diff_levels_list', 'rel_diff_levels_list',
            #For NAME
            'field_attribute_dict', 'z_level_list', 'z_leveltype', 'back',
            'percent',
            #For NAME ensemble
            'num_members', 'percentile_list',
            #For RSMC
            'source_note', 'arrival_time',
            #For field plotting...
            'extent_list', 'levels_list', 'levels_dict', 'nlevels', 'cmap',
            'cbar_orientation', 'cbar_label', 'cbar_num_fmt', 'cbar',
            'title', 'titleprefix',
            'suptitle', 'annote_location', 'annote', 'mapping', 'projection',
            'plottype', 'colorscale', 'marker_size',
            #For ensemble field plotting...
            'threshold_list',
            #For trajectory plotting
            'marker_interval', 'mo_branding', 'plotname',
            'release_info_list',
            #For cross-section plotting
            'cross_section', 'waypoints_list', 'max_height',
            #For mass retrievals...
            'models_op_list', 'models_ps_list', 'mass_filenames',
            'models_id_list', 'models_runtime_list', 'models_mass_dir_list',
            'models_moose_id_list', 'mass_retry_attempts', 'mass_retry_delay',
            #For difference plots
            'suptitle_case1', 'suptitle_case2', 'fail_if_diff',
            'models_dir_case1_list', 'models_dir_case2_list', 'reldiff_tolerance',
            #For mmr2dobsonunits
            'plot_single_time', 'p_target_top', 'mass_retrieve',
            'sat_pass_time_list']

    def __init__(self, filename=None, defaultfilename=None):
        """
        Initialise inifile class

        :param filename: String giving filename of ini file.
        :param defaultfilename: String selecting which inifile is used by
                                default if filename=None.

        Method:
          * Calculate filename if not set
          * Read file
          * Calculate extra variables
          * Convert some formats
        """

        filefmt = 'ini'

        if filename is None:
            filename = get_filename(defaultfilename)

        if filename is None:
            #No filename so cannot automatically read file.close
            return

        print('Reading inifile', filename)

        # Read ini file
        self.__readfile(filename, filefmt)

        # Check against known keys
        self.check_keys()

        # Format conversion if required
        self.check_formats()

        # Calculate dates and convert to datetime format
        self.calculate_dates()


    def __readfile(self, filename, filefmt):
        """
        Reads inifile and converts into a dictionary
        """

        if not os.path.exists(filename):
            raise ValueError("Specified filename %s not found" % filename)

        #Check format of file:
        # if ends in '.nl', then 'nl' (namelist) format
        # otherwise ini format

        if filename.endswith('.nl'):
            filefmt = 'nl'

        fin = open(filename, "r")
        data = fin.readlines()
        fin.close()

        line_continued = False
        previous_line = ''

        for line in data:

            # Skip blank lines or any beginning with '#'...
            if line_continued:
                line = previous_line+line

            if len(line) <= 1:
                continue
            if line[0] == '#':
                continue
            if filefmt == 'nl':
                if line[0] == '&':
                    continue
                if line[0] == '/':
                    continue

            #If line ends in '\' then treat as a line continuation with line
            #continued on following line. (Note extra \ needed to escape meaning
            #of \).
            if line.strip()[-1] == "\\":
                line_continued = True
                previous_line = line.strip()[:-1]
                continue
            else:
                line_continued = False

            #If namelist file, strip off ',' at end of line
            if filefmt == 'nl':
                line = line.strip() #first remove any trailing spaces
                line = line.rstrip(',')

            # Each param in the file becomes a dictionary key...
            arr = line.split("=", 1)
            if len(arr) != 2:
                continue
            keyword = arr[0].strip()

            #Expand any environment variables
            value = os.path.expandvars(arr[1].strip())

            #Replace any config values
            for k, v in config.__dict__.items():
                if isinstance(v, str) and k.isupper():
                    if value.find(k) != -1:
                        value = value.replace(k, v)

            #Check that values are readable Python literals and convert
            try:
                value = ast.literal_eval(value)
            except:
                warnings.warn("Inifile value for {} cannot be evaluated".format(keyword))
                continue

            #Check for quoted non-string variables
            try:
                value = ast.literal_eval(value)
            except:
                pass

            self[keyword] = value

    def check_keys(self):
        """
        Check keys passed in against known list of keys -
        raise a UserWarning if not a known key
        (May be using a new key for a temporary reason)
        """

        for key in self:
            if key not in self.KEYS:
                warnings.warn("Inifile - key not found: {}".format(key))

    def check_formats(self):

        """
        Does some format conversion if required:

          * Converts any keys which end in _list to a list
          * Converts any keys which end in _dict to a dictionary
        """

        removed = []

        for key, value in self.items():

            if key.split('_')[-1] == 'list':
                #Check for falsey value
                if not value:
                    removed.append(key)
                    warnings.warn("{} is an empty list, it will not be added "
                                  "to the ini_dict".format(key))
                #Convert to a list if necessary (eg if just a string read in).
                if not isinstance(value, list):
                    value = list([value])

                self[key] = value

            elif key.split('_')[-1] == 'dict':
                #Check value is a dictionary.
                if not isinstance(self[key], dict):
                    removed.append(key)
                    warnings.warn("Dictionary expected for {}, {} will not "
                                  "be added to the ini_dict".format(key, key))

        for key in removed:
            self.pop(key)

    def calculate_dates(self):
        """
        Calculate unknown dates from known variables
        (range_days, start_date, end_date), also convert to date-time formats.
        If none of range_days, start_date or end_date are given a warning is
        raised and no date information is added to dictionary.
        Otherwise, requires at least two from range_days, start_date, end_date
        """

        if 'start_date' in self:
            self['start_date'] = str(self['start_date'])
            self['start_datetime'] = datetime.strptime(self['start_date'],
                                                       "%Y%m%d")
        if 'end_date' in self:
            self['end_date'] = str(self['end_date'])
            self['end_datetime'] = datetime.strptime(self['end_date'], "%Y%m%d")

        if 'range_days' in self:
            if 'start_date' in self:
                #Check for consistency - calculate expected end date first
                end_datetime = self['start_datetime'] \
                               + timedelta(days=self['range_days'])
                if 'end_date' in self:
                    if end_datetime != self['end_datetime']:
                        raise ValueError("end_date is not consistent with"
                                         "start_date+range_days")
                else:
                    #Set end_date, end_datetime
                    self['end_datetime'] = end_datetime
                    self['end_date'] = end_datetime.strftime("%Y%m%d")

            elif 'end_date' in self:
                #Calculate start_date:
                start_datetime = self['end_datetime'] \
                                 - timedelta(days=self['range_days'])
                self['start_datetime'] = start_datetime
                self['start_date'] = start_datetime.strftime("%Y%m%d")

            else:
                warnings.warn("start and end dates not set: Need 2 of"
                              "range_days, start_date and end_date")

        elif 'start_date' in self:
            if 'end_date' in self:
                self['range_days'] = (self['end_datetime'] \
                                      - self['start_datetime']).days
            else:
                warnings.warn("end date not set: Need 2 of "
                              "range_days, start_date and end_date")

        elif 'end_date' in self:
            warnings.warn("start date not set: Need 2 of "
                          "range_days, start_date and end_date")


def get_inidict(inifilename=None, defaultfilename=None):
    """
    Simple wrapper to return only ini_dict dictionary

    :param inifilename: String giving filename of ini file.
    :param defaultfilename: String selecting which inifile is used by
                            default if inifilename=None.
    :return: ini_dict - Dictionary of a :class:`inifile` object.

    >>> ini_dict = get_inidict(defaultfilename='adaqcode/inifile.ini')
    ... # doctest: +ELLIPSIS
    Reading inifile .../adaqcode/inifile.ini
    """

    ini_dict = Inifile(inifilename, defaultfilename=defaultfilename)

    return ini_dict


if __name__ == '__main__':

    import doctest
    doctest.testmod()
