"""
sites_info.py - Contains SitesInfo class
"""
from __future__ import print_function

from six.moves.builtins import str
from six.moves.builtins import object

import os
import glob
import warnings
import numpy as np

import iris.cube
import matplotlib.pyplot as plt
import cartopy.feature
import cartopy.crs as ccrs

def format_string(data):
    """
    Generate a format string for printing a single row of an ndarray

    :returns: a format string

    Set up an example numpy record array:

    >>> data = np.array([(1.2345, 2, 'ab'), (20.125, 4, 'cde')],
    ... dtype=[('x', float), ('y', int), ('z', str, 4)])

    Generate the format string for this data:

    >>> fmtstr = format_string(data)
    >>> print(fmtstr)
    {:8.3f},{:2d},{:>3}

    Use this format string to print each row of the data in turn.
    Note the \*row syntax is used to pass each element of the row into each
    element of the format string.

    >>> for row in data:
    ...     print(fmtstr.format(*row))
       1.234, 2, ab
      20.125, 4,cde

    """

    if not isinstance(data, np.ndarray):
        raise ValueError('Input data is not a numpy ndarray')

    #Get headers
    headers = list(data.dtype.names)

    #Get format string for each column
    format_str = ""
    for header in headers:
        dtype_str = data.dtype.fields[header][0].str[1:]

        if header == 'site_id':
            #Need to keep full precision
            format_str += '{:17.8f},'
        elif (dtype_str[0] == 'S') or (dtype_str[0] == 'U'):
            #Right aligned string of length of maximum string in data
            maxlen = str(len(max(data[header], key=len)))
            format_str += '{:>' + maxlen + '},'
        elif dtype_str[0] == 'f':
            maxval = np.max(np.abs(data[header]))
            if maxval < 1.:
                format_str += '{:9.6f},'
            elif maxval < 10.:
                format_str += '{:7.3f},'
            elif maxval < 100.:
                format_str += '{:8.3f},'
            elif maxval < 1000.:
                format_str += '{:9.3f},'
            elif maxval < 10000.:
                format_str += '{:10.1f},'
            else:
                format_str += '{:10f},'
        elif dtype_str[0] == 'i':
            maxval = np.max(np.abs(data[header]))
            if maxval < 10:
                format_str += '{:2d},'
            elif maxval < 100:
                format_str += '{:3d},'
            elif maxval < 1000:
                format_str += '{:4d},'
            else:
                format_str += '{:10d},'

    #Remove final comma
    format_str = format_str[:-1]

    return format_str


class SitesInfo(object):

    """
    Holds site information.
    The data attribute (once filled) from this class is an ndarray, which
    can be treated as a list of dictionaries. Here each dictionary contains
    information about a single site. The keys are taken from the column headings
    from a file (converted to lower-case and tidied where necessary), or from
    coordinate names.

    >>> si = SitesInfo()

    Can read data in from a file:

    >>> import config
    >>> sample_datadir=config.SAMPLE_DATADIR+'AURN_obs/'
    >>> sites_data = si.read_from_file(sample_datadir+'aq_sites_DEFRA.txt',
    ... allsites=False, obsdir=sample_datadir)
    Number of sites:  15
    >>> print(sites_data) # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    [('GB0002', 'AH', 52.50..., -3.03..., 370, 'RURAL', 'Aston_Hill')
     ('GB0015', 'BOT', 52.93..., -0.81..., 32, 'SUBURBAN', 'Bottesford')
     ('GB0025', 'BUSH', 55.86..., -3.20..., 180, 'RURAL', 'Bush_Estate')
     ('GB0035', 'ESK', 55.31..., -3.20..., 269, 'RURAL', 'Eskdalemuir')
     ('GB0041', 'GLAZ', 53.4..., -2.47..., 21, 'SUBURBAN', 'Glazebury')
     ('GB0046', 'HM', 54.33..., -0.80..., 267, 'RURAL', 'High_Muffles')
     ('GB0045', 'HAR', 51.57..., -1.32..., 137, 'RURAL', 'Harwell')
     ('GB0050', 'LB', 53.40..., -1.75..., 420, 'RURAL', 'Ladybower')
     ('GB0076', 'LH', 50.79..., 0.18..., 125, 'RURAL', 'Lullington_Heath')
     ('GB0075', 'LN', 54.43..., -7.89..., 130, 'REMOTE', 'Lough_Navar')
     ('GB0097', 'ROCH', 51.45..., 0.63..., 14, 'RURAL', 'Rochester')
     ('GB0104', 'SIB', 52.29..., 1.46..., 46, 'REMOTE', 'Sibton')
     ('GB0113', 'SV', 57.73..., -4.77..., 270, 'REMOTE', 'Strath_Vaich')
     ('GB0128', 'YW', 50.59..., -3.71..., 119, 'RURAL', 'Yarner_Wood')
     ('GB0123', 'WFEN', 52.29..., 0.29..., 5, 'RURAL', 'Wicken_Fen')]

    To allow nice printing of this data, use the method formatstr to generate
    a format string specific to this set of data:

    >>> fmtstr = format_string(sites_data)
    >>> print(fmtstr)
    {:>6},{:>4},{:8.3f},{:7.3f},{:4d},{:>8},{:>16}
    >>> for site in sites_data:
    ...    print(fmtstr.format(*site))
    GB0002,  AH,  52.504, -3.034, 370,   RURAL,      Aston_Hill
    GB0015, BOT,  52.930, -0.815,  32,SUBURBAN,      Bottesford
    GB0025,BUSH,  55.862, -3.206, 180,   RURAL,     Bush_Estate
    GB0035, ESK,  55.315, -3.206, 269,   RURAL,     Eskdalemuir
    GB0041,GLAZ,  53.460, -2.473,  21,SUBURBAN,       Glazebury
    GB0046,  HM,  54.334, -0.809, 267,   RURAL,    High_Muffles
    GB0045, HAR,  51.571, -1.327, 137,   RURAL,         Harwell
    GB0050,  LB,  53.403, -1.752, 420,   RURAL,       Ladybower
    GB0076,  LH,  50.795,  0.182, 125,   RURAL,Lullington_Heath
    GB0075,  LN,  54.439, -7.899, 130,  REMOTE,     Lough_Navar
    GB0097,ROCH,  51.455,  0.634,  14,   RURAL,       Rochester
    GB0104, SIB,  52.294,  1.464,  46,  REMOTE,          Sibton
    GB0113,  SV,  57.732, -4.776, 270,  REMOTE,    Strath_Vaich
    GB0128,  YW,  50.597, -3.716, 119,   RURAL,     Yarner_Wood
    GB0123,WFEN,  52.299,  0.291,   5,   RURAL,      Wicken_Fen

    To get key information - firstly print the list of column names:

    >>> print(sites_data.dtype.names)
    ('gems_code', 'abbrev', 'latitude', 'longitude', 'site_altitude', \
'site_type', 'site_name')

    To get list of all site abbreviations:

    >>> print(sites_data['abbrev']) # doctest: +NORMALIZE_WHITESPACE
     ['AH' 'BOT' 'BUSH' 'ESK' 'GLAZ' 'HM' 'HAR' 'LB' 'LH' 'LN' 'ROCH' 'SIB' \
'SV' 'YW' 'WFEN']

    To get information about first site:

    >>> print(sites_data[0]) # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    ('GB0002', 'AH', 52.50..., -3.03..., 370, 'RURAL', 'Aston_Hill')


    .. keywords:
       allsites - set to False to only include sites where there is an
                  observation file - must have obsdir set as well.
       site_types - only returns sites where 'type' is included in this list

    Alternatively, data can be extracted from a previously generated cube
    in a sites_cube_list:

    >>> si = SitesInfo()
    >>> sample_datadir = config.SAMPLE_DATADIR+'sites_cube_list/'
    >>> import adaq_data
    >>> md = adaq_data.ADAQData()
    >>> md_scl = md.load_ts(sample_datadir+'aqum_oper_1days.nc')
    >>> sites_data = si.read_from_sites_cube(md_scl[0])
    >>> fmtstr = format_string(sites_data)
    >>> for site in sites_data:
    ...    print(fmtstr.format(*site))
    35628361.14059750,  YW, -1.896,-0.772509,  50.597, -3.716, \
119,     Yarner_Wood,           RURAL,  210.025
    35665278.14588334,ACTH,  3.386,-0.476009,  55.883, -3.347, \
260,Auchencorth_Moss,           RURAL,  160.710
    35696583.14250361,  AH,  0.005,-0.325150,  52.504, -3.034, \
370,      Aston_Hill,           RURAL,  235.209
    35790611.14715750, ABD,  4.658, 0.220976,  57.158, -2.094, \
 20,        Aberdeen,URBAN_BACKGROUND,   38.889
    35867333.14157111, HAR, -0.923, 0.729340,  51.571, -1.327, \
137,         Harwell,           RURAL,  104.795


    Finally, the site data can be written back out to file:

    >>> outputfile = config.CODE_DIR+'adaqdocs/figures/sites_file.txt'
    >>> si.write_to_file(outputfile) # doctest: +ELLIPSIS
    Written to file .../adaqdocs/figures/sites_file.txt
    >>> with open(outputfile,"r") as fin:
    ...     for line in fin:
    ...         print(line) # doctest: +NORMALIZE_WHITESPACE
    site_id abbrev grid_latitude grid_longitude latitude longitude \
site_altitude site_name site_type surface_altitude
    <BLANKLINE>
      35628361.1405974999 YW     -1.895994455   -0.772508508   50.597499850   \
-3.716388941      119 Yarner_Wood      RURAL            210.024750
    <BLANKLINE>
      35665278.1458833367 ACTH    3.386127197   -0.476009092   55.883335110   \
-3.347222328      260 Auchencorth_Moss RURAL            160.710464
    <BLANKLINE>
      35696583.1425036117 AH      0.004814994   -0.325150323   52.503612520   \
-3.034166574      370 Aston_Hill       RURAL            235.209167
    <BLANKLINE>
      35790611.1471574977 ABD     4.658122523    0.220975519   57.157501220   \
-2.093888760       20 Aberdeen         URBAN_BACKGROUND 38.888828
    <BLANKLINE>
      35867333.1415711120 HAR    -0.922965487    0.729340354   51.571109770   \
-1.326666594      137 Harwell          RURAL            104.795067
    <BLANKLINE>

    """

    def __init__(self):

        self.data = None
        self.filename = None
        self.allsites = None
        self.obsdir = None
        self.site_types = None
        self.sites_cube = None
        self.fmtstr = None



    def __read_header(self, filename):
        """
        Read in file to generate required headers and to alter the dtype of
        the data to ensure strings are read in as unicode data instead of bytes.

        :param filename: filename to read
        :returns: (headers, dtypes) - Where headers is a list of required header
                  names with those from the file reduced to lower case and some
                  keys ones renamed for consistency (eg from 'lat' to
                  'latitude'). dtypes is a list of tuples - one tuple for each
                  header, where the first value is the header name and the
                  second value is the data type, where strings are modified to
                  be unicode (instead of bytes, which is the default for
                  python3).

        """

        data = np.genfromtxt(filename, dtype=None, names=True)

        headers = []
        dtypes = []

        for name in data.dtype.names:


            #Convert header names to lowercase and modify to required names
            #where obvious
            header = name.lower()
            if header == 'lat':
                header = 'latitude'
            elif header == 'lon':
                header = 'longitude'
            elif header == 'altit':
                header = 'site_altitude'
            elif header == 'name':
                header = 'site_name'
            elif header == 'type':
                header = 'site_type'
            headers.append(header)

            #Check if loading a string - if so, will need to convert to unicode
            dtype = str(data[name].dtype)
            if dtype[:2] == '|S':
                #Set up new dtype, which is a tuple for each header name,
                #(header, dtype). Replace any strings str of same length.
                #Note str is bytes for python2 and unicode for python3.
                dtypes.append((header, str, int(dtype[2:])))
            else:
                dtypes.append((header, dtype))

        return headers, dtypes

    def read_from_file(self, filename=None, allsites=False,
                       obsdir=None, site_types=None):
        """
        Read site information from file with headers

        :param filename: File containing list of sites to consider
        :param allsites: If set to True then uses all sites in sites file,
                         If set to False, then only includes sites where there
                         is an observation file for this site
        :param obsdir: Directory, or list of directories containing observations
        :param site_types: List of site types to consider - if a site does not
                           have this site type then it is not included in the
                           returned list of sites.

        :returns: numpy ndarray containing site information for required sites.

        """

        self.filename = filename
        self.allsites = allsites
        self.obsdir = obsdir
        self.site_types = site_types
        if isinstance(self.obsdir, str):
            #Convert to a list if only a single directory passed in
            self.obsdir = [self.obsdir]

        if not os.path.exists(self.filename):
            raise ValueError("Specified filename %s not found" % self.filename)

        #Read in file once initially to get headers and dtypes.
        # - Modify headers to expected names
        # - Modify dtypes to ensure unicode data is read in instead of bytes.
        headers, dtypes = self.__read_header(self.filename)
        # Then read data in again with required headers and types.
        self.data = np.genfromtxt(self.filename, dtype=dtypes, names=headers,
                                  skip_header=1)

        # If single site in sites file,
        # then ensure array shape is still (nsites,)
        if self.data.ndim == 0:
            self.data = np.array([self.data])

        #Check bare-minimum information:
        if 'latitude' not in self.data.dtype.names:
            warnings.warn('lat not included in names')
        if 'longitude' not in self.data.dtype.names:
            warnings.warn('lon not included in names')
        if 'abbrev' not in self.data.dtype.names:
            warnings.warn('abbrev not included in names')

        if not allsites and not self.obsdir:
            raise ValueError("Obs directories not set - required as allsites=False")

        indices = []
        for i, abbrev in enumerate(self.data['abbrev']):
            indice = i

            if not allsites:
                #Don't use site if obs file does not exist
                allfilenames = []
                for obsdir in self.obsdir:
                    filenames = glob.glob(obsdir+'/'+abbrev+'_*')
                    if filenames:
                        allfilenames.append(filenames)
                if not allfilenames:
                    indice = None

            if site_types:
                #Only keep sites which are of required site type
                if self.data['site_type'][i] not in site_types:
                    indice = None

            if indice is not None:
                indices.append(indice)

        self.data = self.data[indices]
        print('Number of sites: ', len(self.data))

        if not self.data.size:
            raise IOError('No sites found - obsdir='+str(self.obsdir)+
                          ' , sites filename='+filename)

        return self.data

    def read_from_sites_cube(self, sites_cube=None):
        """
        Get the sites data from a sites_cube_list cube.
        This will generate the required ndarray with all required information.
        Adds any values from coordinates on the site_id dimension.
        """

        if sites_cube is None:
            raise IOError("No sites_cube given")
        if isinstance(sites_cube, iris.cube.CubeList):
            #If a cubelist has been given, just take first cube
            sites_cube = sites_cube[0]
        if not sites_cube.coords('site_id'):
            raise IOError("No site_id coordinate in given cube")

        self.sites_cube = sites_cube

        #Find which coordinates to use
        site_dim = sites_cube.coord_dims('site_id')
        site_dim_coords = [str(coord.name())
                           for coord in sites_cube.coords(dimensions=site_dim)]

        #Get the dtype of the data of these coordinates
        dtype = []
        for coord_name in site_dim_coords:
            dtype.append((coord_name, sites_cube.coord(coord_name).dtype))

        #Put into ndarray
        self.data = np.zeros((len(sites_cube.coord('site_id').points),),
                             dtype=dtype)
        for coord_name in site_dim_coords:
            self.data[coord_name] = sites_cube.coord(coord_name).points

        return self.data


    def write_to_file(self, outputfile):
        """
        Write information from sites cube out to a sites file.
        Will write out in the same (or very similar) format as the file
        that is read in.

        :param outputfile: File to write out to.
        """

        if self.data is None:
            raise IOError("Object does not have any data")

        #Get headers
        headers = list(self.data.dtype.names)
        header_str = len(headers)*"%s " % tuple(headers)

        #Get format string for each column
        format_str = ""
        for header in headers:
            dtype_str = self.data.dtype.fields[header][0].str[1:]

            if header in ['latitude', 'longitude', 'grid_latitude', 'grid_longitude']:
                format_str += '%14.9f'
            elif header == 'site_id':
                format_str += '%21.10f'
            elif (dtype_str[0] == 'S') or (dtype_str[0] == 'U'):
                format_str += '%-'+dtype_str[1:]+'s'
            elif dtype_str[0] == 'f':
                format_str += '%.'+str(2+int(dtype_str[1:]))+'f'
            elif dtype_str[0] == 'i':
                format_str += '%'+dtype_str[1:]+'d'

            format_str += ' '
        format_str += '\n'

        #Now output to file
        with open(outputfile, 'w') as fout:

            fout.write(header_str + "\n")

            for site in self.data:
                fout.write(format_str % tuple(site))

        print("Written to file " + outputfile)


    def plot_location_map(self, filename, label=False, title_prefix=None,
                          extent=None, labelsize='medium',
                          colours_attribute=None,
                          colours_dict=None):
        """
        Plot site locations on a map.

        :param filename: Output filename
        :param label: If set True, then annotates each location with its abbrev.
        :param extent: List of [lonmin, lonmax, latmin, latmax] to set
                       extent of map
        :param labelsize: Size of font used for labels:
                          Either an relative value of 'xx-small', 'x-small'
                          'small', 'medium', 'large', 'x-large', 'xx-large'
                          or an absolute font size, e.g., 12
        :param colours_attribute: String, attribute from data to use to define
                                  different colours for sites, eg 'site_type'
        :param colours_dict: Dictionary, defining possible variables from
                             required colour attribute and the colour they
                             should map to.

        >>> import config
        >>> si = SitesInfo()
        >>> sample_datadir=config.SAMPLE_DATADIR+'AURN_obs/'
        >>> sites_data = si.read_from_file( sample_datadir+'aq_sites_DEFRA.txt',
        ... allsites=False, obsdir=sample_datadir )
        Number of sites:  15
        >>> filename = config.CODE_DIR + "/adaqdocs/figures/site_locations.png"
        >>> si.plot_location_map(filename,label=True,
        ... title_prefix='Site Location Map', extent=[-10.,2.5,49.5,61.])
        ... # doctest: +ELLIPSIS
        Saved to file  .../adaqdocs/figures/site_locations.png

        .. image:: ../adaqdocs/figures/site_locations.png
           :scale: 50%

        Example of setting colours based on a site attribute, in this case
        the site type:

        >>> colours_attribute = 'site_type'
        >>> colours_dict = {'REMOTE': 'green',
        ... 'RURAL': 'orange', 'SUBURBAN': 'red'}
        >>> filename = config.CODE_DIR + "/adaqdocs/figures/site_locations_legend.png"
        >>> si.plot_location_map(filename,label=True,
        ... title_prefix='Site Location Map', extent=[-10.,2.5,49.5,61.],
        ... colours_attribute=colours_attribute, colours_dict=colours_dict)
        ... # doctest: +ELLIPSIS
        Saved to file  .../adaqdocs/figures/site_locations_legend.png

        .. image:: ../adaqdocs/figures/site_locations_legend.png
           :scale: 50%

        .. todo:: Once contouring capability is available, this routine may
                  want to be updated to use this capability if appropriate.

        """

        plt.figure()
        ax = plt.axes(projection=ccrs.PlateCarree())

        if colours_attribute in self.data.dtype.names:
            unique_att = sorted(set(self.data[colours_attribute]))
            if colours_dict is None:
                colours_dict = {}
            for att in unique_att:
                indices = self.data[colours_attribute] == att
                if att in colours_dict:
                    colour = colours_dict[att]
                else:
                    colour = 'black'
                ax.scatter(self.data['longitude'][indices], self.data['latitude'][indices],
                           c=colour, transform=ccrs.PlateCarree(),
                           edgecolors="k",
                           label=att.replace('_', ' '))
                ax.legend(scatterpoints=1, fontsize=labelsize)
        else:
            ax.scatter(self.data['longitude'], self.data['latitude'], c='r',
                       transform=ccrs.PlateCarree(),
                       edgecolors="k")


        ax.coastlines(resolution='50m')
        ax.add_feature(cartopy.feature.BORDERS, linestyle=':')
        if title_prefix is None:
            title_prefix = 'Site locations'
        title = title_prefix + '\n(%d sites)' % (len(self.data['longitude']))
        ax.set_title(title)

        if label:
            for site in self.data:
                if 'abbrev' in site.dtype.names:
                    ax.annotate(site['abbrev'],
                                (site['longitude'], site['latitude']),
                                transform=ccrs.PlateCarree(),
                                fontsize=labelsize)
                else:
                    ax.annotate(site['site_name'],
                                (site['longitude'], site['latitude']),
                                transform=ccrs.PlateCarree(),
                                fontsize=labelsize)

        if extent is not None:
            ax.set_extent(extent)

        plt.tight_layout()
        plt.savefig(filename)
        print("Saved to file ", filename)
        plt.close()


def get_siteinfo(ini_dict, sites_cube=None):
    """
    Wrapper function, to get site information from sites file or sites cube,
    as defined in ini_dict.

    :param ini_dict: Dictionary of a :class:`inifile` object.
                     Could (optionally) contain:

                      * 'sites_types_list'
                      * 'sites_file'
                      * 'obs_dir'

    :return: numpy ndarray containing site information data
             from a :class:`sites_info.SitesInfo` object

    Gets data from either a sites_file, or a sites_cube: if neither
    of these are set then None is returned.

    >>> import inifile
    >>> ini_dict = inifile.get_inidict(
    ... defaultfilename='adaqcode/inifile.ini') # doctest: +ELLIPSIS
    Reading inifile .../adaqcode/inifile.ini
    >>> sites_data = get_siteinfo(ini_dict)
    Number of sites:  5

    >>> fmtstr = format_string(sites_data)
    >>> for site in sites_data:
    ...    print(fmtstr.format(*site))
    GB0001, ABD,  57.158, -2.094,  20,URBAN_BACKGROUND,        Aberdeen
    GB0003,ACTH,  55.883, -3.347, 260,           RURAL,Auchencorth_Moss
    GB0002,  AH,  52.504, -3.034, 370,           RURAL,      Aston_Hill
    GB0045, HAR,  51.571, -1.327, 137,           RURAL,         Harwell
    GB0128,  YW,  50.597, -3.716, 119,           RURAL,     Yarner_Wood

    """


    #If site_types_list not in ini_dict, then set to None (etc)
    site_types = ini_dict.get('site_types_list', None)
    sites_file = ini_dict.get('sites_file', None)
    #Try and get obs_dir_list from ini_dict - if not available,
    #use obs_dir,
    # set to '' and set allsites=True to get every site, not
    # just those where obs data are available.
    obs_dir_list = ini_dict.get('obs_dir_list', None)
    if not obs_dir_list:
        if 'obs_dir' in ini_dict:
            obs_dir_list = [ini_dict['obs_dir']]
    obs_fmt = ini_dict.get('obs_fmt', None)
    if obs_dir_list and obs_fmt == 'aurn':
        #Can only limit sites if obs filename is of format abbrev*,
        #which corresponds to aurn data only.
        allsites = False
    else:
        allsites = True

    sites = SitesInfo()

    if sites_file is not None:
        #Read from file
        sites_data = sites.read_from_file(sites_file,
                                          allsites=allsites,
                                          obsdir=obs_dir_list,
                                          site_types=site_types)
    elif sites_cube is not None:
        #Extract from sites_cube
        sites_data = sites.read_from_sites_cube(sites_cube=sites_cube)
    else:
        sites_data = None

    return sites_data


if __name__ == "__main__":

    import doctest
    doctest.testmod()
    #doctest.run_docstring_examples(get_siteinfo, globals())
