"""
Class for producing a html website,
and other functions to produce specific designs of websites,
eg for showing output from aq_plot.
"""
from __future__ import print_function

from six.moves.builtins import str
from six.moves.builtins import object

import os
import glob
import shutil
import warnings
import datetime

import numpy as np

import cube_statistics
import config
import timeseries_stats


class Website(object):
    """
    Class to build a simple html website page.

    >>> import config
    >>> html_dir = config.CODE_DIR+'/adaqdocs/figures/'

    The class should be initialised with the directory of where it will based.

    >>> web = Website(html_dir)

    A title and section title could then be added:

    >>> web.add_title('Example Webpage')
    >>> web.add_section_header('Example text file')

    An ascii text file could be added to appear in a frame (so if the frame
    is too small it can be scrolled around to see the rest of the file):

    >>> web.add_text_file('sites_file.txt', width=800, height=100)

    An image file could be added, giving it a fixed width:

    >>> web.add_section_header('Example single plot')
    >>> web.add_image('Harwell_O3.png', width=200)

    Alternatively, lots of images can be added together in a table.
    This is particularly useful if want for example each row to display
    a separate site and each column then to display a different species.
    In this simple example we will select fields from the
    figures directory under adaq_plotting/gridded_fields and have the first
    row to be figures that match *AQ* in the filename and the second row
    to be figures that match *NAME* in their filename.
    We will also set the width of the figures to be 30% of the webpage width.

    >>> web.add_section_header('Example plots table')
    >>> web.add_images_table('adaq_plotting/gridded_fields/',
    ... row_searchstr=['AQ','NAME'], width='30%')

    A footer and header will automatically be added to the website when it
    is initialised. However these could be overridden to add extra options.
    For example, the header could be given a title which will appear at the
    top of a browser window:

    >>> web.add_header('Test Webpage')

    For the footer, address details can be added, for example a telephone
    number to contact the page owner:

    >>> web.add_footer(telephone='xxxx')

    Lastly this webpage can be written to a html file. This will be
    generated in the directory set up on initialising the class.

    >>> web.write_file('test_webpage.html')  # doctest: +ELLIPSIS
    Saved website  .../adaqdocs/figures/test_webpage.html

    To see this webpage, look at adaqdocs/figures/test_webpage.html
    in a browser.
    To see the text behind this webpage:

    >>> with open(html_dir + 'test_webpage.html',"r") as fin:
    ...    for line in fin:
    ...        print(line)
    <html>
    <BLANKLINE>
    <head>
    <BLANKLINE>
    <title>Test Webpage</title>
    <BLANKLINE>
    </head>
    <BLANKLINE>
    <body>
    <BLANKLINE>
    <BLANKLINE>
    <BLANKLINE>
    <BLANKLINE>
    <BLANKLINE>
    <h1>Example Webpage</h1>
    <BLANKLINE>
    <BLANKLINE>
    <BLANKLINE>
    <BLANKLINE>
    <BLANKLINE>
    <h2>Example text file</h2>
    <BLANKLINE>
    <iframe width=800 height=100 src="sites_file.txt"></iframe>
    <BLANKLINE>
    <BLANKLINE>
    <BLANKLINE>
    <h2>Example single plot</h2>
    <BLANKLINE>
    <a href="Harwell_O3.png"><img width=200 src="Harwell_O3.png"> </a>
    <BLANKLINE>
    <BLANKLINE>
    <BLANKLINE>
    <h2>Example plots table</h2>
    <BLANKLINE>
    <BLANKLINE>
    <BLANKLINE>
    <table border="0" cellspacing="0" cellpadding="0">
    <BLANKLINE>
    <tr>
    <BLANKLINE>
    <td><a href="adaq_plotting/gridded_fields//Fieldplot_aqum_oper_O3_\
201404020000_basicAQ.png"><img width=30% src="adaq_plotting/gridded_fields//\
Fieldplot_aqum_oper_O3_201404020000_basicAQ.png"> </a></td>
    <BLANKLINE>
    <td><a href="adaq_plotting/gridded_fields//Fieldplot_aqum_oper_O3_\
201404020000_fullAQ.png"><img width=30% src="adaq_plotting/gridded_fields//\
Fieldplot_aqum_oper_O3_201404020000_fullAQ.png"> </a></td>
    <BLANKLINE>
    </tr>
    <BLANKLINE>
    <tr>
    <BLANKLINE>
    <td><a href="adaq_plotting/gridded_fields//Fieldplot_aqum_oper_O3_\
201404020000_basicNAME.png"><img width=30% src="adaq_plotting/gridded_fields//\
Fieldplot_aqum_oper_O3_201404020000_basicNAME.png"> </a></td>
    <BLANKLINE>
    </tr>
    <BLANKLINE>
    </table>
    <BLANKLINE>
    <BLANKLINE>
    <BLANKLINE>
    <p><hr><p>
    <BLANKLINE>
    <br><I>Modified:</I> <B><!--#config timefmt="%d %b %Y"-->\
<!--#echo var="LAST_MODIFIED"--> </B>
    <BLANKLINE>
    <p><hr><p>
    <BLANKLINE>
    <address>
    <BLANKLINE>
    Telephone: xxxx<br></address>
    <BLANKLINE>
    </body>
    <BLANKLINE>
    </html>

    """


    def __init__(self, html_dir='./'):

        self.header = None
        self.footer = None
        self.html_dir = html_dir
        self.html_file = ''
        self.body = ''

        #Make html directory
        if not os.path.exists(self.html_dir):
            os.makedirs(self.html_dir)

        #Ensure html directory has / on end
        if self.html_dir[-1] != '/':
            self.html_dir += '/'

        #Add header and footer so not forgotten by user,
        #although can be overwritten by calling these
        #routines later.
        self.add_header()
        self.add_footer()

    def add_header(self, pagename=None):
        """
        Add header to website.
        This must be included at the start of any webpage.

        :param pagename: Name to give to page (appears in tab/browser header).

        ..Note:: This can be overwritten by setting self.header, but should
                 still include '<html>\\n<body>'

        """

        self.header = '<html>\n<head>\n'
        if pagename is not None:
            self.header += '<title>'+pagename+'</title>\n'
        self.header += '</head>\n<body>\n'

    def add_footer(self, operational_warning=False,
                   name=None,
                   desk_number=None,
                   telephone=None,
                   email=None):
        """
        Add footer to website.
        This must be included at the end of any webpage.

        :param operational_warning: Display warning in red, "WARNING: This page
                                    is not supported operationally."
        :param name: Real name to add to address section and to
                     who modified last.
        :param desk_number: Desk number to add to address section
        :param telephone: Telephone number to add to address section
        :param email: Email address to add to address section.


        ..Note:: This can be overwritten by setting self.footer, but should
                 still end with '</body>\\n</html>'

        """

        footer = '\n<p><hr><p>\n'

        if operational_warning:
            footer += ('<font color="#990000"> WARNING: This page is not '
                       'supported operationally.</font>\n')

        footer += ('<br><I>Modified:</I> <B><!--#config timefmt="%d %b %Y"'
                   '--><!--#echo var="LAST_MODIFIED"--> </B>')
        if name is not None:
            footer += '<i>by</i><b> '+name+'</b>'
        footer += '\n<p><hr><p>\n'

        footer += '<address>\n'
        if name is not None:
            footer += name + '<br>'
        if desk_number is not None:
            footer += 'Location: ' + desk_number + '<br>'
        if telephone is not None:
            footer += 'Telephone: ' + str(telephone) + '<br>'
        if email is not None:
            footer += 'Email: <a href="mailto:' + email + '">'
            footer += email + '</a>'
        footer += '</address>\n'

        footer += '</body>\n</html>'

        self.footer = footer

    def add_title(self, title):
        """
        Add a title to the page (uses <h1>)
        """

        self.body += '\n\n<h1>'+title+'</h1>\n'

    def add_section_header(self, header):
        """
        Add a section header to the page (uses <h2>)
        """

        self.body += '\n\n<h2>'+header+'</h2>\n'

    def add_symlink(self, input_dir, sym_dir):
        """
        Set up a symbolic link from input_dir to sym_dir.
        This can be useful to provide links to places not
        under public_html etc. Note sym_dir is a
        relative directory compared to self.html_dir.
        """

        sym_dir = self.html_dir + sym_dir
        #First remove symbolic link if it already exists
        if os.path.islink(sym_dir):
            os.unlink(sym_dir)
        #Then add the symbolic link
        os.symlink(input_dir, sym_dir)

    def add_image(self, plotfilename, width=None):
        """
        Add an image - clicking on it will open it up larger

        :param plotfilename: Filename of plot to display.
                             Note this should be as a web address, or
                             as a relative link to self.html_dir
        :param width: Width of image on page. Can be in pixels (eg 500),
                      or as percentage of the page width (eg 100%).
        """

        self.body += '<a href="'+plotfilename+'"><img '
        if width is not None:
            self.body += 'width='+str(width)+' '
        self.body += 'src="'+plotfilename+'"> </a>'

    def add_images_table(self, sym_dir='./', row_searchstr=None,
                         col_searchstr=None, width=None, wildcards=True,
                         reverse_row_cols=False, headers=None,
                         ignore_string=None):
        """
        Add a selection of images in a table.

        :param sym_dir: Link to plots directory relative to webpage
                        (generally as a symbolic link)
        :param row_searchstr: list of strings to search filenames for to put on
                              rows
        :param col_searchstr: list of strings to search filenames for to put on
                              columns. Note, the col_searchstr should always be
                              after the row_searchstr in filenames
        :param width: Width of each plot on page. Can be in pixels (eg 500),
                      or as percentage of the page width (eg 100%).
        :param wildcards: Add wildcards (*) in between directory, row_searchstr
                          and col_searchstr. If false, then these are all
                          concatenated together.
        :param reverse_row_cols: If True, then row_searchstr should instead be
                                 after the col_searchstr in the filename search.
        :param headers: List of column headers.
        :param ignore_string: String which if included in filename indicates this
                              file should not be used.

        This does a loop over all the images in first rows whose filenames match
        the required row_searchstr and columns whose filenames match
        col_searchstr. Note if no filenames are matched, a table is still built
        (and appears within the html), but will not contain any images, so won't
        appear on the displayed webpage.

        """

        #Set search strings to list containing an empty string.
        #This is then the equivalent of a search for '*'
        if row_searchstr is None:
            row_searchstr = ['']
        if col_searchstr is None:
            col_searchstr = ['']

        self.body += '\n<table border="0" cellspacing="0" cellpadding="0">\n'

        if headers is not None:
            self.body += '<tr>\n'
            for header in headers:
                self.body += '<td align="center"><b>' + header + '</b></td>\n'
            self.body += '</tr>\n'

        for row_name in row_searchstr:
            self.body += '<tr>\n'
            for col_name in col_searchstr:
                if wildcards:
                    search_string = self.html_dir + '/' + sym_dir + '/*' + \
                                    row_name + '*' + col_name + '*'
                    if reverse_row_cols:
                        search_string = self.html_dir + '/' + sym_dir + '/*' + \
                                        col_name + '*' + row_name + '*'

                else:
                    search_string = self.html_dir + '/' + sym_dir + '/' + \
                                    row_name + col_name
                    if reverse_row_cols:
                        search_string = self.html_dir + '/' + sym_dir + '/' + \
                                        col_name + row_name

                filenames = [os.path.basename(filename)
                             for filename in sorted(glob.glob(search_string))]

                for filename in filenames:
                    #Only consider files - ignore directories
                    if os.path.isfile(self.html_dir+'/'+sym_dir+'/'+filename):
                        if ignore_string is not None:
                            if ignore_string in filename:
                                continue
                        self.body += '<td>'
                        plotfilename = sym_dir + '/' + filename
                        self.add_image(plotfilename, width=width)
                        self.body += '</td>\n'
            self.body += '</tr>\n'

        self.body += '</table>\n'

    def add_text(self, text):
        """
        Add any text to webpage. Can include html formatting.

        :param text: Text to add to page.

        """

        self.body += text


    def add_text_file(self, filename, width='100%', height='50%'):
        """
        Include a text file in webpage in a frame (allows scrolling).

        :param filename: Filename of text file to be added
        :param width: Width of frame on page. Can be in pixels (eg 500),
                      or as percentage of the page width (eg 100%).
        :param height: Height of frame on page. Can be in pixels (eg 500),
                      or as percentage of the page height (eg 100%).

        """

        self.body += '<iframe width='+str(width)+' height='+str(height)+' '
        self.body += 'src="'+filename+'"></iframe>'


    def write_file(self, filename):
        """
        Write out to html file.

        :param filename: File to output to, given relative to self.html_dir

        """

        self.html_file = self.html_dir + filename

        with open(self.html_file, 'w') as fout:
            fout.write(self.header)
            fout.write(self.body)
            fout.write(self.footer)

        print('Saved website ', self.html_file)


def write_aq_plot_tseries_website(ini_dict, sites_data, site_types=None,
                                  short_names=None, filename=None):
    """
    Writes a website to display timeseries and location of sites.

    :param ini_dict: Dictionary of a :class:`inifile` object. Should contain:

                     * 'plot_dir' - directory containing plots
                     * 'html_dir' - directory to create website in.

    :param sites_data: numpy ndarray containing site information data
                       from a :class:`sites_info.SitesInfo` object.
                       Note this is required to extract the names of sites,
                       plus the types of sites.
    :param site_types: List of different site types to display. Other site
                       types will not be shown.
                       If set to None, then displays all sites.
    :param short_names: List of short_names to include. If set to None, then
                        uses ini_dict['short_name_list']
    :param filename: Filename of output website. Will be placed in html_dir.

    Note a symbolic link ('plots/') is set up in html_dir to point to plot_dir

    >>> import adaq_functions
    >>> ini_dict, sites_data, od, md_list = adaq_functions.get_exampledata()
    ... # doctest: +ELLIPSIS
    Reading inifile .../adaqcode/example_data_1days.ini
    Number of sites:  5

    For this example only, manually change plot and html directory to
    save to gallery,
    and ensure that timeseries is set to True in the ini_dict:

    >>> import config
    >>> ini_dict['plot_dir'] = "./adaq_plotting"
    >>> ini_dict['html_dir'] = config.CODE_DIR + \
    "/adaqdocs/figures/"
    >>> ini_dict['timeseries'] = True

    Now can write out to html website:

    >>> write_aq_plot_tseries_website(ini_dict, sites_data) # doctest: +ELLIPSIS
    Saved website  .../adaqdocs/figures/aq_plot_tseries.html

    Note this webpage be viewed by opening adaqdocs/figures/aq_plot_tseries.html
    in a browser.

    Finally check contents of this file:

    >>> with open(ini_dict['html_dir'] + 'aq_plot_tseries.html', "r") as fin:
    ...    for line in fin:
    ...        print(line.strip())  # doctest: +ELLIPSIS
    <html>
    <head>
    <title>AQ plot timeseries</title>
    </head>
    <body>
    <BLANKLINE>
    <BLANKLINE>
    <h1>AQ plot timeseries for 20140402 - 20140403</h1>
    <BLANKLINE>
    <BLANKLINE>
    <h2>Timeseries</h2>
    <BLANKLINE>
    <table border="0" cellspacing="0" cellpadding="0">
    <tr>
    <td><a href="plots/timeseries//Aberdeen_O3.png"><img width=500 \
src="plots/timeseries//Aberdeen_O3.png"> </a></td>
    </tr>
    <tr>
    </tr>
    <tr>
    </tr>
    <tr>
    </tr>
    <tr>
    </tr>
    </table>
    <BLANKLINE>
    <p><hr><p>
    <font color="#990000"> WARNING: This page is not supported \
operationally.</font>
    <br><I>Modified:</I> <B><!--#config timefmt="%d %b %Y"--><!--#echo \
var="LAST_MODIFIED"--> </B>
    <p><hr><p>
    <address>
    </address>
    </body>
    </html>

    """


    #Set up symbolic link to plot directory
    input_dir = ini_dict['plot_dir']
    html_dir = ini_dict.get('html_dir', './')
    if not os.path.exists(html_dir):
        os.makedirs(html_dir)
    sym_dir = 'plots'

    web = Website(html_dir)
    web.add_symlink(input_dir, sym_dir)

    #Start html page
    web.add_header('AQ plot timeseries')
    title = 'AQ plot timeseries for ' + ini_dict['start_date'] + \
            ' - ' + ini_dict['end_date']
    web.add_title(title)

    #If available, add site locations map
    if os.path.isfile(html_dir + '/' + sym_dir + '/site_locations.png'):
        web.add_section_header('Site Locations')
        web.add_image(sym_dir +'/site_locations.png')

    #Add timeseries plots
    if ini_dict.get('timeseries', False) or \
       ini_dict.get('timeseries_multiple_short_names', False):
        #Pick out sites - eg only display rural sites
        if site_types is None:
            site_names = sorted(sites_data['site_name'])
        else:
            site_names = []
            for site_type in site_types:
                indices = np.where(sites_data['site_type'] == site_type)
                if indices[0].size:
                    site_names += list(sites_data['site_name'][indices[0]])
            site_names = sorted(site_names)
        #All filenames have _ after site name
        site_names = [sn + '_' for sn in site_names]
        if not site_names:
            warnings.warn('No timeseries sites for ' + ' '.join(site_types)
                          +' - webpage not written')
            #No sites to display, so don't write a webpage
            return

        web.add_section_header('Timeseries')
        #Each row is for a single site, ordered in terms of short_name_list
        if short_names is None:
            short_names = ini_dict['short_name_list']
        col_names = [name+'.png' for name in short_names]
        web.add_images_table(sym_dir+'/timeseries/', site_names, col_names,
                             width=500, wildcards=False)

    #Add footer
    web.add_footer(operational_warning=True)

    #Output to file
    if filename is None:
        filename = 'aq_plot_tseries'
        if site_types is not None:
            filename += '_'.join(site_types)
        filename += '.html'
    web.write_file(filename)



def write_aq_plot_stats_website(ini_dict):
    """
    Writes a website to display statistics and statisitical plots generated
    by aq_plot.

    :param ini_dict: Dictionary of a :class:`inifile` object. Should contain:

                     * 'plot_dir' - directory containing plots
                     * 'html_dir' - directory to create website in.
                     * 'calc_stats' - if True, then displays statistics file

    Note a symbolic link ('plots/') is set up in html_dir to point to plot_dir

    >>> import adaq_functions
    >>> ini_dict, sites_data, od, md_list = adaq_functions.get_exampledata()
    ... # doctest: +ELLIPSIS
    Reading inifile .../adaqcode/example_data_1days.ini
    Number of sites:  5

    For this example only, manually change plot and html directory to
    save to gallery.

    >>> import config
    >>> ini_dict['plot_dir'] = "./adaq_plotting"
    >>> ini_dict['html_dir'] = config.CODE_DIR + \
    "/adaqdocs/figures/"

    Now can write out to html website:

    >>> write_aq_plot_stats_website(ini_dict) # doctest: +ELLIPSIS
    Saved website  .../adaqdocs/figures/aq_plot_stats.html

    Note this webpage be viewed by opening adaqdocs/figures/aq_plot_stats.html
    in a browser.

    Finally check contents of this file:

    >>> with open(ini_dict['html_dir'] + 'aq_plot_stats.html', "r") as fin:
    ...    for line in fin:
    ...        print(line.strip())  # doctest: +ELLIPSIS
    <html>
    <head>
    <title>AQ plot statistics</title>
    </head>
    <body>
    <BLANKLINE>
    <BLANKLINE>
    <h1>AQ plot statistics for 20140402 - 20140403</h1>
    <BLANKLINE>
    <BLANKLINE>
    <h2>Statistics Plots</h2>
    <BLANKLINE>
    <table border="0" cellspacing="0" cellpadding="0">
    <tr>
    <td><a href="plots/Histogram_O3.png"><img width=500 \
src="plots/Histogram_O3.png"> </a></td>
    </tr>
    <tr>
    <td><a href="plots/Soccer_Plot_O3.png"><img width=500 \
src="plots/Soccer_Plot_O3.png"> </a></td>
    </tr>
    <tr>
    <td><a href="plots/Quantile-Quantile_O3.png"><img width=500 \
src="plots/Quantile-Quantile_O3.png"> </a></td>
    </tr>
    <tr>
    <td><a href="plots/Diurnal_O3.png"><img width=500 \
src="plots/Diurnal_O3.png"> </a></td>
    </tr>
    ...
    <tr>
    </tr>
    <tr>
    <td><a href="plots/Timeseries_of_bias_O3.png">\
<img width=500 src="plots/Timeseries_of_bias_O3.png"> </a></td>
    </tr>
    ...
    <tr>
    <td><a href="plots/Timeseries_of_rmse_O3.png">\
<img width=500 src="plots/Timeseries_of_rmse_O3.png"> </a></td>
    </tr>
    ...
    <p><hr><p>
    <font color="#990000"> WARNING: This page is not supported \
operationally.</font>
    <br><I>Modified:</I> <B><!--#config timefmt="%d %b %Y"-->\
<!--#echo var="LAST_MODIFIED"--> </B>
    <p><hr><p>
    <address>
    </address>
    </body>
    </html>
    """



    #Set up symbolic link to plot directory
    input_dir = ini_dict['plot_dir']
    html_dir = ini_dict.get('html_dir', './')
    if not os.path.exists(html_dir):
        os.makedirs(html_dir)
    sym_dir = 'plots'

    web = Website(html_dir)
    web.add_symlink(input_dir, sym_dir)

    #Start html page
    web.add_header('AQ plot statistics')
    title = 'AQ plot statistics for ' + ini_dict['start_date'] + \
            ' - ' + ini_dict['end_date']
    web.add_title(title)

    #Add statistics csv file
    if ini_dict.get('calc_stats', False):
        #csv file does not display nicely, so add symbolic link
        #to rename to a text file
        web.add_symlink(html_dir+'/'+sym_dir+'/stats.csv',
                        sym_dir+'/stats.txt')
        web.add_text_file(sym_dir+'/stats.txt', width=700, height=500)

    #Add statistics plots

    web.add_section_header('Statistics Plots')
    #Each row is for a single stat type, ordered in terms of short_name_list
    #Note if any/all of stats are missing, a table is created,
    #but nothing displays on webpage
    stats_filenames = ['Histogram_', 'Soccer_Plot_',
                       'Quantile-Quantile_', 'Diurnal_']
    for stat in cube_statistics.CUBE_AGGREGATORS:
        stats_filenames.append('Timeseries_of_sites_'+stat.lower()+'_')
    for stat in timeseries_stats.STATS_INFO:
        stats_filenames.append('Timeseries_of_'+stat+'_')

    col_names = [name+'.png' for name in ini_dict['short_name_list']]
    web.add_images_table(sym_dir, stats_filenames, col_names, width=500,
                         wildcards=False)

    #Add footer
    web.add_footer(operational_warning=True)
    web.write_file('aq_plot_stats.html')

def write_aq_plot_contours_website(ini_dict):
    """
    Produce a contours / gridded fields animated html browser for each model.

    :param ini_dict: Dictionary of a :class:`inifile` object. Should contain:

                     * 'plot_dir' - directory containing plots
                     * 'html_dir' - directory to create website in.
                     * 'short_name_list' - list of short_names to include on
                       website.
                     * 'models_list' - list of model names to produce websites
                       for.

    Note a symbolic link ('plots/') is set up in html_dir to point to plot_dir

    .. note:: Assumes files are of the format
              gridded_fields/Fieldplot_<model>_<short_name>_<yyyymmddHHMM>_*.png

    Example using images already in figures directory:

    >>> ini_dict = {
    ... 'plot_dir': config.CODE_DIR + "/adaqdocs/figures/adaq_plotting",
    ... 'html_dir': config.CODE_DIR + "/adaqdocs/figures/",
    ... 'short_name_list': ['O3'],
    ... 'models_list': ['aqum_oper']}

    Now can write out to html website:
    >>> write_aq_plot_contours_website(ini_dict) # doctest: +ELLIPSIS
    Saved website  /.../aq_plot_gridded_fields_aqum_oper.html

    Note this webpage be viewed by opening
    adaqdocs/figures/aq_plot_gridded_fields_aqum_oper.html in a browser.

    Finally check contents of this file:

    >>> webfile = ini_dict['html_dir'] + 'aq_plot_gridded_fields_aqum_oper.html'
    >>> with open(webfile, "r") as fin:
    ...    for line in fin:
    ...        print(line.strip())  # doctest: +ELLIPSIS
    <html>
    <head>
    <title>aqum_oper- gridded fields</title>
    </head>
    <body>
    <BLANKLINE>
    <BLANKLINE>
    <h1>Gridded Fields for aqum_oper</h1>
    <BLANKLINE>
    <SCRIPT LANGUAGE="JavaScript">
    loadicons();
    var FilenamesDict = {'O3': \
['./plots/gridded_fields/Fieldplot_aqum_oper_O3_201404020000_basicAQ.png', \
'./plots/gridded_fields/Fieldplot_aqum_oper_O3_201404020000_basicNAME.png', \
'./plots/gridded_fields/Fieldplot_aqum_oper_O3_201404020000_fullAQ.png']}
    mydiag   = "O3";
    tag = ['O3'];
    <BLANKLINE>
    //**************************************************************************
    //First part of include text for contour browser
    //Contains main part of javascript code.
    ...
    <!--  End of first part for use in contour browser
    **************************************************************
    -->
    <BLANKLINE>
    <-- List of options for drop-down selection-->
    <OPTION value="0" SELECTED>O3</OPTION>
    <BLANKLINE>
    <!--
    ************************************************************
    Second part of fixed text for contour browser
    -->
    ...
    <!--
    End of Second part of fixed text for contour browser
    **************************************************************
    -->
    <BLANKLINE>
    <p><hr><p>
    <font color="#990000"> WARNING: This page is not supported operationally.\
</font>
    <br><I>Modified:</I> <B><!--#config timefmt="%d %b %Y"--><!--#echo \
var="LAST_MODIFIED"--> </B>
    <p><hr><p>
    <address>
    </address>
    </body>
    </html>
    """

    #Set up symbolic link to plot directory
    input_dir = ini_dict['plot_dir']
    html_dir = ini_dict.get('html_dir', './')
    if not os.path.exists(html_dir):
        os.makedirs(html_dir)
    sym_dir = 'plots'


    #Copy required icon images to html directory
    icons_dir = html_dir+'/icons/'
    if not os.path.exists(icons_dir):
        os.makedirs(icons_dir)
    icon_list = ['back.jpg', 'faster.jpg', 'forward.jpg', 'play.jpg',
                 'slower.jpg', 'stop.jpg',
                 'back2.jpg', 'faster2.jpg', 'forward2.jpg', 'play2.jpg',
                 'slower2.jpg', 'stop2.jpg', 'spacer.gif']
    for icon_file in icon_list:
        shutil.copy(config.CODE_DIR + '/adaqdocs/figures/html_source/'
                    + icon_file, icons_dir)

    #Loop over models creating a website for each model
    for model in ini_dict['models_list']:

        #Start html page
        web = Website(html_dir)
        web.add_symlink(input_dir, sym_dir)
        web.add_header(model + '- gridded fields')
        web.add_title('Gridded Fields for '+model)

        #--------------

        #Start writing contour browser using javascript
        js = '\n<SCRIPT LANGUAGE="JavaScript">' #Setup javascript string

        js += '\nloadicons();'

        #Generate dictionary of contour filenames
        filenames = {}
        for short_name in ini_dict['short_name_list']:
            #Generate list of files by pointing to real location
            filename_list = sorted(glob.glob(
                html_dir + '/' + sym_dir + '/gridded_fields/Fieldplot_' +
                model.replace(' ', '_') + '_' + short_name + '_*.png'))
            #Now modify to point to location relative to html_dir
            filenames[short_name] = []
            for fname in filename_list:
                if 'DAQI' not in short_name:
                    if short_name + '_DAQI' in fname:
                        continue #Don't include in list of filenames
                filenames[short_name].append(
                    './' + sym_dir + '/gridded_fields/' +
                    os.path.basename(fname))

        # Dictionary of filenames, labelled by short_name
        js += '\nvar FilenamesDict = '+str(filenames)
        # Initial short_name to display
        js += '\nmydiag   = "'+ ini_dict['short_name_list'][0] +'"; '
        # Full list of short_names to display
        js += '\ntag = '+str(ini_dict['short_name_list'])+';'
        js += '\n\n'


        web.body += js

        #Add in all the other javascript code and functions,
        #plus start main body of html, by reading from fixed text file

        with open(config.CODE_DIR +
                  '/adaqdocs/figures/html_source/contour_browser_part1.html',
                  'r') as fin:
            for line in fin:
                web.body += line

        #Add drop-down options list (started in above fixed text)

        web.body += '\n<-- List of options for drop-down selection-->'
        web.body += ('\n<OPTION value="0" SELECTED>'
                     + ini_dict['short_name_list'][0] + '</OPTION>')
        for i, short_name in enumerate(ini_dict['short_name_list'][1:]):
            web.body += ('\n<OPTION value="' + str(i+1) + '">'
                         + short_name + '</OPTION>')
        web.body += '\n\n'

        #Now finish off required text to actually call the animations and
        #to change short_name when requested.

        with open(config.CODE_DIR +
                  '/adaqdocs/figures/html_source/contour_browser_part2.html',
                  'r') as fin:
            for line in fin:
                web.body += line


        #--------------
        #Add footer
        web.add_footer(operational_warning=True)
        web.write_file('aq_plot_gridded_fields_'
                       + model.replace(' ', '_')+'.html')

def write_aq_plot_sites_website(ini_dict, daqi_only=True):
    """
    Produce a website for each short_name which displays the model data at
    observation sites on a map, alongside the observations and finally the
    contour plots for the same time/short_name.
    By default only plots DAQI and components of DAQI (eg O3_DAQI).
    If there are no models available (ie ini_dict['models_list'] is empty),
    then only observations are plotted.

    :param ini_dict: Dictionary of a :class:`inifile` object. Should contain:

                     * 'plot_dir' - directory containing plots
                     * 'html_dir' - directory to create website in.
                     * 'short_name_list' - list of short_names to produce
                        websites for
                     * 'models_list' - list of model names to include on
                        websites.
                     * 'start_datetime' - start of dates to display, in datetime
                        format.
                     * 'end_datetime' - end of dates to display, in datetime
                        format.
    :param daqi_only: If True, only produces webpages for DAQI and components of
                      DAQI. If False, produces webpages for all species with
                      hourly data.

    Note a symbolic link ('plots/') is set up in html_dir to point to plot_dir

    .. note:: Assumes files are of the format
              gridded_fields/Fieldplot_sites_<model>_
              <short_name>_*<yyyymmddHHMM>.png

    Example using images already in figures directory
    (for O3 at 0Z and 1Z on April 2nd):
    >>> ini_dict = {
    ...    'plot_dir': config.CODE_DIR + "/adaqdocs/figures/adaq_plotting",
    ...    'html_dir': config.CODE_DIR + "/adaqdocs/figures/",
    ...    'short_name_list': ['O3'],
    ...    'models_list': ['aqum_oper','aqum_casestudy'],
    ...    'start_datetime': datetime.datetime(2014,4,2),
    ...    'end_datetime': datetime.datetime(2014,4,2,1)}

    >>> write_aq_plot_sites_website(ini_dict, daqi_only=False)
    ... # doctest: +ELLIPSIS
    Saved website  .../adaqdocs/figures/aq_plot_site_fields_O3.html

    Note this webpage be viewed by opening
    adaqdocs/figures/aq_plot_site_fields_O3.html in a browser.

    Finally check contents of this file:

    >>> webfile = ini_dict['html_dir'] + 'aq_plot_site_fields_O3.html'
    >>> with open(webfile, "r") as fin:
    ...    for line in fin:
    ...        print(line.strip())  # doctest: +ELLIPSIS
    <html>
    <head>
    <title>O3 - site maps</title>
    </head>
    <body>
    <BLANKLINE>
    <BLANKLINE>
    <h1>O3 - site maps</h1>
    <BLANKLINE>
    <table border="0" cellspacing="0" cellpadding="0">
    <tr>
    <td align="center"><b> aqum_oper at sites </b></td>
    <td align="center"><b> aqum_casestudy at sites </b></td>
    <td align="center"><b> Obs at sites </b></td>
    <td align="center"><b> aqum_oper contours </b></td>
    <td align="center"><b> aqum_casestudy contours </b></td>
    </tr>
    <tr>
    <td><a href="plots/gridded_fields/Fieldplot_sites_aqum_oper_O3_\
201404020000.png"><img width=500 src="plots/gridded_fields/\
Fieldplot_sites_aqum_oper_O3_201404020000.png"> </a></td>
    <td><a href="plots/gridded_fields/Fieldplot_sites_aqum_casestudy_O3_\
201404020000.png"><img width=500 src="plots/gridded_fields/\
Fieldplot_sites_aqum_casestudy_O3_201404020000.png"> </a></td>
    <td><a href="plots/gridded_fields/Fieldplot_sites_Obs_O3_\
201404020000.png"><img width=500 src="plots/gridded_fields/\
Fieldplot_sites_Obs_O3_201404020000.png"> </a></td>
    </tr>
    <tr>
    <td>.../Fieldplot_sites_aqum_oper_O3_201404020100.png"> </a></td>
    <td>.../Fieldplot_sites_aqum_casestudy_O3_201404020100.png"> </a></td>
    <td>.../Fieldplot_sites_Obs_O3_201404020100.png"> </a></td>
    </tr>
    </table>
    <BLANKLINE>
    <p><hr><p>
    <font color="#990000"> WARNING: This page is not supported operationally.</font>
    <br><I>Modified:</I> <B><!--#config timefmt="%d %b %Y"--><!--#echo var="LAST_MODIFIED"--> </B>
    <p><hr><p>
    <address>
    </address>
    </body>
    </html>

    """

    #Set up symbolic link to plot directory
    input_dir = ini_dict['plot_dir']
    html_dir = ini_dict.get('html_dir', './')
    if not os.path.exists(html_dir):
        os.makedirs(html_dir)
    sym_dir = 'plots'

    use_models = bool(ini_dict.get('models_list', None))

    #Loop over short_names (one website per short_name)
    for short_name in ini_dict['short_name_list']:
        if daqi_only and ('DAQI' not in short_name):
            #Only produce webpages for DAQI species
            continue

        #Start page
        web = Website(html_dir)
        web.add_symlink(input_dir, sym_dir)
        web.add_header(short_name + ' - site maps')
        web.add_title(short_name + ' - site maps')

        #Main body - table of images
        #Column headers
        headers = []
        if use_models:
            headers = [' ' + label + ' at sites '
                       for label in ini_dict['models_list']]
        headers += [' Obs at sites ']
        if use_models:
            headers += [' ' + label + ' contours '
                        for label in ini_dict['models_list']]

        #Columns: model at sites, obs at sites, model contours
        col_names = []
        if use_models:
            col_names += ['Fieldplot_sites_' + label + '_' + short_name
                          for label in ini_dict['models_list']]
        col_names += ['Fieldplot_sites_Obs_' + short_name]
        if use_models:
            col_names += ['Fieldplot_' + label + '_' + short_name
                          for label in ini_dict['models_list']]

        #Rows: loop over dates (for DAQI this is at midday only)
        row_names = []
        dt = ini_dict['start_datetime']

        if 'DAQI' in short_name:
            dt = dt.replace(hour=12)

        while dt <= ini_dict['end_datetime']:
            row_names.append('*' + dt.strftime("%Y%m%d%H%M.png"))
            if 'DAQI' in short_name:
                dt += datetime.timedelta(days=1)
            else:
                dt += datetime.timedelta(hours=1)

        #Don't include DAQI plots if not requested
        #eg don't allow 'O3_DAQI' if short_name='O3'
        if 'DAQI' not in short_name:
            ignore_string = 'DAQI'
        else:
            ignore_string = None

        #Populate table with filenames of images.
        web.add_images_table(sym_dir+'/gridded_fields', row_names, col_names, width=500,
                             wildcards=False, reverse_row_cols=True,
                             headers=headers, ignore_string=ignore_string)

        #Finish page
        web.add_footer(operational_warning=True)
        web.write_file('aq_plot_site_fields_'
                       + short_name +'.html')


if __name__ == '__main__':

    import doctest
    doctest.testmod()
