"""
Basic class and code for line plots
"""
from __future__ import print_function

from six.moves.builtins import str
from six.moves.builtins import object

import os
import warnings

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import iris.plot as iplt

import cube_time
import plotting_functions

if not os.getenv('DISPLAY'):
    #Enable DISPLAY if running under cron
    mpl.use('Agg')


def add_mobranding():
    """ Add met office branding.
    (Added as a function in line_plot so monty does not
    have to be imported elsewhere).
    .. note::
        The import statement is included in the subroutine and surrounded
        by a try/except so that the ADAQ python code can be run on
        computers where Monty is not available.
    """
    try:
        import monty
        monty.brand()
    except:
        warnings.warn('Monty not available, branding will not be added')

#--- Class ---
class LinePlot(object):
    """
    Class for holding information about line plots
    """

    def __init__(self):
        """
        Initialise class.
        """

        self.lines = []
        self.fig = None #Figure object

        #Information for labelling plot
        self.title = None #Figure title
        self.ylabel = None
        self.xlabel = None
        self.units = None
        self.phenomena_name = None #Standard name if possible
        self.phenomena_short_name = None
        self.site_name = None
        self.site_type = None
        self.legend = True #Add legend

        #Other information
        self.plotdir = None #Output directory
        self.figformat = 'png'
        self.mobrand = False
        self.semilogy = False #Log y-axis (?)
        self.xlim = None #Override automatically set x-axis
        self.ylim = None #Override automatically set y-axis
        self.gridlines = True #Add gridlines
        self.wrap = False #Wrap colours?

    def set(self, key, value):
        """
        Use this to add any attribute to the class.
        """

        self.__dict__[key] = value


    def add_line(self, cube, x=None, y=None, label=None, colour=None,
                 linestyle='-', linewidth=1.5, marker=None, wrap=False,
                 **kwargs):
        """
        Add a line to the TimeSeriesPlot, based on a cube.
        The cube must have a 'time' dimension. Ideally it should be
        a sites_cube from the ADAQData class.
        Adds intelligently chosen colours and labels if not set.
        Also adds some defaults to the class if not already set.

        Firstly, get a cube:

        >>> import adaq_data
        >>> import config
        >>> sample_data_path = config.SAMPLE_DATADIR+'sites_cube_list/'
        >>> OD = adaq_data.ADAQData()
        >>> scl = OD.load_ts(sample_data_path+'aurn_1days.nc')
        >>> obs_cube = OD.extract(short_name='O3', abbrev='HAR',
        ... singlecube=True)

        When first initalised, the lines attribute is an empty list:

        >>> LP = LinePlot()
        >>> print(LP.lines)
        []

        Now add a line, using the cube:

        >>> LP.add_line(obs_cube)

        Have a look at the line (first element in TSD.lines) that has
        just been added - many defaults have been set from this:

        >>> line = LP.lines[0]
        >>> print(line['cube']) # doctest: +ELLIPSIS
        mass_concentration_of_ozone_in_air / (ug/m3) (time: 25)
             Dimension coordinates:
                  time                                    x
             Scalar coordinates:
                  abbrev: HAR
                  latitude: 51.57110977 degrees
                  longitude: -1.326666594 degrees
                  site_altitude: 137 m
                  site_id: 35867333.14...
                  site_name: Harwell
                  site_type: RURAL
             Attributes:
                  Conventions: CF-1.5
                  label: Obs
                  short_name: O3
                  source: AURN
             Cell methods:
                  mean: time (1 hour)
        >>> print(line['colour'])
        r
        >>> print(line['linestyle'])
        -
        >>> print(line['label'])
        Obs
        """

        #Set intelligent defaults for input keywords if not given
        iline = len(self.lines)
        if label is None:
            if 'label' in cube.attributes.keys():
                #Try and get it from cube first
                label = cube.attributes['label']
            else:
                #Give default
                label = 'Line '+str(iline+1)
        if colour is None:
            #Find first colour not already used
            usedcolours = [line['colour'] for line in self.lines]
            icolour = 1 #Dont include black as a default line colour
            colour = plotting_functions.COLOURS[icolour]
            while colour in usedcolours:
                icolour += 1
                if icolour < len(plotting_functions.COLOURS):
                    colour = plotting_functions.COLOURS[icolour]
                elif wrap:
                    icolour = np.mod(icolour, len(plotting_functions.COLOURS))
                    colour = plotting_functions.COLOURS[icolour]
                else:
                    raise ValueError('Too many colours used: '\
                                     'COLOURS is out of range')
        if x is None:
            if cube.dim_coords:
                #By default, plot 1st dimension coordinate
                x = cube.coord(cube.dim_coords[0].name())
            else:
                #If cube doesn't have any dimension coordinates,
                # eg a scalar cube
                #Just set x up to match size of cube data.
                x = np.arange(cube.data.size)
        if y is None:
            #By default, plot cube data values from cube
            y = cube

        #Set some other defaults for class from cube information
        if self.units is None:
            self.units = cube.units
        if self.phenomena_name is None:
            self.phenomena_name = cube.name()
        if self.phenomena_short_name is None:
            if 'short_name' in cube.attributes:
                self.phenomena_short_name = cube.attributes['short_name']
        if self.site_name is None:
            if cube.coords('site_name'):
                if len(cube.coord('site_name').points) == 1:
                    #Only makes sense to add a site name
                    # if for a single site
                    self.site_name = cube.coord('site_name').points[0]
            elif 'Location' in cube.attributes.keys():
                self.site_name = cube.attributes['Location']
        if self.site_type is None:
            if cube.coords('site_type'):
                if len(cube.coord('site_type').points) == 1:
                    #Only makes sense to add a site type
                    # if for a single site (? - may want all rural sites?)
                    self.site_type = cube.coord('site_type').points[0]


        #Put information into line dictionary
        line = {'cube': cube,
                'x': x,
                'y': y,
                'label': label,
                'colour': colour, #Matplotlib colour code
                'linestyle': linestyle, #Matplotlib linestyle code
                'linewidth': linewidth, #Matplotlib line thickness (weight)
                'marker': marker, #Matplotlib marker style
               }

        #Add any extra keywords into the add_settings dictionary
        # in the line dictionary
        add_settings = {}
        if kwargs is not None:
            for key, value in kwargs.items():
                add_settings[key] = value
                line['add_settings'] = add_settings

        self.lines.append(line)


    def plot(self):
        """
        Basic plotting function.
        Generally should be overwritten by subclasses

        >>> import adaq_data
        >>> import config
        >>> sample_data_path = config.SAMPLE_DATADIR+'sites_cube_list/'
        >>> OD = adaq_data.ADAQData()
        >>> scl = OD.load_ts(sample_data_path+'aurn_1days.nc')
        >>> obs_cube = OD.extract(short_name='O3', abbrev='HAR',
        ... singlecube=True)
        >>> LP = LinePlot()
        >>> LP.add_line(obs_cube)
        >>> LP.add_line(obs_cube*0.5,label='obs*0.5')
        >>> LP.title = 'Basic line plot'
        >>> fig = LP.plot()
        >>> plotdir = config.CODE_DIR + "/adaqdocs/figures"
        >>> LP.save_fig(plotdir) # doctest: +ELLIPSIS
        Saved figure  .../lineplot.png

        .. image:: ../adaqdocs/figures/lineplot.png
          :scale: 75%

        .. Note:: To display figure to screen, just use plt.show()
        """

        if not self.lines:
            raise ValueError("LinePlot: no lines have been added")

        self.fig = plt.figure()
        ax = plt.gca()

        # Setup log y-axis if requested
        if self.semilogy:
            plt.semilogy()

        for line in self.lines:

            iplt.plot(line['x'], line['y'], color=line['colour'],
                      label=line['label'], linestyle=line['linestyle'],
                      linewidth=line['linewidth']
                     )

        #Apply axis limits
        if self.xlim is not None:
            ax.set_xlim(self.xlim)
        if self.ylim is not None:
            ax.set_ylim(self.ylim)

        #Add legend
        if self.legend:
            plotting_functions.add_legend_belowaxes()

        #Add title and axis labels
        if self.title is not None:
            self.title = ax.set_title(self.title)
        if self.xlabel is not None:
            ax.set_xlabel(self.xlabel)
        if self.ylabel is not None:
            ax.set_ylabel(self.ylabel)

        #Add gridlines
        if self.gridlines:
            plt.grid()

        # Apply branding
        if self.mobrand:
            add_mobranding()

        return self.fig

    def get_startend_str(self):
        """
        Get start and end times as a string.
        Only returns valid data if same dates for all cubes in lines.
        (Otherwise returns None)
        """
        startstr = None
        endstr = None
        for line in self.lines:
            if line['cube'].coords('time'):
                startdt, enddt = cube_time.get_startenddt(
                    self.lines[0]['cube'])
                startstr_tmp = startdt.strftime("%d/%m/%Y %H:%M")
                endstr_tmp = enddt.strftime("%d/%m/%Y %H:%M")
                if startstr is None:
                    startstr = startstr_tmp
                else:
                    if startstr != startstr_tmp:
                        print('Lines have different starttimes')
                        startstr = None
                        break
                if endstr is None:
                    endstr = endstr_tmp
                else:
                    if endstr != endstr_tmp:
                        print('Lines have different endtimes')
                        endstr = None
                        break

        return startstr, endstr


    def gen_filename(self):
        """"
        Generate automatic filename
        Note generally should be overwritten by other subclasses functions
        """
        filename = 'lineplot.' + self.figformat
        return filename

    def save_fig(self, plotdir=None, filename=None, verbose=1):
        """
        Save figure to file.

        :param plotdir: Directory to save plot in. If not set,
                        will use plotdir set previously,
                        or default to current directory.
        :param filename: Filename to save plot as.
                         If filename not given, will generate automatically.
        :param verbose: 0 (no extra print output) or
                        1 (standard level of print output)
                        2 (extra level of print output for debugging)

        """

        if plotdir is None:
            if self.plotdir is None:
                self.plotdir = './'
            plotdir = self.plotdir
        #Check finishes with /
        if plotdir[-1] != '/':
            plotdir += '/'

        if filename is None:
            filename = self.gen_filename()

        if not os.path.isdir(plotdir):
            if verbose > 1:
                print('Creating output directory:', plotdir)
            os.makedirs(plotdir)

        self.fig.savefig(plotdir + filename)
        if verbose >= 1:
            print('Saved figure ', plotdir+filename)
        plt.close() #Clear current figure to avoid opening too many at once

if __name__ == '__main__':

    import doctest
    doctest.testmod()
