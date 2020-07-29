"""
Functions which can be used for a mixture of python plots
"""
from __future__ import division
from __future__ import print_function

from six.moves.builtins import str

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

from name_functions import name_default_annote, rsmc_annote
from name_functions import rsmc_traj_annote
import config

#: List of distinctive colours (see
#: `full list <http://rawgit.com/pelson/7248780/raw
#:/8e571ff02a02aeaacc021edfa7d899b5b0118ea8/colors.html/>`_
#: of available names)
COLOURS = ['k', 'r', 'b', 'g', 'orange', 'c', 'm',
           'lawngreen', 'slategrey', 'y', 'limegreen', 'purple',
           'lightgrey', 'indigo', 'darkblue', 'plum',
           'teal', 'violet', 'saddlebrown', 'lightpink']

#: Dictionary of conversion from units to latex-style units for prettier
#: labelling.
LATEX = {' > ': r'$>$',
         ' < ': r'$<$',
         '>=': r'$\geq$',
         '_': r'$\_$',
         '%': r'$\%$',
         'ug/m3': r'$\mu g\ m^{-3}$',
         'mg/m3': r'$mg\ m^{-3}$',
         'g/m3': r'$g\ m^{-3}$',
         'g / m^3': r'$g\ m^{-3}$',
         'g s / m^3': r'$g\ s\ m^{-3}$',
         'g/m2': r'$g\ m^{-2}$',
         'g / m^2': r'$g\ m^{-2}$'}


def annotate_plot(cube=None, annote='NAME_defaults',
                  annote_location='right', run_time=None):
    """
    Add annotation to field plots.

    :param cube: Cube from which to obtain information to be used in
                 the annotation.

    :param annote: Text to add to plot. Can also be "NAME_defaults"
                   in which case information about the release will
                   be added.

    :param annote_location: Use 'right' to place text to right of
                            plot and 'below' to place text below
                            plot.

    :param run_time: Run time to be added to annotation of RSMC plots.
    """
    if annote_location == 'right':
        plt.subplot2grid((1, 3), (0, 2), frameon=False)
    elif annote_location == 'below':
        plt.subplot2grid((3, 1), (2, 0), frameon=False)
    else:
        plt.subplot2grid((1, 3), (0, 2), frameon=False)

    # Remove the tick marks and labels from the axes
    plt.tick_params(axis='both', which='both',
                    bottom=False, top=False, labelbottom=False,
                    left=False, right=False, labelleft=False)

    if annote == 'NAME_defaults':
        if cube is None:
            raise ValueError('NAME_defaults annotation requires a cube')
        name_default_annote(cube, annote_location=annote_location)

    elif annote == 'RSMC':
        if cube is None:
            raise ValueError('RSMC annotation requires a cube')
        rsmc_annote(cube)

    elif annote == 'RSMC_trajectories':
        if cube is None:
            raise ValueError('RSMC annotation requires a cube')
        rsmc_traj_annote(cube, run_time=run_time)

    else:
        plt.annotate(annote, xy=(0.1, 0.9), xycoords='axes fraction',
                     fontsize=12, verticalalignment='top')


def _dummy_axes(rect, label, frame=False):
    'Set up dummy axes'
    figure = plt.gcf()
    axes = figure.add_axes(rect, label=label, navigate=False, frame_on=frame)
    axes.set_xticks([])

    axes.set_yticks([])
    return axes


def insert_logo():
    '''
    Insert a logo into the current figure

    If the logo file is not found, this routine just returns silently, otherwise
    the logo file is read and the graphic inserted into the parent plot in the
    top left corner.
    '''
    # Save curent axes
    fig = plt.gcf()
    original_axes = plt.gca()
    # Layout parameters
    inset = 0.01  # Distance from edge to logo
    logo_size = 0.16  # Size of logo (height)
    # Source of logo file
    logo_path = config.CODE_DIR + '/adaqcode/MO_Master_B.jpg'
    # Plot logo if found
    if os.path.exists(logo_path):
        logo = plt.imread(logo_path)
        ypix, xpix = np.shape(logo)[0:2]
        # Plot the logo within dummy axes
        left = inset
        bottom = 1.0 - logo_size + inset
        height = logo_size
        width = height * ypix / xpix
        rect = (left, bottom, width, height)
        logo_axes = _dummy_axes(rect, 'logo')
        logo_axes.imshow(logo, zorder=20)
        # Reset the current axes
        fig.sca(original_axes)
    else:
        print("Logo not found")


def combine_extents(extent_list_list):
    """
    Compute the maximum extent from a list of extents.
    Provided a list of latitude and longitude extents the code
    will return the maximum extent.

    .. note:: Code is only suitable for use with latitude and
              longitude coordinates

    >>> extent_list_list = [[59.3, 72.3, 54.6, 63.0],
    ...     [50.9, 77.3, 55.3, 72.9],]
    >>> combine_extents(extent_list_list)
    [50.9, 77.3, 54.6, 72.9]

    """

    # Find maximum extent
    max_lon = -180.
    min_lon = 360.
    max_lat = -90.
    min_lat = 90.

    for extent_i in extent_list_list:
        min_lon = min(extent_i[0], min_lon)
        max_lon = max(extent_i[1], max_lon)
        min_lat = min(extent_i[2], min_lat)
        max_lat = max(extent_i[3], max_lat)

    # Ensure that new extent doesn't exceed globe
    max_lat = min(90.0, max_lat)
    min_lat = max(-90.0, min_lat)
    max_lon = min(360, max_lon)
    min_lon = max(-180., min_lon)

    if min_lon < 0 and max_lon > 180:
        # Try converting all extents to -180:180
        new_max_lon = -180.
        new_min_lon = 360
        for extent_i in extent_list_list:
            if extent_i[0] > 180:
                new_min_lon = min(extent_i[0] - 360., new_min_lon)
                new_max_lon = max(extent_i[1] - 360., new_max_lon)
            else:
                new_min_lon = min(extent_i[0], new_min_lon)
                new_max_lon = max(extent_i[1], new_max_lon)
        if new_max_lon < 180:
            extent = [new_min_lon, new_max_lon, min_lat, max_lat]
        else:
            # Use a global extent
            extent = [-180, 180, -90, 90]
    else:
        extent = [min_lon, max_lon, min_lat, max_lat]

    return extent


def plot_cmap(cmap, filename=None):
    """
    Plot a colormap.

    :param cmap: valid colormap name, or colormap instance
    :param filename: file to save image to. If None, shows to screen.

    >>> cmap = 'YlGnBu'
    >>> import config
    >>> filename = config.CODE_DIR + "/adaqdocs/figures/cmap.png"
    >>> plot_cmap(cmap, filename=filename)  # doctest: +ELLIPSIS
    Saved figure .../cmap.png

    .. image:: ../adaqdocs/figures/cmap.png
       :scale: 50%

    """

    if isinstance(cmap, str):
        cmap = plt.get_cmap(cmap)
    if isinstance(cmap, colors.ListedColormap):
        ncol = cmap.N
    elif isinstance(cmap, colors.LinearSegmentedColormap):
        ncol = 256
    else:
        raise IOError('Unknown type for colormap: ' + type(cmap))

    x = np.linspace(0, 1, ncol)
    norm = colors.BoundaryNorm(x, ncol - 1)
    gradient = np.vstack((x, x))

    plt.figure(figsize=(8, 3))
    ax = plt.gca()
    ax.imshow(gradient, aspect='auto', cmap=cmap, norm=norm,
              interpolation='nearest')
    ax.set_title('Colormap: ' + cmap.name)
    if filename is None:
        plt.show()
    else:
        plt.savefig(filename)
        print('Saved figure', filename)
        plt.close()


# --- Generic functions that could be applied to any plot ---
def add_legend_belowaxes(ncol=5, fontsize='medium', **kwargs):
    """
    Add legend below current axes.
    With sensible (ie approx < 15 character) text strings in legend,
    5 columns and 3 rows can be easily used.
    :param ncol: Number of columns in legend
    :param fontsize: Fontsize of legend
    """
    ax = plt.gca()
    box = ax.get_position()
    # Adjust axes box by shrinking height by 15% and moving upwards
    ax.set_position([box.x0, box.y0 + box.height * 0.15,
                     box.width, box.height * 0.85])

    # Add legend.
    # Nb bbox_to_anchor(x0,y0,width,height),
    # relative to axes lower left position
    # ncol - allow a maximum of ncol columns in each legend.
    ax.legend(bbox_to_anchor=(-0.1, -0.35 * box.height,
                              1.2, 0.15 * box.height),
              loc='center', columnspacing=0.5,
              ncol=ncol, fontsize=fontsize,
              **kwargs)


def units_str(units):
    """
    Converts units to nice formatted units string.
    """
    unitsstr = str(units)
    if unitsstr in LATEX:
        unitsstr = LATEX[unitsstr]
    return unitsstr


def resize_fig(fig, figsize):
    """
    Force a figure to resize to new figsize.


    :param fig: matplotlib figure object to resize.

    :param figsize: list of numbers to resize figure to.
    """
    fig.set_figheight(figsize[1])
    fig.set_figwidth(figsize[0])


def marker_spacing(time, marker_int):
    '''
    Determine marker spacing and start point based on time
    points

    Used in trajectory plotting
    '''

    time_interval = time.units.num2date(time.points[1]) \
        - time.units.num2date(time.points[0])

    marker_interval_seconds = float(marker_int) * 3600
    marker_interval = marker_interval_seconds / time_interval.total_seconds()

    for it, timepoint in enumerate(time.points):
        thispoint = time.units.num2date(timepoint)
        if np.mod(thispoint.hour, float(marker_int)) == 0:
            marker_start = it
            break

    return marker_start, marker_interval


if __name__ == '__main__':
    import doctest

    doctest.testmod()
