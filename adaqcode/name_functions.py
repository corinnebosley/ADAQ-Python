'''
Contains functions aimed specifically at NAME data
'''
import warnings

import iris
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import (OffsetImage, AnnotationBbox)
from PIL import Image

import config

def height_text_fmt(cube, text_type='title'):
    '''
    Determines the vertical coordinate and presents it in a
    suitable text format

    Two options are available for the text format
       * title - which produces text suitable for a title or legend
       * filename - which produces text suitable for use in a
         filename. i.e. no spaces or special characters are used

    >>> cubes = iris.load(config.SAMPLE_DATADIR+
    ... 'name/Fields_grid2_201304172000.txt')
    >>> cube = cubes[0]
    >>> height_text_fmt(cube)
    'From 0.0 to 100.0 m agl'
    >>> height_text_fmt(cube, text_type='filename')
    'From0_0to100_0magl'

    .. note:: Code only accepts cubes where the z-coordinate consists
              of a single point

    '''
    z_coord_names = ['z', 'flight_level', 'height', 'altitude']

    z_coord = None

    # First determine name of z coordinate
    for coord in cube.coords():
        if iris.util.guess_coord_axis(coord) == 'Z':
            z_coord = coord
        elif coord.name() in z_coord_names:
            z_coord = coord

    # Next need to verify that the z_coord exists and contains only 1 point
    if z_coord is None:
        warnings.warn('No vertical coordinate stored in cube')
        height_text = ''
        return height_text

    z_points = len(z_coord.points)
    if z_points > 1:
        raise ValueError('Vertical coordinate contains more than 1 point')

    # Finally need to create a useful text string
    # Note that z and flight_levels have no units in IRIS, although FL needs to
    # appear in text, and z levels and NAME III format files have no bounds

    if z_coord.name() == 'flight_level':
        if z_coord.bounds is None:
            height_text = 'FL{:03.0f}'.format(z_coord.points[0])
        else:
            height_text = 'From FL{:03.0f} to FL{:03.0f}'.format(
                z_coord.bounds[0][0],
                z_coord.bounds[0][1])

    elif z_coord.bounds is None:
        height_text = '{}'.format(z_coord.points[0])
        if z_coord.units != 'no_unit':
            height_text += ' {}'.format(z_coord.units)
    else:
        height_text = 'From {} to {} {}'.format(z_coord.bounds[0][0],
                                                z_coord.bounds[0][1],
                                                z_coord.units)

    # Add reference level (agl or asl) if appropriate
    if z_coord.long_name == 'altitude above sea level':
        height_text += ' asl'
    elif z_coord.long_name == 'height above ground level':
        height_text += ' agl'

    if text_type == 'filename':
        height_text = height_text.replace(' ', '').replace('.', '_')

    return height_text


def name_default_annote(ncl, annote_location='right'):

    """
    Add a set of default annotation to NAME plots describing
    the release
    """

    # Find max value
    max_loc = np.unravel_index(ncl.data.argmax(),
                               ncl.data.shape)
    maxvalue = float(ncl[max_loc].data)

    # Build some suitable text
    name_title_txt = 'Release Data'

    name_txt = 'Location: {}\n'.format(ncl.attributes['Release location'])
    name_txt += 'Start: {}\n'.format(ncl.attributes['Start of release'])
    name_txt += 'End: {}\n'.format(ncl.attributes['End of release'])
    if 'Release rate' in ncl.attributes:
        name_txt += 'Rate: {}\n'.format(ncl.attributes['Release rate'])
    elif 'Source strength' in ncl.attributes:
        name_txt += 'Rate: {}\n'.format(ncl.attributes['Source strength'])
    name_txt += 'Height: {}\n'.format(ncl.attributes['Release height'])
    name_txt += 'Pollutant: {}\n\n'.format(ncl.attributes['Species'])
    name_txt += 'Run time: {}\n\n'.format(ncl.attributes['Run time'])
    name_txt += 'Maximum Value = {:.3e} {}'.format(maxvalue, ncl.units)

    # These axes fractions relate to the annotation axes, not the plot axes.
    # This means that the xycoords for the annotation start at the RHS
    # of the plot axes (or colorbar axes if that is on the RHS of the plot).
    if annote_location == 'right':
        xytitle = [0.2, 0.87]
        xybody = [0.2, 0.81]
    elif annote_location == 'below':
        xytitle = [0.1, 1.0]
        xybody = [0.1, 0.81]
    else:
        xytitle = [0.2, 0.87]
        xybody = [0.2, 0.81]

    plt.annotate(name_title_txt, xy=xytitle,
                 xycoords='axes fraction', fontsize=12,
                 weight='bold', verticalalignment='top')
    plt.annotate(name_txt, xy=xybody,
                 xycoords='axes fraction', fontsize=12,
                 verticalalignment='top')


def rsmc_annote(ncl):

    """
    Add a set of default annotation to NAME plots describing
    the release
    """

    # Find max value
    max_loc = np.unravel_index(ncl.data.argmax(),
                               ncl.data.shape)
    maxvalue = float(ncl[max_loc].data)

    # Build some modelling centre text
    centre_txt = 'Issuing Centre: Met Office\n'
    centre_txt += 'Dispersion Model: NAME'

    # Build some suitable text
    name_title_txt = 'Release Data'

    name_txt = 'Location: {}\n'.format(ncl.attributes['Release location'])
    name_txt += 'Start: {}\n'.format(ncl.attributes['Start of release'])
    name_txt += 'End: {}\n'.format(ncl.attributes['End of release'])
    if 'Release rate' in ncl.attributes:
        name_txt += 'Rate: {}\n'.format(ncl.attributes['Release rate'])
    elif 'Source strength' in ncl.attributes:
        name_txt += 'Rate: {}\n'.format(ncl.attributes['Source strength'])
    name_txt += 'Height: {}\n'.format(ncl.attributes['Release height'])
    name_txt += 'Pollutant: {}\n\n'.format(ncl.attributes['Species'])
    name_txt += 'Run time: {}\n\n'.format(ncl.attributes['Run time'])

    if 'source_note' in ncl.attributes:
        name_txt += '{}\n\n'.format(ncl.attributes['source_note'])

    name_txt += 'Contour values may change from chart to chart\n\n'
    if ncl.attributes['short_name'] != 'TimeOfArrival':
        name_txt += 'Maximum Value = {:.3e} {}\n'.format(maxvalue, ncl.units)
    else:
        threshold_coord = ncl.coord('Threshold')
        threshold = threshold_coord.points[0]
        threshold_units = threshold_coord.units
        name_txt += 'Threshold for arrival time: {} {}'.format(threshold,
                                                               threshold_units)

    xycentre = [0.0, 0.88]
    xytitle = [0.0, 0.75]
    xybody = [0.0, 0.71]

    plt.annotate(centre_txt, xy=(xycentre),
                 xycoords='axes fraction', fontsize=11,
                 verticalalignment='top')

    plt.annotate(name_title_txt, xy=(xytitle),
                 xycoords='axes fraction', fontsize=10,
                 weight='bold', verticalalignment='top')
    plt.annotate(name_txt, xy=(xybody),
                 xycoords='axes fraction', fontsize=10,
                 verticalalignment='top')


    # Create custom artists for legend
    source_loc = plt.Line2D((0, 1), (0, 0), color='k', marker='s', linestyle='')
    max_loc = plt.Line2D((0, 1), (0, 0), color='k', marker='*', linestyle='')

    # Create legend from custom artist/label lists
    if ncl.attributes['short_name'] != 'TimeOfArrival':
        plt.gca().legend([source_loc, max_loc],
                         ['Source Location', 'Location of Maximum Value'],
                         numpoints=1, fontsize=10,
                         bbox_to_anchor=(0., 0.22), loc=2)
    else:
        plt.gca().legend([source_loc], ['Source Location'],
                         numpoints=1, fontsize=10,
                         bbox_to_anchor=(0., 0.22), loc=2)

    # Add a logo
    img = Image.open(config.CODE_DIR + '/adaqcode/MO_MASTER_black_mono.png')
    img.thumbnail((277, 87), Image.ANTIALIAS)

    imagebox = OffsetImage(img, zoom=1)
    imagebox.image.axes = plt.gca()
    ab = AnnotationBbox(imagebox, [0.2, 0.95],
                        xybox=[30, 0], xycoords='axes fraction',
                        boxcoords="offset points", frameon=False)
    plt.gca().add_artist(ab)


def rsmc_traj_annote(ncl, run_time=None):

    """
    Add a set of default annotation to NAME plots describing
    the release
    """

    # Build some modelling centre text
    centre_txt = 'Issuing Centre: Met Office\n'
    centre_txt += 'Dispersion Model: NAME'

    # Build some suitable text
    name_title_txt = 'Release Data'

    rel_lat = ncl[0].coord('latitude').points[0]
    if rel_lat > 0:
        rel_lat = '{}N'.format(rel_lat)
    else:
        rel_lat = '{}S'.format(abs(rel_lat))
    rel_lon = ncl[0].coord('longitude').points[0]
    if rel_lon > 0:
        rel_lon = '{}E'.format(rel_lon)
    else:
        rel_lon = '{}W'.format(abs(rel_lon))
    name_txt = 'Location: {}  {}\n'.format(rel_lon, rel_lat)

    t_coord = ncl[0].coord('time')
    start = t_coord.units.num2date(t_coord.points[0])
    start_time = start.strftime('%H:%M UTC %d/%m/%Y ')
    name_txt += 'Start: {}\n'.format(start_time)
    name_txt += 'Heights: 500, 1500, 3000 m agl\n\n'

    if run_time is not None:
        name_txt += 'Run time: {}\n\n'.format(run_time)

    xycentre = [0.0, 0.88]
    xytitle = [0.0, 0.75]
    xybody = [0.0, 0.71]

    plt.annotate(centre_txt, xy=(xycentre),
                 xycoords='axes fraction', fontsize=11,
                 verticalalignment='top')

    plt.annotate(name_title_txt, xy=(xytitle),
                 xycoords='axes fraction', fontsize=10,
                 weight='bold', verticalalignment='top')
    plt.annotate(name_txt, xy=(xybody),
                 xycoords='axes fraction', fontsize=10,
                 verticalalignment='top')


    # Create custom artists for legend
    source_loc = plt.Line2D((0, 1), (0, 0), color='k', marker='s', linestyle='')

    # Create legend from custom artist/label lists
    plt.gca().legend([source_loc], ['Source Location'],
                     numpoints=1, fontsize=10,
                     bbox_to_anchor=(0., 0.22), loc=2)

    # Add a logo
    img = Image.open(config.CODE_DIR + '/adaqcode/MO_MASTER_black_mono.png')
    img.thumbnail((277, 87), Image.ANTIALIAS)

    imagebox = OffsetImage(img, zoom=1)
    imagebox.image.axes = plt.gca()
    ab = AnnotationBbox(imagebox, [0.2, 0.95],
                        xybox=[30, 0], xycoords='axes fraction',
                        boxcoords="offset points", frameon=False)
    plt.gca().add_artist(ab)


if __name__ == '__main__':

    import doctest
    doctest.testmod()
