"""
Module to hold bespoke colormaps

A number of bespoke colormaps are defined in this module.

They can be used with other ADAQ python code by importing
this module and calling the required cmap. for example:

>>> cmap = rsmcdep_cmap()

Or they can be added to an ini_dict as follows:
>>> ini_dict = {}
>>> ini_dict['cmap'] = rsmcac_cmap()

"""
from __future__ import division

import numpy as np
from matplotlib import colors, cm

def rsmcdep_cmap():
    '''
    Generate colour map based on Toulouse deposition map

    >>> cmap = rsmcdep_cmap()

    To plot and display colour map:

    >>> import config
    >>> import plotting_functions
    >>> filename = config.CODE_DIR + "/adaqdocs/figures/rsmcdep_cmap.png"
    >>> plotting_functions.plot_cmap(cmap, filename=filename)
    ... # doctest: +ELLIPSIS
    Saved figure .../rsmcdep_cmap.png

    .. image:: ../adaqdocs/figures/rsmcdep_cmap.png
       :scale: 50%
    '''
# pylint: disable=bad-whitespace
    colours = [[254, 254, 100],
               [154, 254, 154],
               [102, 153, 254],
               [153,   0, 153]]
# pylint: enable=bad-whitespace
    colours = np.array(colours)/255.
    cmap = colors.ListedColormap(colours, name='RSMCDEP')

    #Register colour map, so can be called
    # eg plt.get_cmap('RSMCDEP')
    cm.register_cmap(cmap=cmap)

    return cmap


def rsmcac_cmap():
    '''
    Generate colour map based on Toulouse air concentration map

    >>> cmap = rsmcac_cmap()

    To plot and display colour map:

    >>> import config
    >>> import plotting_functions
    >>> filename = config.CODE_DIR + "/adaqdocs/figures/rsmcac_cmap.png"
    >>> plotting_functions.plot_cmap(cmap, filename=filename)
    ... # doctest: +ELLIPSIS
    Saved figure .../rsmcac_cmap.png

    .. image:: ../adaqdocs/figures/rsmcac_cmap.png
       :scale: 50%
    '''
# pylint: disable=bad-whitespace
    colours = [[254, 252, 106],
               [250, 205,  87],
               [200, 150,  29],
               [149,  89,   9]]
# pylint: enable=bad-whitespace
    colours = np.array(colours)/255.
    cmap = colors.ListedColormap(colours, name='RSMCAC')

    #Register colour map, so can be called
    # eg plt.get_cmap('RSMCAC')
    cm.register_cmap(cmap=cmap)

    return cmap


def rsmctoa_cmap():
    '''
    Generate colour map for arrival times.
    First color must be grey

    >>> cmap = rsmctoa_cmap()

    To plot and display colour map:

    >>> import config
    >>> import plotting_functions
    >>> filename = config.CODE_DIR + "/adaqdocs/figures/rsmctoa_cmap.png"
    >>> plotting_functions.plot_cmap(cmap, filename=filename)
    ... # doctest: +ELLIPSIS
    Saved figure .../rsmctoa_cmap.png

    .. image:: ../adaqdocs/figures/rsmctoa_cmap.png
       :scale: 50%
    '''
# pylint: disable=bad-whitespace
    colours = [[224, 224, 224],
               [215,  48,  31],
               [252, 141,  89],
               [253, 204, 138],
               [254, 240, 217]]
# pylint: enable=bad-whitespace
    colours = np.array(colours)/255.
    cmap = colors.ListedColormap(colours, name='RSMCTOA')

    #Register colour map, so can be called
    # eg plt.get_cmap('RSMCTOA')
    cm.register_cmap(cmap=cmap)

    return cmap

def rsmcgrey_cmap():
    '''
    Generate colour map for faxing.
    All colours are shades of grey

    >>> cmap = rsmcgrey_cmap()

    To plot and display colour map:

    >>> import config
    >>> import plotting_functions
    >>> filename = config.CODE_DIR + "/adaqdocs/figures/rsmcgrey_cmap.png"
    >>> plotting_functions.plot_cmap(cmap, filename=filename)
    ... # doctest: +ELLIPSIS
    Saved figure .../rsmcgrey_cmap.png

    .. image:: ../adaqdocs/figures/rsmcgrey_cmap.png
       :scale: 50%
    '''
# pylint: disable=bad-whitespace
    colours = [[192, 192, 192],
               [128, 128, 128],
               [ 64,  64,  64],
               [  0,   0,   0]]
# pylint: enable=bad-whitespace
    colours = np.array(colours)/255.
    cmap = colors.ListedColormap(colours, name='RSMCGREY')

    #Register colour map, so can be called
    # eg plt.get_cmap('RSMCGREY')
    cm.register_cmap(cmap=cmap)

    return cmap


if __name__ == '__main__':

    import doctest
    doctest.testmod()
