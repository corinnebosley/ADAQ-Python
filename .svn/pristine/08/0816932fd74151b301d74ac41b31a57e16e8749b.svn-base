"""
trial_winddirgram.py

Trial code to plot an ensemble of wind directions.
Code plots a series of wind sticks to show the variation
in wind direction in the ensemble.

Susan, 19 April 2017
"""
from __future__ import division
from six.moves.builtins import zip

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
import iris
import iris.coords as icoords

_MET_SHORT_NAME = {
    '_TEMPERATURE_(C)': 'temperature',
    '_RELATIVE_HUMIDITY_(%)': 'relative_humidity',
    '_PRECIPITATION_RATE_(MM/HR)': 'precipitation',
    '_CLOUD_AMOUNT_(OKTAS)': 'cloud',
    '_WIND_DIRECTION_(DEGREES)': 'wind_direction',
    '_WIND_SPEED': 'wind_speed',
    '_BOUNDARY_LAYER_DEPTH': 'boundary_layer_depth',
    '_PASQUILL_STABILITY': 'pasquill_stability',
    '_PRESSURE_(PA)': 'pressure',
    '_POTENTIAL_TEMPERATURE_(K)': 'potential_temperature',
    '_SEA_LEVEL_PRESSURE_(PA)': 'sea_level_pressure'
    }

def my_callback(cube, field, filename):

    """
    Function to:
     * Add an ensemble coordinate
     * Rename cubes
     * Remove unwanted coordinates
    """

    unwanted_keys = ['Number of field cols',
                     'Number of preliminary cols',
                     'Run time',
                     'Met data']
    for key in unwanted_keys:
        if key in cube.attributes:
            del cube.attributes[key]

    if cube.long_name in _MET_SHORT_NAME:
        cube.attributes['short_name'] = _MET_SHORT_NAME[cube.long_name]
        cube.rename(cube.attributes['short_name'])

    if not cube.coords('realization'):
        ensemble_number = filename.strip('.txt').split('_')[-1]
        realization = ensemble_number[1:]
        ensemble_coord = icoords.AuxCoord(realization,
                                          standard_name='realization')
        cube.add_aux_coord(ensemble_coord)


def plot_windgram(cube):
    '''
    Module for plotting a wind metogram from a cube
    Currently assumes that cube is correctly formatted (i.e
    contains only wind direction data with two dimensions,
    realization and time). The code also currently assumes
    that there are 8 time steps.

    The code also assumes that the zeroth ensemble member is
    the control and plots this in black.
    '''

    sitelat = [0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25]
    lsize = 5.75

    fig = plt.figure(figsize=[16, 4])
    ax = plt.gca()
    plt.xlim([-3, 45])
    plt.ylim([0, 0.5])
    plt.yticks([])
    plt.xticks([0, 6, 12, 18, 24, 30, 36, 42])
    plt.xlabel('Hours since source start')

    t_coord = cube.coord('time')
    time = t_coord.units.num2date(t_coord.points)
    t_delta = [(t - time[0]).total_seconds()/3600. for t in time]
    t_hours = t_delta[0:48:6]

    for x, y in zip(t_hours, sitelat):
        # calculate the x and y size of the bounding box depending
        # on the total
        piesize = lsize

        # now set up a BBox for the pie chart
        # using bounds so that centred at x, y with
        # a total width of piesize
        bb_data = Bbox.from_bounds(x-piesize/2.,
                                   y-piesize/2., piesize, piesize)
        # transform these data coordinates to figure cordinates
        disp_coords = ax.transData.transform(bb_data)
        fig_coords = fig.transFigure.inverted().transform(disp_coords)
        # now add axes based on these figure cordinates
        ax1 = fig.add_axes(Bbox(fig_coords), projection='polar')

        # now do the polar line plot on this axes
        icube = cube[:, int(x)]
        for ireal, winddir in enumerate(icube.data):
            windrad = winddir * np.pi/180.0
            if ireal == 0:
                ax1.plot([windrad, windrad], [0.0, 1.0],
                         'k', linewidth=2, zorder=4)
            else:
                ax1.plot([windrad, windrad], [0.0, 1.0],
                         'r', linewidth=2)

        ax1.set_theta_zero_location("N")
        ax1.set_theta_direction(-1)
        ax1.set_rticks([])
        ax1.set_xticklabels([])

    plt.show()

def main():
    '''
    Main part of script - change workdir to point at your own data
    '''

    workdir = '/data/users/apdg/python_sample_data/name_ensemble/'
    cube = iris.load_cube(workdir + 'Met_Data_*.txt',
                          'wind_direction',
                          callback=my_callback)

    plot_windgram(cube)

if __name__ == '__main__':
    main()


