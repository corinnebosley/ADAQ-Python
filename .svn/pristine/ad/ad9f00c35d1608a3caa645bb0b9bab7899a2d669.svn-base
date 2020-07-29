"""
radiological.py

Includes a module for decay correcting data on a cube

Example of use::

   from cube_tools.radiological import decay_correct
   new_cube = decay_correct( cube, half_life, correct_date )

Where cube is a suitable Iris cube, half_life is the half_life of the radionuclide
in seconds and correct_date is a python datetime object for the date to which the 
cube should be corrected (e.g. datetime.datetime(2011,03,11))

.. note::

   This only works for cubes without a time dimension. I.e. where all the data has the same
   valid time

"""
from __future__ import division

import numpy as np
import datetime
import iris

def decay_correct( cube, half_life, correct_date ):

    lamda = np.log(2.0)/half_life

    cubetime = cube.coord('time')

    validtime = cubetime.units.num2date(cubetime.points[0])

    deltaT = correct_date - validtime

    deltaTseconds = deltaT.seconds + (86400.0)*deltaT.days

    decayFactor = np.exp(-lamda*deltaTseconds)

    decayCube = iris.analysis.maths.multiply(cube,decayFactor)
    decayCube.name = cube.name
    decayCube.attributes = cube.attributes
    decayCube.attributes['Decay Correction'] = 'Decay corrected to %s' % correct_date.strftime('%H%M%Z %d/%m/%Y')
    decayCube.convert_units('kBq / m^2')


    return decayCube


def main( ):

    depositfile = '/project/NAME/apnm/NAMEIIIout/PythonHackathon/DataFiles/D3Local_36Day_C1_T1_201104160000.txt'

    deposit = iris.load_cube(depositfile,'I131P_DEPOSITION')

    correct_date = datetime.datetime(2011,6,14,0,0)

    half_life = 6.95e5

    decayCube = decay_correct( deposit, half_life, correct_date )

if __name__ == '__main__':
    main()
