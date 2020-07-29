"""
Class for reading NAME trajectory data files into an ADAQData class.
"""
import iris
import cf_units

import adaq_data

# A dictionary of units for some of the variables that might be read
# in on trajectories. These are not provided in the trajectory file
# (except for cloud cover)
_VARIABLE_UNITS = {'Temperature': 'K',
                   'Pressure': 'Pa',
                   'Potential Temp': 'K',
                   'BL Depth': 'm',
                   'Cloud (oktas)': '1',
                   'Rel Humidity': '%',
                   'Wind Speed': 'm/s',
                   'Wind Direction': 'degrees',
                   'U Turb': 'm/s',
                   'V Turb': 'm/s',
                   'W Turb': 'm/s',
                   'U Ambient': 'm/s',
                   'V Ambient': 'm/s',
                   'W Ambient': 'm/s'}

class TrajData(adaq_data.ADAQData):

    """
    Subclass of ADAQData, which contains extra functionality specific to
    NAME trajectory data

    **Example:**

    >>> import config
    >>> directory = config.SAMPLE_DATADIR+'name_trajectory/'

    Initialise class:

    >>> traj = TrajData()

    >>> filenames = [directory + 'Data_Traj_C1*.txt']

    >>> tcl = traj.readdata(filenames)

    Check data is trajectory data and contains all necessary coordinates

    >>> traj.check_trajectory_cube()
    short_name - OK
    label - OK
    X axis - OK - coord= longitude
    Y axis - OK - coord= latitude
    time  dim_coord - OK
    short_name - OK
    label - OK
    X axis - OK - coord= longitude
    Y axis - OK - coord= latitude
    time  dim_coord - OK

    Examine data:

    >>> print(traj.trajectory_cube_list)
    0: U Turb / (m/s)                      (time: 193)
    1: U Turb / (m/s)                      (time: 193)

    Have a look at a single trajectory:

    >>> traj1 = traj.trajectory_cube_list[0]
    >>> print(traj1) # doctest: +ELLIPSIS
    U Turb / (m/s)                      (time: 193)
         Dimension coordinates:
              time                           x
         Auxiliary coordinates:
              Travel Time                    x
              altitude                       x
              latitude                       x
              longitude                      x
         Scalar coordinates:
              PP Index: 2.0
              Release Time: 2016-12-02 08:00:00
              Source: Source_2
         Attributes:
              Met data: NWP Flow.Global_PT1_flow; NWP Flow.Global_PT2_flow; ...
              NAME Version: NAME III (version 6.5.1)
              Run name: London 02122016 0800
              Run time: 08/12/2016 08:27:38.939 UTC
              Trajectory direction: Forward trajectories
              label: TRAJ
              short_name: U Turb

    By default the 'U Turb' variable is loaded because it is present in
    all trajectory files. However, it is possible
    to select any variable from the trajectory file to load. e.g.:

    >>> tcl = traj.readdata(filenames, short_name='Temperature')
    >>> print(tcl)
    0: Temperature / (K)                   (time: 193)
    1: Temperature / (K)                   (time: 193)

    .. note::
        Note that each cube contains only one trajectory so each element
        of the cubelist tcl contains a different trajectory. Here we
        have loaded two trajectories.

    """

    def __init__(self, label='NAME', short_name='U Turb'):
        """
        Initiates class as a subset of adaq_data.ADAQData, plus other
        model-specific data
        """

        adaq_data.ADAQData.__init__(self)

        self.filenames = None #List of NAME filenames
        self.label = label #Label
        self.short_name = short_name #Field to load

    def readdata(self, filenames=None, label='TRAJ', short_name='U Turb'):
        """
        Create an trajectory_cube_list from the filename(s) specified.

        :param filenames: list of files to read, can include wildcards etc
        :param label: label for the cubes

        """

        if filenames is not None:
            self.filenames = filenames
        if short_name is not None:
            self.short_name = short_name
        if label is not None:
            self.label = label

        if self.filenames is None:
            raise ValueError("ModelCube: no filename(s) specified")

        # Load cubes in
        self.trajectory_cube_list = iris.load(self.filenames,
                                              short_name,
                                              callback=self.__callback)

        if not self.trajectory_cube_list:
            raise ValueError('No cubes found')

        return self.trajectory_cube_list


    def __callback(self, cube, field, filename):

        """
        Private method to
            Add a 'short_name' for easier identification
            Add a 'label' required in all ADAQ data objects
            Remove unnecessary coordinates
            Try to add units if possible
        """

        cube.attributes['short_name'] = cube.long_name
        cube.attributes['label'] = self.label
        cube.remove_coord('Puff?')

        # Add some units
        if cube.long_name in _VARIABLE_UNITS:
            cube.units = cf_units.Unit(_VARIABLE_UNITS[cube.long_name])

        # If z_coord is pressure add units
        for coord in cube.coords():
            if coord.name() == 'Z (Pa)':
                cube.coord('Z (Pa)').units = cf_units.Unit('Pa')

if __name__ == "__main__":

    import doctest
    doctest.testmod()
