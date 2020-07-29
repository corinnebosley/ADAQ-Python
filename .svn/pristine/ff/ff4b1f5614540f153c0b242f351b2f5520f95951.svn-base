"""
Class for reading NAME horizontally gridded files and NAME time series files
into an ADAQData class.
"""
from __future__ import print_function
from __future__ import absolute_import

import datetime
import os
import glob
import re
import warnings

from six.moves.builtins import str

import numpy as np
import cf_units
import iris

import adaq_data
import cube_chemistry
import shell_commands

def add_euconc_to_short_names(short_name_dict):
    """
    Add Eulerian Concentration and Concentration (=Lagrangian+Eulerian)
    to short name dictionary.

    :param short_name_dict: Dictionary whose keys are long names as determined
                             by iris on reading NAME input files and whose
                             value is a a dictionary of an appropriate
                             short-name, units and standard-name to match to.
    :returns: short_name_dict, as input but with extra long-name keys
              of <species>_EULERIAN_CONCENTRATION and <species>_CONCENTRATION
              with the same value dictionary as <species>_AIR_CONCENTRATION
              from the input dictionary.

    >>> short_name_dict = {
    ... 'O3_AIR_CONCENTRATION' : {
    ... 'short_name': 'O3', 'units': cf_units.Unit('ug/m3'),
    ... 'standard_name': 'mass_concentration_of_ozone_in_air'}}
    >>> short_name_dict = add_euconc_to_short_names(short_name_dict)
    >>> short_name_dict == {
    ... 'O3_AIR_CONCENTRATION' : {
    ... 'short_name': 'O3', 'units': cf_units.Unit('ug/m3'),
    ... 'standard_name': 'mass_concentration_of_ozone_in_air'},
    ... 'O3_EULERIAN_CONCENTRATION' : {
    ... 'short_name': 'O3', 'units': cf_units.Unit('ug/m3'),
    ... 'standard_name': 'mass_concentration_of_ozone_in_air'},
    ... 'O3_CONCENTRATION' : {
    ... 'short_name': 'O3', 'units': cf_units.Unit('ug/m3'),
    ... 'standard_name': 'mass_concentration_of_ozone_in_air'}}
    True
    """
    long_names = list(short_name_dict.keys())
    for long_name in long_names:
        sndict = short_name_dict[long_name]
        match = re.search(r'(.+)_AIR_CONCENTRATION', long_name)
        if match:
            eu_long_name = match.group(1)+'_EULERIAN_CONCENTRATION'
            short_name_dict[eu_long_name] = sndict
            tot_long_name = match.group(1)+'_CONCENTRATION'
            short_name_dict[tot_long_name] = sndict

    return short_name_dict


class NAMEData(adaq_data.ADAQData):

    """
    Subclass of ADAQData, which contains extra functionality specific to
    NAME data

    **Gridded Data Example:**

    This example steps through the process of loading a gridded NAME output
    file.

    >>> import config
    >>> directory = config.SAMPLE_DATADIR+'name/'

    Initialise class:

    >>> name = NAMEData()

    Determine which filenames should be read in according to required
    start and end dates:

    >>> filenames = name.get_filenames(directory,
    ...     file_pattern = 'Fields_grid1',
    ...     start_datetime=datetime.datetime(2011, 3, 13),
    ...     end_datetime=datetime.datetime(2011, 3, 18))
    >>> print(filenames[0]) # doctest: +ELLIPSIS
    /.../Fields_grid1_C1_T2_201103130000.txt

    Now read the data from name files into a gridded_cube_list.
    Note filenames, start_datetime, end_datetime and have been saved
    in the name object so no need to give these again.

    >>> filenames = [directory + 'Fields_grid1*.txt']

    >>> gcl = name.readdata(filenames,
    ... field_attributes = {'Species':'CAESIUM-137'})

    Check data is gridded data not time series data

    >>> name.check_gridded_cube()
    short_name - OK
    label - OK
    T axis - OK - coord= time
    X axis - OK - coord= longitude
    Y axis - OK - coord= latitude

    Examine data:

    >>> print(name.gridded_cube_list)
    0: CAESIUM-137_AIR_CONCENTRATION / (Bq / m^3) \
(time: 6; latitude: 90; longitude: 180)

    Have a look at the gridded cube:

    >>> Cs137 = name.extract(gridded=True, singlecube=True,
    ... short_name='CAESIUM-137_AIR_CONCENTRATION')
    >>> print(Cs137)
    ... # doctest: +ELLIPSIS
    CAESIUM-137_AIR_CONCENTRATION / (Bq / m^3) (time: 6; latitude: 90; longitude: 180)
         Dimension coordinates:
              time                                  x            -              -
              latitude                              -            x              -
              longitude                             -            -              x
         Scalar coordinates:
              height: 1000.0 m
         Attributes:
              End of release: 20/04/2011 00:00 UTC
              Ensemble Av: No ensemble averaging
              Horizontal Av or Int: No horizontal averaging
              Met data: NWP Flow.Global_PT1_flow; NWP Flow.Global_PT2_flow;...
              NAME Version: NAME III (version 6.3)
              Name: Unnamed Field Req 5
              Quantity: Air Concentration
              Release height: Multiple Sources
              Release location: Multiple Sources
              Run duration: 91day 9hr 0min
              Run name: Stohl2degFukushimaTake2
              Run time: 10/02/2014 11:37:15.277 UTC
              Source strength: Multiple Sources
              Sources: All sources
              Species: CAESIUM-137
              Species Category: RADIONUCLIDE
              Start of release: 10/03/2011 12:00 UTC
              Time Av or Int: 1day 0hr 0min average
              Vertical Av or Int: No vertical averaging
              label: NAME
              short_name: CAESIUM-137_AIR_CONCENTRATION
         Cell methods:
              mean: time

    **Time series Data Example:**

    This example steps through the process of loading
    a NAME time series file. The data is automatically read into a
    sites_cube_list.

    >>> import config
    >>> directory = config.SAMPLE_DATADIR+'name/'

    Initialise class:

    >>> name_ts = NAMEData()

    Now read the data from name files into a sites_cube_list.

    >>> filenames = [directory + 'Time_series_grid1.txt']

    >>> scl = name_ts.readdata(filenames)

    Examine data:

    >>> print(name_ts.sites_cube_list)
    0: VOLCANIC_ASH_AIR_CONCENTRATION / (g/m3) (site_id: 2; time: 132)

    >>> print(name_ts.sites_cube_list[0])
    VOLCANIC_ASH_AIR_CONCENTRATION / (g/m3) (site_id: 2; time: 132)
         Dimension coordinates:
              site_id                               x        -
              time                                  -        x
         Auxiliary coordinates:
              latitude                              x        -
              longitude                             x        -
              site_abbrev                           x        -
              site_name                             x        -
         Scalar coordinates:
              site_type: Unknown
              source_latitude: 63.63 degrees
              source_longitude: -19.05 degrees
              z: Boundary layer
         Attributes:
              End of release: 0830UTC 28/03/2022
              Forecast duration: 87648 hours
              Met data: NWP Flow.Global_PT1_flow; NWP Flow.Global_PT2_flow
              NAME Version: NAME III (version 6.0 for VAAC)
              Quantity: Air Concentration
              Release height: 8376.000m asl +/- 6864.000m
              Release location: 19.0500W   63.6300N
              Release rate: 2.5890222E+08g/s
              Run time: 1358UTC 30/04/2012
              Species: VOLCANIC_ASH
              Species Category: VOLCANIC
              Start of release: 0830UTC 28/03/2012
              Title: KATLA_20120401_12Z
              label: NAME
              short_name: VOLCANIC_ASH_AIR_CONCENTRATION

    **Back Run Example:**

    The validity time of the output data in NAME needs to be treated
    differently for back runs. Iris assumes that the data is valid from
    the validity time minus the averaging or integrating period to the
    validity time but for back runs the data is valid from the validity
    time to the validity time plus the averaging or integrating period.
    The code in name_data corrects the validity time in the Iris cube.

    >>> import config
    >>> directory = config.SAMPLE_DATADIR+'name/'

    Initialise class:

    >>> name = NAMEData()

    >>> filenames = [directory + 'Fields_grid2*.txt']

    Read data setting back to True:

    >>> gcl = name.readdata(filenames, back=True,
    ... field_attributes = {'Quantity':'Air Concentration'})

    Check time bounds and time point noting that the point should
    now be equal to the lower bound

    >>> np.set_printoptions(formatter={'float':lambda x: '{:5.2f}'.format(x)})
    >>> print(gcl[0].coord('time').bounds)
    [[379508.00 379514.00]]

    >>> print(gcl[0].coord('time').points)
    [379508.00]
    >>> np.set_printoptions()

    **Ensemble Example:**

    name_data also includes a functionality for creating an ensemble
    of all cubes which share the same species and quantity and averaging
    or integrating periods. The ensemble is specified by providing a
    dictionary of the directories containing each member together with
    a ensemble member number to use in the cube's realization coordinate.
    In an Iris cube all dimension coordinates must be numeric. This example
    steps through the loading of an ensemble of NAME output.

    Initialise class:

    >>> name = NAMEData()
    >>> directory = config.SAMPLE_DATADIR+'name_ensemble/member[012]/'
    >>> filenames = [directory + 'Fields_grid4*.txt']

    Read data using an ensemble dictionary:

    >>> ensemble_dict={config.SAMPLE_DATADIR+'name_ensemble/member0/': ['member0', 0],
    ... config.SAMPLE_DATADIR+'name_ensemble/member1/': ['member1', 1],
    ... config.SAMPLE_DATADIR+'name_ensemble/member2/': ['member2', 2]}

    >>> gcl = name.readdata(filenames, ensemble_dict=ensemble_dict,
    ... field_attributes = {'Field Name':'Unnamed Field Req 8'})

    >>> print(gcl[0].coord('realization'))
    DimCoord(array([0, 1, 2]), \
standard_name='realization', units=Unit('1'))

    >>> print(gcl[0].coord('ensemble_member').points)
    ['member0' 'member1' 'member2']

    .. note::
        Requires iris 1.7+

    """

    # Dictionary for translating the long names of NAME cubes being loaded to
    # generic (model independent) short names. Phenomena present in this
    # table can be requested by specifying their short or long names in the
    # short_name_list argument passed to 'readdata'. See '__callback' function
    # for use of units and standard names.
    # Entries are not required for null translations.
    _std_name = 'mass_concentration_of_{}_in_air'
    _SHORT_NAME_DICT = {
        'O3_AIR_CONCENTRATION' : {
            'short_name': 'O3',
            'units': cf_units.Unit('ug/m3'),
            'standard_name': _std_name.format('ozone')},
        'NO2_AIR_CONCENTRATION' : {
            'short_name': 'NO2',
            'units': cf_units.Unit('ug/m3'),
            'standard_name': _std_name.format('nitrogen_dioxide')},
        'SULPHUR-DIOXIDE_AIR_CONCENTRATION' : {
            'short_name': 'SO2',
            'units': cf_units.Unit('ug/m3'),
            'standard_name': _std_name.format('sulfur_dioxide')},
        'SULPHATE_AIR_CONCENTRATION' : {
            'short_name': 'SO4',
            'units': cf_units.Unit('ug/m3'),
            'standard_name': _std_name.format(
                'sulfate_dry_aerosol_particles')},
        'AMMONIA_AIR_CONCENTRATION' : {
            'short_name': 'NH3',
            'units': cf_units.Unit('ug/m3'),
            'standard_name': _std_name.format('ammonia')},
        'NH42SO4_AIR_CONCENTRATION' : {
            'short_name': 'NH42SO4',
            'units': cf_units.Unit('ug/m3'),
            'standard_name': _std_name.format(
                'ammonium_sulfate_dry_aerosol_particles')},
        'NH4NO3_AIR_CONCENTRATION' : {
            'short_name': 'NH4NO3',
            'units': cf_units.Unit('ug/m3'),
            'standard_name': _std_name.format(
                'ammonium_nitrate_dry_aerosol_particles')},
        'PMC_AIR_CONCENTRATION' : {
            'short_name': 'PMC_primary',
            'units': cf_units.Unit('ug/m3'),
            'standard_name': _std_name.format(
                'primary_pm_coarse_dry_aerosol_particles')},
        'PM25_AIR_CONCENTRATION' : {
            'short_name': 'PM2p5_primary',
            'units': cf_units.Unit('ug/m3'),
            'standard_name': _std_name.format(
                'primary_pm2p5_dry_aerosol_particles')},
        'PM10_AIR_CONCENTRATION' : {
            'short_name': 'PM10_primary',
            'units': cf_units.Unit('ug/m3'),
            'standard_name': _std_name.format(
                'primary_pm10_dry_aerosol_particles')},
        'CO_AIR_CONCENTRATION' : {
            'short_name': 'CO',
            'units': cf_units.Unit('ug/m3'),
            'standard_name': _std_name.format('carbon_monoxide')},
        'NO_AIR_CONCENTRATION' : {
            'short_name': 'NO',
            'units': cf_units.Unit('ug/m3'),
            'standard_name': _std_name.format('nitrogen_monoxide')},
        'NAER_AIR_CONCENTRATION' : {
            'short_name': 'NAER',
            'units': cf_units.Unit('ug/m3'),
            'standard_name': _std_name.format(
                'nitrate_coarse_mode_dry_aerosol')},
        'NO3_AIR_CONCENTRATION' : {
            'short_name': 'NO3',
            'units': cf_units.Unit('ug/m3'),
            'standard_name': _std_name.format('nitrate_radical')},
        'HNO3_AIR_CONCENTRATION' : {
            'short_name': 'HNO3',
            'units': cf_units.Unit('ug/m3'),
            'standard_name': _std_name.format('nitric_acid')},
        'OH_AIR_CONCENTRATION' : {
            'short_name': 'OH',
            'units': cf_units.Unit('ug/m3'),
            'standard_name': _std_name.format('hydroxyl_radical')},
        'BD_AIR_CONCENTRATION' : {
            'short_name': 'BD',
            'units': cf_units.Unit('ug/m3'),
            'standard_name': _std_name.format('1,3-butadiene')},
        'C2H4_AIR_CONCENTRATION' : {
            'short_name': 'C2H4',
            'units': cf_units.Unit('ug/m3'),
            'standard_name': _std_name.format('ethene')},
        'C3H6_AIR_CONCENTRATION' : {
            'short_name': 'C3H6',
            'units': cf_units.Unit('ug/m3'),
            'standard_name': _std_name.format('propene')},
        'C4H10_AIR_CONCENTRATION' : {
            'short_name': 'C4H10',
            'units': cf_units.Unit('ug/m3'),
            'standard_name': _std_name.format('butane')},
        'C5H8_AIR_CONCENTRATION' : {
            'short_name': 'C5H8',
            'units': cf_units.Unit('ug/m3'),
            'standard_name': _std_name.format('isoprene')},
        'HCHO_AIR_CONCENTRATION' : {
            'short_name': 'HCHO',
            'units': cf_units.Unit('ug/m3'),
            'standard_name': _std_name.format('formaldehyde')},
        'CH3CHO_AIR_CONCENTRATION' : {
            'short_name': 'CH3CHO',
            'units': cf_units.Unit('ug/m3'),
            'standard_name': _std_name.format('acetaldehyde')},
        'OXYL_AIR_CONCENTRATION' : {
            'short_name': 'OXYL',
            'units': cf_units.Unit('ug/m3'),
            'standard_name': _std_name.format('o-xylene')},
        'TOLUEN_AIR_CONCENTRATION' : {
            'short_name': 'TOLUENE',
            'units': cf_units.Unit('ug/m3'),
            'standard_name': _std_name.format('toluene')},
        'JNO2_AIR_CONCENTRATION' : {
            'short_name': 'jNO2',
            'units': cf_units.Unit('s-1'),
            'standard_name': 'photolysis_rate_of_nitrogen_dioxide'},
        'JO1D_AIR_CONCENTRATION' : {
            'short_name': 'jO1D',
            'units': cf_units.Unit('s-1'),
            'standard_name': 'photolysis_rate_of_ozone_to_1D_oxygen_atom'},
        'OH_CHEMISTRY_FIELD' : {
            'short_name': 'OH',
            'units': cf_units.Unit('ug/m3'),
            'standard_name': _std_name.format('hydroxyl_radical')},
        'OP_CHEMISTRY_FIELD' : {
            'short_name': 'O3P',
            'units': cf_units.Unit('ug/m3'),
            'standard_name': _std_name.format('3P_oxygen_atom')},
        'OD_CHEMISTRY_FIELD' : {
            'short_name': 'O1D',
            'units': cf_units.Unit('ug/m3'),
            'standard_name': _std_name.format('1D_oxygen_atom')},
        # Remaining uncertainty in names below, may require changing in future.
        'TOLP1_CHEMISTRY_FIELD' : {
            'short_name': 'TOLP1',
            'units': cf_units.Unit('ug/m3'),
            'standard_name': _std_name.format('toluene_peroxy_radical')},
        # Remaining uncertainty in names below, may require changing in future.
        'CH2O2C_CHEMISTRY_FIELD' : {
            'short_name': 'CH2O2C',
            'units': cf_units.Unit('ug/m3'),
            'standard_name': _std_name.format('ethylene_peroxy_radical')},
        # Remaining uncertainty in names below, may require changing in future.
        'HO2NO2_CHEMISTRY_FIELD' : {
            'short_name': 'HO2NO2',
            'units': cf_units.Unit('ug/m3'),
            'standard_name': _std_name.format('peroxynitric_acid')},
        # Remaining uncertainty in names below, may require changing in future.
        'BDPEROXY_CHEMISTRY_FIELD' : {
            'short_name': 'BDPEROXY',
            'units': cf_units.Unit('ug/m3'),
            'standard_name': _std_name.format('1,3-butadiene_peroxy_radical')},
        'HO2_CHEMISTRY_FIELD' : {
            'short_name': 'HO2',
            'units': cf_units.Unit('ug/m3'),
            'standard_name': _std_name.format('hydroperoxy_radical')},
        'GRASS_POLLEN_AIR_CONCENTRATION': {
            'short_name': 'grass_pollen',
            'units': cf_units.Unit('g/m3'),
            'standard_name': 'grain_concentration_of_poaceae_pollen_in_air'},
        'BIRCH_POLLEN_AIR_CONCENTRATION': {
            'short_name': 'birch_pollen',
            'units': cf_units.Unit('g/m3'),
            'standard_name': 'grain_concentration_of_betula_pollen_in_air'},
        'OAK_POLLEN_AIR_CONCENTRATION': {
            'short_name': 'oak_pollen',
            'units': cf_units.Unit('g/m3'),
            'standard_name': 'grain_concentration_of_quercus_pollen_in_air'},
        'GRASS_POLLEN_HEAT_SUM': {
            'short_name': 'grass_pollen_heat_sum',
            'units': cf_units.Unit('degree days'),
            'standard_name': 'integral_wrt_time_of_air_temperature_excess'},
        'GRASS_POLLEN_REVISED_SOURCE_STRENGTH': {
            'short_name': 'grass_pollen_emission',
            'units': cf_units.Unit('g'),
            'standard_name': ('tendency_of_atmosphere_mass_content_of_'
                              'poaceae_pollen_due_to_emission')}

        }
    #Also include eulerian and total concentrations in _SHORT_NAME_DICT
    add_euconc_to_short_names(_SHORT_NAME_DICT)

    # Dictionary of derived quantities referenced by short name.
    # Each is the sum of 2 or more NAME species. Species names listed
    # as components are the species short names after any translation
    # from _SHORT_NAME_DICT has been applied.
    DERIVED_SHORT_NAME_DICT = {
        'PM2p5' : {
            'components' : ('PM2p5_primary', 'SO4', 'NH42SO4', 'NH4NO3'),
            'standard_name' : 'mass_concentration_of_pm2p5_dry_aerosol_in_air'
            },
        'PM10' : {
            'components' : ('PM2p5_primary', 'SO4', 'NH42SO4', 'NH4NO3',
                            'PMC_primary', 'NAER'),
            'standard_name' : 'mass_concentration_of_pm10_dry_aerosol_in_air'
            },
        'NOx': {
            'components': ['NO', 'NO2'],
            'standard_name': 'mole_fraction_of_nox_in_air',
            'function': cube_chemistry.calculate_nox
            },
        'Ox': {
            'components': ['O3', 'NO2'],
            'standard_name': 'mole_fraction_of_ox_in_air',
            'function': cube_chemistry.calculate_ox
            }
        }


    def __init__(self, label='NAME'):
        """
        Initiates class as a subset of adaq_data.ADAQData, plus other
        model-specific data
        """

        adaq_data.ADAQData.__init__(self)

        self.filenames = None  # List of NAME filenames
        self.directory = None  # Directory containing NAME files
        self.file_pattern = None  # Pattern for matching files within directory
        self.label = label  # Label
        self.z_levels = None  # Vertical levels
        self.z_leveltype = None  # Vertical coordinate
        self.field_attributes = None  # Other options incl. species, field
        self.start_datetime = None  # Starting datetime
        self.end_datetime = None  # End datetime
        self.back = False  # Back run?
        self.sites_data = None  # sites_data object
        self.short_name_list = None
            # List of short names requested (can include phenomena to be
            # derived from multiple fields)
        self.short_name_list_to_load = None
            # List of short names required to satisfy request
        self.ensemble_dict = None
            # Dictionary matching directories containing ensemble members
            # to lables.


    def readdata(self, filenames=None, z_levels=None, z_leveltype=None,
                 field_attributes=None, start_datetime=None,
                 end_datetime=None, back=False, label='NAME',
                 short_name_list=None, ensemble_dict=None):

        """
        Create an gridded_cube_list from the filename(s) specified.

        :param filenames: list of files to read, can include wildcards etc
        :param z_level: output levels
        :param z_leveltype: units of z_level (e.g. flight_level)
        :param field_attributes: dictionary specifying column headers
                                 (e.g. {'Species','CAESIUM-137'}).
                                 A full description of how to specify
                                 field_attributes can be found
                                 in the documentation.
                                 See :ref:`field_attributes_ref`
        :param start_datetime: datetime object start date
        :param end_datetime: datetime object end date
        :param back: is the output from a back run? - default is False
        :param label: label for the cubes
        :param short_name_list: list of short names to read (can include
            long names of NAME fields, short names associated with these
            and/or associated with phenomena derived from multiple fields)
        :param ensemble_dict: dictionary matching ensemble directories to
                              ensemble labels
        """

        # Change input values to self. values
        if filenames is not None:
            self.filenames = filenames
        if label is not None:
            self.label = label
        if z_levels is not None:
            try:
                z_levels = [float(v) for v in z_levels]
            except:
                z_levels = [str(v) for v in z_levels]
            self.z_levels = z_levels
        if z_leveltype is not None:
            self.z_leveltype = z_leveltype
        if field_attributes is not None:
            self.field_attributes = field_attributes
        if start_datetime is not None:
            self.start_datetime = start_datetime
        if end_datetime is not None:
            self.end_datetime = end_datetime
        if back is not None:
            self.back = back
        if short_name_list is not None:
            self.short_name_list = short_name_list
        if ensemble_dict is not None:
            self.ensemble_dict = ensemble_dict

        if self.filenames is None:
            raise ValueError("ModelCube: no filename(s) specified")

        # Determine attribute constraints
        att_spec = ({} if self.field_attributes is None
                    else self.field_attributes)
        att_constraint = iris.AttributeConstraint(**att_spec)

        constraints = att_constraint

        #Set up short name constraint
        #If short names represent items derived from more than one field
        #in the NAME file, the requested short name is first converted to a
        #list of component short names
        #Names of any derived items requested are saved for later reference
        if self.short_name_list is not None:
            specifications = self.get_derived_sn_specs(self.short_name_list)
            dependencies = set(self.short_name_list)
            for spec in specifications:
                if 'components' in spec:
                    dependencies.update(spec['components'])
            constraints &= iris.AttributeConstraint(
                short_name=lambda c: c in dependencies)
            self.short_name_list_to_load = list(dependencies)

        # Determine z_level constraints
        if self.z_levels is not None:
            # Set z_level constraint based on vertical coordinate
            if 'flight_level' in self.z_leveltype:
                z_constraint = iris.Constraint(flight_level=self.z_levels)
                constraints = constraints & z_constraint

            if 'height' in self.z_leveltype:
                z_constraint = iris.Constraint(height=self.z_levels)
                constraints = constraints & z_constraint

            if 'altitude' in self.z_leveltype:
                z_constraint = iris.Constraint(altitude=self.z_levels)
                constraints = constraints & z_constraint

            if 'z' in self.z_leveltype:
                z_constraint = iris.Constraint(z=self.z_levels)
                constraints = constraints & z_constraint

        # Time Constraints
        if self.start_datetime is not None:
            time_constraint = iris.Constraint(
                time=lambda c: c.point >= self.start_datetime)

            constraints = constraints & time_constraint

        if self.end_datetime is not None:
            time_constraint = iris.Constraint(
                time=lambda c: c.point <= self.end_datetime)

            constraints = constraints & time_constraint

        # Load cubes in
        cube_list = iris.load(self.filenames,
                              constraints,
                              callback=self.__callback)

        if not cube_list:
            raise ValueError('No cubes found')

        # Test if cube is gridded or sites and set appropriately
        try:
            self.check_gridded_cube(gridded_cube_list=cube_list, verbose=False)
            self.gridded_cube_list = cube_list
        except:
            self.sites_cube_list = cube_list

        if self.short_name_list is not None:
            #Create any derived items requested
            self.derive(*specifications)

            #Remove any unwanted cubes used in derivation
            #Note that self.gridded_cube_list or self.sites_cube_list
            # was modified in place, so cube_list is still pointing to it
            sname_constraint = iris.AttributeConstraint(
                short_name=lambda c: c in self.short_name_list)
            cube_list = cube_list.extract(sname_constraint)

            #Check final cube list against requested short name list
            short_names_available = [cube.attributes['short_name']
                                     for cube in cube_list]
            if sorted(short_names_available) != sorted(self.short_name_list):
                warnings.warn(
                    'Missing or duplicate data loaded: '
                    'available cube short names are {}, '
                    'expected {}'.format(
                        sorted(short_names_available),
                        sorted(self.short_name_list)))

        return cube_list

    def __callback(self, cube, field, filename):

        """
        Private method to
            Add a 'short_name' for easier identification
            Remove unnecessary attributes
            Move the source location to coordinates
        """

        # Phenomena present in the 'long_name' to 'short_name' translation
        # table can be requested by their short or long names.
        # If the phenomenon is requested by short name, the cube is renamed
        # (from long name) to its standard name, given as an additional
        # attribute in the translation table, and the short name is recorded
        # as a cube attribute.
        # If the phenomenon is requested by long name, the cube is not renamed
        # and it's 'short_name' attribute will match its long name.
        # N.B. the phenomenon represented by a cube is only treated as
        # present in the translation table if its long name is present AND
        # its units are consistent with those expected.
        cube.attributes['short_name'] = cube.long_name
        translation_active = (self.short_name_list_to_load is not None)
        if translation_active and cube.long_name in self._SHORT_NAME_DICT:
            short_name = self._SHORT_NAME_DICT[cube.long_name]['short_name']
            if short_name in self.short_name_list_to_load:
                # Only change name, and short_name if units are consistent with
                # those expected for this particular short name
                # (Must be convertible, e.g. from g/m3 to ug/m3,
                # rather than in ppb or some other unit)
                if short_name in ['jNO2', 'jO1D']:
                    cube.units = self._SHORT_NAME_DICT[cube.long_name]['units']

                units_are_consistent = cube.units.is_convertible(cf_units.Unit(
                    self._SHORT_NAME_DICT[cube.long_name]['units']))
                if units_are_consistent:
                    cube.attributes['short_name'] = short_name
                    cube.rename(self._SHORT_NAME_DICT[cube.long_name]
                                ['standard_name'])

        # Remove unwanted attributes
        unwanted_keys = ['Number of field cols', 'Number of preliminary cols']
        for key in unwanted_keys:
            if key in cube.attributes:
                del cube.attributes[key]

        cube.attributes['label'] = self.label

        # Check that longitudes are in range -180:360
        if cube.coords('longitude'):
            lon_points = cube.coord('longitude').points
            lon_bounds = cube.coord('longitude').bounds
            minlon = min(lon_points)
            maxlon = max(lon_points)
            if minlon < -180:
                cube.coord('longitude').points = lon_points + 360
                cube.coord('longitude').bounds = lon_bounds + 360
            elif maxlon > 360:
                cube.coord('longitude').points = lon_points - 360
                cube.coord('longitude').bounds = lon_bounds - 360
            # Recompute maximum and minimum longitude
            minlon = min(cube.coord('longitude').points)
            maxlon = max(cube.coord('longitude').points)
            if minlon < -180 or maxlon > 360:
                print("WARNING: Longitudes extend beyond the range -180:360")
                print("         Plots may not appear as expected")

        # Convert source location from attribute to coordinate
        if 'Release location' in cube.attributes:
            rel_loc = cube.attributes['Release location']
            if 'multiple' not in rel_loc.lower():
                sourcelon, sourcelat = get_source_location(rel_loc)

                new_coord = iris.coords.AuxCoord(
                    np.array(sourcelat,
                             dtype=np.dtype('float64')),
                    long_name='source_latitude',
                    units=cf_units.Unit('degrees'))
                cube.add_aux_coord(new_coord)

                new_coord = iris.coords.AuxCoord(
                    np.array(sourcelon,
                             dtype=np.dtype('float64')),
                    long_name='source_longitude',
                    units=cf_units.Unit('degrees'))
                cube.add_aux_coord(new_coord)

        if 'X-Y location' in cube.attributes:
            cube.attributes['Location'] = cube.attributes['X-Y location']
            del cube.attributes['X-Y location']

        if 'Location' in cube.attributes:
            # Add site data to cube in ADAQ sites data format
            # NB. Assumes that always has longitude and latitude coordinates

            # Generate and then add site_id to cube as dim_coord
            site_id = adaq_data.generate_siteids(
                [cube.coord('longitude').points],
                [cube.coord('latitude').points])[0]
            cube.add_aux_coord(iris.coords.DimCoord(
                site_id, long_name='site_id'))

            # Note - assumes maximum string length for location
            # is 30 characters
            new_coord = iris.coords.AuxCoord(
                np.array(cube.attributes['Location'],
                         dtype=np.dtype('|U30')),
                long_name='site_name',
                units='no_unit')
            cube.add_aux_coord(new_coord)
            new_coord = iris.coords.AuxCoord(
                np.array(cube.attributes['Location'][0:5],
                         dtype=np.dtype('|U5')),
                long_name='site_abbrev',
                units='no_unit')
            cube.add_aux_coord(new_coord)
            new_coord = iris.coords.AuxCoord(
                np.array('Unknown', dtype=np.dtype('|U7')),
                long_name='site_type',
                units='no_unit')
            cube.add_aux_coord(new_coord)

            del cube.attributes['Location']

        if self.back:
            # Correct time bounds for back run
            time_bounds = cube.coord('time').bounds
            time_diff = time_bounds[:, 1] - time_bounds[:, 0]
            new_bounds = np.ones(time_bounds.shape)
            new_bounds[:, 0] = time_bounds[:, 1]
            new_bounds[:, 1] = time_bounds[:, 1] + time_diff
            cube.coord('time').bounds = new_bounds

        if self.ensemble_dict is not None and not cube.coords('realization'):
            # Now ensemble members are stored in separate directories several
            # attributes need to be removed from the cube or the cubes won't merge
            unwanted_keys = ['Number of field cols',
                             'Number of preliminary cols',
                             'Run time',
                             'Met data',
                             'Run name',
                             'Title']
            for key in unwanted_keys:
                if key in cube.attributes:
                    del cube.attributes[key]

            # Add an ensemble coordinate - note that the ensemble number is in
            # different locations in different filenames
            # Currently assuming the ensemble name is the same as the folder name
            dir_path = os.path.dirname(filename)

            if dir_path in self.ensemble_dict:
                ensemble_member = self.ensemble_dict[dir_path][0]
                realization = self.ensemble_dict[dir_path][1]
            elif dir_path + '/' in self.ensemble_dict:
                ensemble_member = self.ensemble_dict[dir_path + '/'][0]
                realization = self.ensemble_dict[dir_path + '/'][1]
            else:
                err_m = 'No label has been provided for files in {} {}'\
		         .format(dir_path, self.ensemble_dict.keys())
                raise ValueError(err_m)

            ensemble_coord = iris.coords.AuxCoord(ensemble_member,
                                                  long_name='ensemble_member')
            cube.add_aux_coord(ensemble_coord)
            realization_coord = iris.coords.AuxCoord(realization,
                                                     standard_name='realization')
            cube.add_aux_coord(realization_coord)


    def extract_sites(self, sites_data=None):
        """
        Extract site-specific data from gridded_cube_list into
        sites_cube_list, given

        :param sites_data: site information from :class:`sites_info.SitesInfo`
        """

        if self.sites_data is None:
            self.sites_data = sites_data
        if sites_data is None:
            if self.sites_data is None:
                raise ValueError("No sites requested")
            else:
                sites_data = self.sites_data

        #Assume all cubes have the same coordinate system
        if not self.gridded_cube_list:
            raise ValueError("gridded_cube_list has no cubes")
        gridded_cube = self.gridded_cube_list[0]

        xcoord_name = 'longitude' #Default value
        ycoord_name = 'latitude'

        for axis in ['X', 'Y']:
            for coord in gridded_cube.coords():
                if iris.util.guess_coord_axis(coord) == axis:
                    if coord in gridded_cube.dim_coords:
                        if axis == 'X':
                            xcoord_name = coord.name()
                        if axis == 'Y':
                            ycoord_name = coord.name()

        self.sites_cube_list = self.extract_scl_from_gridded(sites_data,
                                                             xcoord_name,
                                                             ycoord_name)

        return self.sites_cube_list


    def get_fileinfo(self, filename):
        """
        Get validity time information from filename,
        returns information in a dictionary.
        Works on files of format ..._yyyymmddhhmm.txt
        Assumes each file only contains one time.

        >>> import config
        >>> name = NAMEData()
        >>> filename = config.SAMPLE_DATADIR + \
        'name/Fields_grid1_C1_T4_201103150000.txt'
        >>> fileinfo = name.get_fileinfo(filename)
        >>> fileinfo == {'data_valid': datetime.datetime(2011, 3, 15, 0, 0)}
        True
        """

        fileinfo = {}
        # Remove path if attached
        filename = os.path.basename(filename)
        split_file = filename.split('_')
        dt_string = split_file[-1][:-4]
        if dt_string.isdigit():
            valid_datetime = datetime.datetime.strptime(
                dt_string, "%Y%m%d%H%M")
            fileinfo['data_valid'] = valid_datetime
        else:
            fileinfo['data_valid'] = None

        return fileinfo


    def get_filenames(self, directory, file_pattern=None,
                      start_datetime=None, end_datetime=None):

        """
        Get file names from a directory, possibly limited by validity
        times/dates (in datetime format)

        :param directory: string, directory containing files
        :param file_pattern: string, to limit filenames (e.g. Fields_grid1)
        :param start_datetime: datetime object start date
        :param end_datetime: datetime object end date

        >>> import config
        >>> name = NAMEData()
        >>> filenames = name.get_filenames(config.SAMPLE_DATADIR+'name/',
        ...     file_pattern = 'Fields_grid1',
        ...     start_datetime=datetime.datetime(2011,3,14,0,0),
        ...     end_datetime=datetime.datetime(2011,3,15,0,0))
        >>> print(filenames[0]) # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
        /.../python_sample_data/name/Fields_grid1_C1_T3_201103140000.txt

        >>> name = NAMEData()
        >>> filenames = name.get_filenames(config.SAMPLE_DATADIR+'name/Fields_grid1*',
        ...     start_datetime=datetime.datetime(2011,3,14,0,0),
        ...     end_datetime=datetime.datetime(2011,3,15,0,0))
        >>> print(filenames[0]) # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
        /.../python_sample_data/name/Fields_grid1_C1_T3_201103140000.txt

        >>> name = NAMEData()
        >>> filenames = name.get_filenames(config.SAMPLE_DATADIR+'name/',
        ...     file_pattern = 'Fields_grid1')
        >>> print(filenames[0]) # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
        /.../python_sample_data/name/Fields_grid1_C1_T1_201103120000.txt

        >>> name = NAMEData()
        >>> filenames = name.get_filenames(config.SAMPLE_DATADIR+'name')
        >>> print(filenames[0]) # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
        /.../python_sample_data/name/Fields_grid1_C1_T1_201103120000.txt

        >>> name = NAMEData()
        >>> filenames = name.get_filenames(config.SAMPLE_DATADIR+'name_pol*')
        >>> print(filenames[0]) # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
        /.../python_sample_data/name_pollen/20150405T0000Z_Fields_grid1_201504050000.txt

        >>> name = NAMEData()
	>>> directory = config.SAMPLE_DATADIR+'name_vaac_ensemble/member*/Fields_grid99*'
        >>> filenames = name.get_filenames(directory)
        >>> print(filenames[0]) # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
        /.../python_sample_data/name_vaac_ensemble/member0/Fields_grid99_201905141600.txt

        >>> name = NAMEData()
	>>> directory = config.SAMPLE_DATADIR+'name_vaac_ensemble/member*/Fields_grid99*'
        >>> filenames = name.get_filenames(directory,
	...       start_datetime=datetime.datetime(2019,5,14,0,0),
        ...       end_datetime=datetime.datetime(2019,5,15,0,0))
        >>> print(filenames[0]) # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
        /.../python_sample_data/name_vaac_ensemble/member0/Fields_grid99_201905141600.txt

        >>> name = NAMEData()
        >>> filenames = name.get_filenames(config.SAMPLE_DATADIR+'name_vaac_ensemble/member*',
        ...                                file_pattern = 'Fields_grid99')
        >>> print(filenames[0]) # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
	/.../python_sample_data/name_vaac_ensemble/member0/Fields_grid99_201905141600.txt

        """

        # Set values into class
        if directory is not None:
            self.directory = directory
        if start_datetime is not None:
            self.start_datetime = start_datetime
        if end_datetime is not None:
            self.end_datetime = end_datetime

        if file_pattern is None:
            if not glob.glob(directory):
                raise ValueError("directory doesn't exist: "+ self.directory)
            else:
                listdir = glob.glob(directory)
                if os.path.isdir(listdir[0]):
                    self.file_pattern = '*.txt'
                else: # variable "directory" also includes a file pattern
                    self.directory = os.path.dirname(directory)
                    self.file_pattern = os.path.basename(directory)
        else:
            if file_pattern.endswith('.txt'):
                self.file_pattern = '*' + file_pattern
            else:
                self.file_pattern = '*' + file_pattern + '*.txt'

        if not glob.glob(self.directory):
            raise ValueError("self.directory doesn't exist: "+ self.directory)

        if not self.directory.endswith('/'):
            self.directory += '/'

        listfilename = glob.glob(self.directory + self.file_pattern)
        listfilename.sort()

        self.filenames = None

        # Check whether or not we need to select by date
        if self.start_datetime is None and self.end_datetime is None:
            if listfilename:
                self.filenames = listfilename
        else:
            for filename in listfilename:

                fileinfo = self.get_fileinfo(filename)

                # Discard file if all data is before required start time
                if self.start_datetime is not None:
                    if fileinfo['data_valid'] is not None:
                        if fileinfo['data_valid'] < self.start_datetime:
                            continue

                # Discard file if all data is after required end time
                if self.end_datetime is not None:
                    if fileinfo['data_valid'] is not None:
                        if fileinfo['data_valid'] > self.end_datetime:
                            continue

                if self.filenames is None:
                    self.filenames = []
                self.filenames.append(filename)

        if self.filenames is None:
            raise UserWarning("No files found. "+self.directory+self.file_pattern)

        return self.filenames


def get_source_location(rel_loc):
    '''
    Extract source longitude from cube attribute text

    >>> rel_loc = '3.9497W   52.9252N'
    >>> sourcelon, sourcelat = get_source_location(rel_loc)
    >>> print(sourcelon)
    -3.9497
    '''

    slontext, slattext = rel_loc.split()

    northsouth = slattext[-1]
    eastwest = slontext[-1]

    ns = 1
    ew = 1
    if northsouth == 'S':
        ns = -1
    if eastwest == 'W':
        ew = -1

    sourcelat = ns * float(slattext[:-1])
    sourcelon = ew * float(slontext[:-1])

    return sourcelon, sourcelat


def __mass_retrieve_filenames(start_datetime, end_datetime):

    '''
    Get list of file name patterns to match files for each date in the period
    from start_datetime to end_datetime. It is assumed that the file names are
    of the form <run_time>Z_Fields_grid1_<valid_time>.txt.gz where <run_time>
    is a time in the format yyyymmddT0000Z, <valid_time> is a time of validity
    in the format yyyymmddhhmm and <run_time> is the beginning of the day that
    includes <valid_time>.
    N.B. The forecast period is assumed to be 1 day. This function does
    not currently allow for longer forecast periods.
    '''

    assert start_datetime.time() == datetime.time(0) and \
           end_datetime.time() == datetime.time(0)
    assert start_datetime.date() < end_datetime.date()

    file_patterns = []
    date = start_datetime.date()
    while date < end_datetime.date():
        datestr = date.strftime('%Y%m%d')
        date_next = date + datetime.timedelta(days=1)
        datestr_next = date_next.strftime('%Y%m%d')
        file_patterns.append(datestr + 'T0000Z_Fields_grid1_*' +
                             datestr + '????.txt.gz')
        file_patterns.append(datestr + 'T0000Z_Fields_grid1_*' +
                             datestr_next + '0000.txt.gz')
        date = date_next

    return file_patterns


def name_mass_retrieve(outputdir, start_datetime, end_datetime,
                       massdir=None, mooseid=None, runid=None,
                       massretries=0, massretrydelay=60, retrieve=True):

    '''
    Function to set up and perform mass retrieval for NAME-AQ files.

    :param outputdir: String. Directory to retrieve files into.
                      Will be created if does not already exist.
    :param start_datetime: Datetime format of start date (time 00:00).
    :param end_datetime: Datetime format of end date (time 00:00).
    :param massdir: String. Directory on mass to retrieve from.
    :param mooseid: String. Moose id under which the files are archived in the
                    form firstname.lastname (default is user's id).
    :param runid: String. Rose suite name.
    :param massretries: Number of times to retry mass retrieval.
    :param massretrydelay: Sleep time in seconds between each retry.
    :param retrieve: Logical. If set to True retrieves files from mass.
                     If False, returns the commands to be run manually by user.
    :returns: 2 element list containing the UNIX commands required to perform
              the retrieval and uncompress the retrieved files.

    Example of calling, although here retrieve is set to False.

    >>> commands = name_mass_retrieve(
    ...     '$DATADIR/mass_retrieve',
    ...     start_datetime=datetime.datetime(2014,3,26),
    ...     end_datetime=datetime.datetime(2014,3,28),
    ...     mooseid='john.hemmings', runid='mi-aw675',
    ...     massretries=1, massretrydelay=5,
    ...     retrieve=False) # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    Retrieving for period  2014-03-26 00:00:00 2014-03-28 00:00:00
    Moose Retrieval Command:
    moo get -f \
    moose:/adhoc/users/john.hemmings/mi-aw675/20140326T0000Z_Fields_grid1_*20140326????.txt.gz \
    moose:/adhoc/users/john.hemmings/mi-aw675/20140326T0000Z_Fields_grid1_*201403270000.txt.gz \
    moose:/adhoc/users/john.hemmings/mi-aw675/20140327T0000Z_Fields_grid1_*20140327????.txt.gz \
    moose:/adhoc/users/john.hemmings/mi-aw675/20140327T0000Z_Fields_grid1_*201403280000.txt.gz \
    .../mass_retrieve

    >>> print(commands) # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    ['moo get -f
    moose:/adhoc/users/john.hemmings/mi-aw675/20140326T0000Z_Fields_grid1_*20140326????.txt.gz
    moose:/adhoc/users/john.hemmings/mi-aw675/20140326T0000Z_Fields_grid1_*201403270000.txt.gz
    moose:/adhoc/users/john.hemmings/mi-aw675/20140327T0000Z_Fields_grid1_*20140327????.txt.gz
    moose:/adhoc/users/john.hemmings/mi-aw675/20140327T0000Z_Fields_grid1_*201403280000.txt.gz
    .../mass_retrieve',
    'gunzip -f \
    .../mass_retrieve/20140326T0000Z_Fields_grid1_*20140326????.txt.gz
    .../mass_retrieve/20140326T0000Z_Fields_grid1_*201403270000.txt.gz
    .../mass_retrieve/20140327T0000Z_Fields_grid1_*20140327????.txt.gz
    .../mass_retrieve/20140327T0000Z_Fields_grid1_*201403280000.txt.gz']
    '''

    # A single mass retrieval is used here.
    # multiple-get file-number limit from "moo si -l" is currently 10000
    # If the number of files to be retrieved is greater than this then
    # multiple retrievals will be required.

    commands = [] # Will hold list of shell commands

    # Make output directory if does not already exist
    outputdir = os.path.expandvars(outputdir)
    if not os.path.isdir(outputdir):
        os.makedirs(outputdir)

    # Get info for locating the required data on mass
    if massdir is None:
        if mooseid is None:
            # Default is user's own moose id. Attempt to get this from user's
            # moose file.
            try:
                with open(os.path.expanduser('~/.moosedir/moose')) as inp:
                    for line in inp:
                        if line.startswith('user='):
                            mooseid = line.split('=')[1].strip()
                            break
            except IOError as err_msg:
                warnings.warn(err_msg)
        if mooseid is None:
            raise ValueError('Missing moose id for locating data set on mass')
        if runid is None:
            raise ValueError('Missing run id for locating data set on mass')
        massdir = 'moose:/adhoc/users/' + mooseid + '/' + runid

    file_patterns = __mass_retrieve_filenames(start_datetime, end_datetime)
    selectfile = None

    # Do retrieval
    print("Retrieving for period ", start_datetime, end_datetime)
    cmd = shell_commands.call_mass(selectfile, massdir, outputdir,
                                   file_patterns=file_patterns,
                                   massretries=massretries,
                                   massretrydelay=massretrydelay,
                                   retrieve=retrieve)
    commands.append(cmd)

    # Unzip retrieved files
    file_list = [outputdir + '/' + p for p in file_patterns]
    cmd = 'gunzip -f ' + ' '.join(file_list)
    if retrieve:
        shell_commands.call_shell(cmd)
    commands.append(cmd)

    return commands


if __name__ == "__main__":

    import doctest
    doctest.testmod()
