"""
Contains functions relating to chemistry performed on a cube.
"""
from __future__ import division
from __future__ import print_function

import re
import warnings
import iris
import numpy as np
import cf_units
from scipy.constants import R, Boltzmann, Avogadro

import chemistry
from constants import STD_P, STD_T, MMASS_AIR

def convert_cubelist_units(cl, req_units='ug/m3', at_stp=True,
                           chem_units=None, aerosol_units=None):
    """
    Converts units of a supplied cubelist. Currently only does conversion
    to mass concentrations in ug/m3, or to volume mixing ratios in ppb,
    or to number of molecules concentration in cm-3.
    Desired units can be specified for chemical and aerosol species separately.

    The assumption of standard pressure and temperature (stp) may be made.
    If this assumption is not made (e.g. for conversion of 3D vertical profile
    data), the supplied cubelist must contain pressure and temperature fields
    in addition to the chemical species.

    :param cl: cubelist containing chemical and aerosol species to be converted.
               For non-stp conditions this must also contain p and T cubes.
    :param req_units: units to convert to, defaulting to ug/m3. Overridden by
                      the more specific parameters chem_units and aerosol_units
                      if they are supplied.
    :param chem_units: the units that chemical species are converted to.
                       Default req_units. Can be set to 'ug/m3', 'ppb' or 'cm-3'
    :param aerosol_units: the units that aerosol species are converted to.
                          Default req_units. Currently only 'ug/m3' accepted.
    :param at_stp: Set to True for stp assumption, otherwise set False and
                   ensure cl includes the p and T cubes required to carry out
                   the conversion.
    :returns: cl, the same cubelist that was input, but the chemical
              species have units converted.

    Note that this function is currently called from
    :class:`adaq_functions.unit_conversion` only,
    and so the doctests for that function sufficiently test this one.
    """

    #default to using the same units
    if chem_units is None:
        chem_units = req_units
    assert chem_units in ("ug/m3", "ppb", "cm-3")
    if aerosol_units is None:
        aerosol_units = req_units
    assert aerosol_units == "ug/m3"

    #attempt to extract pressure and temperature data
    pressure = None
    temperature = None
    if not at_stp:
        try:
            pressure = cl.extract('air_pressure', strict=True)
        except iris.exceptions.ConstraintMismatchError:
            pass
        try:
            temperature = cl.extract('air_temperature', strict=True)
        except iris.exceptions.ConstraintMismatchError:
            pass

    #attempt to convert each cube
    for i, cube in enumerate(cl):
        cubename = cube.name()

        #skip these without producing a warning, because we already know their
        #purpose in the list
        if cubename in ("air_pressure", "air_temperature"):
            continue

        #This only makes sense as ppb, due to the mixing of different compounds
        if cubename == 'mass_fraction_of_total_voc_in_air':
            continue

        #ensure we are not expressing a compound as an element in it
        try:
            cube = cube_elemental_to_molecular_mass(cube)
        except ValueError:
            pass

        cubetype, molecule = cube_type(cubename)
        if cubetype is not None:
            #to know which conversion function to use, we need the units
            #according to whether this is chemical or an aerosol
            #most convenient check is to try and find a molecular mass
            if cubetype == "mole_fraction": #ie definitely not an aerosol
                req_units = chem_units
            else:
                try:
                    chemistry.get_molmass(molecule)
                    req_units = chem_units
                except ValueError:
                    req_units = aerosol_units

            #to mass concentration
            if req_units == "ug/m3":
                cube = cube_to_mass_conc_in_air(
                    cube, ounits=req_units,
                    at_stp=at_stp, pressure=pressure, temperature=temperature)

            #to volume mixing ratio as ppb
            elif req_units == "ppb":
                cube = cube_to_vmr_in_air(
                    cube, ounits=req_units,
                    at_stp=at_stp, pressure=pressure, temperature=temperature)

            #to number of molecules concentration
            elif req_units == "cm-3":
                cube = cube_to_molecules_conc_in_air(
                    cube, ounits=req_units,
                    at_stp=at_stp, pressure=pressure, temperature=temperature)

            #any other target units are not currently supported
            else:
                raise ValueError("Cannot convert to " + req_units)

        else:
            warnings.warn('Not converting units for ' + cubename)

        cl[i] = cube

    return cl

def cube_elemental_to_molecular_mass(cube):
    """
    This function takes a cube which holds a mass fraction of
    a compound in air, expressed as the mass of an element in
    that compound and converts it to a cube expressed as
    the molecular mass.

    e.g. convert from
    mass_fraction_of_sulfur_dioxide_expressed_as_sulfur_in_air
    (i.e. expressed as an elemental mass)
    and converts to:
    mass_fraction_of_sulfur_dioxide_in_air
    (i.e. expressed as the molecular mass)

    >>> import adaq_data
    >>> import config
    >>> import datetime
    >>> import pp_data
    >>> import numpy as np
    >>> start_dt = datetime.datetime(2014, 4, 6, 0)
    >>> end_dt = datetime.datetime(2014, 4, 6, 1)
    >>> sample_data_path = config.SAMPLE_DATADIR+'aqum_output/oper/'
    >>> ppfiles = [sample_data_path+'aqum_20140406_T+000_QM18.pp']
    >>> pp = pp_data.PPData ()
    >>> pp.readdata(ppfiles, short_name_list=['SO2','O3'])
    [<iris 'Cube' of mass_fraction_of_sulfur_dioxide_expressed_as_sulfur\
_in_air / (kg/kg) (time: 12; grid_latitude: 182; grid_longitude: 146)>,
    <iris 'Cube' of mass_fraction_of_ozone_in_air \
/ (kg kg-1) (time: 12; grid_latitude: 182; grid_longitude: 146)>]

    >>> so2 = pp.extract(short_name='SO2',gridded=True,singlecube=True)
    >>> print(so2.name())
    mass_fraction_of_sulfur_dioxide_expressed_as_sulfur_in_air

    >>> print('{:.3e} {:.3e}'.format(np.amin(so2.data), np.amax(so2.data)))
    1.922e-14 2.586e-08

    Convert SO2 cube to molcular mass

    >>> so2 = cube_elemental_to_molecular_mass(so2)
    >>> print(so2.name())
    mass_fraction_of_sulfur_dioxide_in_air

    >>> print('{:.3e} {:.3e}'.format(np.amin(so2.data), np.amax(so2.data)))
    3.840e-14 5.166e-08

    If the code tries to convert something which
    is not a mass expressed as elemental mass, return the original cube

    >>> o3 = pp.extract(short_name='O3',gridded=True,singlecube=True)
    >>> o3 = cube_elemental_to_molecular_mass(o3)
    >>> print(o3.name())
    mass_fraction_of_ozone_in_air
    """

    # get the cube name
    name = cube.name()

    # check if the name contains expressed_as
    pattern = r"_expressed_as_"
    match = re.search(pattern, name)

    if match:

        # pass the name to the chemistry
        # function which returns conversion factor
        # and new name
        try:
            convfac, newname = chemistry.elemental_to_molecular(name)
        except:
            raise ValueError('Error in cube_elemental_to_molecular_mass.'
                             + 'Unable to convert: ' + name)

        # finally change the mass of the cube and its name
        cube.data = cube.data * convfac
        cube.rename(newname)

    return cube

def cube_to_mass_conc_in_air(cube, ounits='ug/m3',
                             at_stp=False, pressure=None, temperature=None):
    """
    This function takes a cube as input and, if possible,
    returns a new cube which is transformed to mass concentration
    of the species in air. If not possible, ValueError is raised.
    If the input cube is already in mass concentration, returns a
    copy with the correct units.

    :param cube: the cube to copy and transform to mass concentration. Required.
    :param at_stp: whether to assume standard temperature and pressure when \
    calculating density for conversions. Optional argument. Default is False.
    :param pressure: a cube of pressure, required to calculate density \
    (unless at_stp is True)
    :param temperature: a cube of pressure, required to calculate density \
    (unless at_stp is True)
    :param ounits: output units for cube. Optional. Default is ug/m3 but \
    anything supported by udunits for a mass per unit volume is allowed.

    Examples:

    Load some data:

    >>> import adaq_data
    >>> import config
    >>> import datetime
    >>> import pp_data
    >>> start_dt = datetime.datetime(2014, 8, 31, 0)
    >>> end_dt = datetime.datetime(2014, 9, 4, 0)
    >>> sample_data_path = config.SAMPLE_DATADIR+'aqum_output/oper/'
    >>> ppfiles = [sample_data_path+'prodm_op_aqum_20140901_18.000.pp']
    >>> pp = pp_data.PPData()
    >>> pp.readdata(ppfiles, start_datetime=start_dt,
    ...     end_datetime=end_dt, short_name_list=['p','T','O3'])
    [<iris 'Cube' of air_pressure / (Pa) \
(time: 2; model_level_number: 38; grid_latitude: 182; grid_longitude: 146)>,
    <iris 'Cube' of air_temperature / (K) \
(time: 2; model_level_number: 38; grid_latitude: 182; grid_longitude: 146)>,
    <iris 'Cube' of mass_fraction_of_ozone_in_air / (kg kg-1) \
(time: 2; model_level_number: 38; grid_latitude: 182; grid_longitude: 146)>]

    >>> pressure = pp.gridded_cube_list[0]
    >>> temperature = pp.gridded_cube_list[1]
    >>> o3 = pp.gridded_cube_list[2]

    Default unit for output is ug/m3

    >>> o3_mc = cube_to_mass_conc_in_air(o3, at_stp=True)
    >>> print(o3_mc.name(), o3_mc.units)
    mass_concentration_of_ozone_in_air ug/m3

    If the data is already mass concentration this is not an error
    but units are converted. Any output units with the
    correct dimensions are okay if udunits can convert to them.

    >>> newcube = cube_to_mass_conc_in_air(o3_mc, at_stp=True,
    ... ounits='troy_ounce/ft3')
    >>> print('After: ', newcube.name(), newcube.units)
    After:  mass_concentration_of_ozone_in_air troy_ounce/ft3

    If the output units do not have the correct dimensions this is an error.

    >>> newcube = cube_to_mass_conc_in_air(o3, at_stp=True, ounits='mm')
    Traceback (most recent call last):
      File ".../cube_chemistry.py", line 477, in cube_to_mass_conc_in_air
        raise ValueError('cannot convert to '+ounits)
    ValueError: cannot convert to mm

    If don't want to do conversion at STP then need to
    to pass in p and T as arguments.

    >>> newcube = cube_to_mass_conc_in_air(o3, at_stp=False, pressure=pressure,
    ... temperature=temperature, ounits='ug/m3')
    >>> print('After: ', newcube.name(), newcube.units)
    After:  mass_concentration_of_ozone_in_air ug/m3

    If a cube is not convertible to mass concentration,
    this is an error.

    >>> newcube = cube_to_mass_conc_in_air(temperature)
    Traceback (most recent call last):
      File ".../cube_chemistry.py", line 463, in cube_to_mass_conc_in_air
        + cubename)
    ValueError: Error in cube_to_mass_conc_in_air. Cannot transform: \
air_temperature

    For species with no molecular formula, it may still be possible to
    express the mass concentration in terms of one of the components.
    For example:

    >>> pp = pp_data.PPData()
    >>> pp.readdata(ppfiles, start_datetime=start_dt,
    ...     end_datetime=end_dt, short_name_list=['NOx', 'Ox'])
    [<iris 'Cube' of mole_fraction_of_nox_in_air / (ppb) \
(time: 2; model_level_number: 38; grid_latitude: 182; grid_longitude: 146)>,
    <iris 'Cube' of mole_fraction_of_ox_in_air / (ppb) \
(time: 2; model_level_number: 38; grid_latitude: 182; grid_longitude: 146)>]

    >>> nox = pp.extract(short_name='NOx', gridded=True, singlecube=True)
    >>> print(nox.name(), nox.units)
    mole_fraction_of_nox_in_air ppb
    >>> nox_as_no2 = cube_to_mass_conc_in_air(nox, at_stp=True)
    >>> print(nox_as_no2.name(), nox_as_no2.units)
    mass_concentration_of_nox_expressed_as_nitrogen_dioxide_in_air ug/m3

    >>> ox = pp.extract(short_name='Ox', gridded=True, singlecube=True)
    >>> print(ox.name(), ox.units)
    mole_fraction_of_ox_in_air ppb
    >>> ox_as_o = cube_to_mass_conc_in_air(ox, at_stp=True)
    >>> print(ox_as_o.name(), ox_as_o.units)
    mass_concentration_of_ox_expressed_as_ground_state_atomic_oxygen_in_air ug/m3

    """

    # There are a limited number of phenomena types we are
    # able to transform to mass concentration in air. Test
    # each of these in turn.

    cubename = cube.name()
    cubetype, molecule = cube_type(cubename)

    # first check if it is already mass concentration in air
    if cubetype == 'mass_concentration':
        mc_cube = cube.copy()

    # transform from mass mixing ratio
    elif cubetype == 'mass_fraction':
        mc_cube = cube.copy()

        # convert mass mixing ratio to kg/kg
        mc_cube.convert_units('kg/kg')

        # transform cube data. mass concentration = mmr * density of air

        # density of air = pressure * molar mass of air
        #                    / (universal gas constant * temperature)

        # Two options here - use standard temperature and pressure or
        # use cubes of pressure and temperature

        # Calculate at STP (frequently used for surface obs)
        # note that MMASS_AIR is in g/mol so need to divide by 1000
        if at_stp:
            mc_cube.data = (mc_cube.data * STD_P *
                            (MMASS_AIR/1000.) / (R * STD_T))
        else:
            _ensure_p_t_compatible(pressure, temperature, mc_cube)

            mc_cube.data = (mc_cube.data * pressure.data *
                            (MMASS_AIR/1000.) / (R * temperature.data))

        # fix units - have used MKS so output is in kg/m3
        mc_cube.units = 'kg/m3'

        # work out new new name for cube:
        # "mass_concentration_of_X_in_air"
        new_name = "mass_concentration_of_"+molecule+"_in_air"

        # rename the cube
        mc_cube.rename(new_name)

    elif cubetype == 'mole_fraction':
        mc_cube = cube.copy()

        # convert volume mixing ratio to ppb
        mc_cube.convert_units('ppb')

        #With standard units, we have
        # mc = vmr (pressure * molmass) / (R * temperature)
        #With our preferred units:
        #- 10^-3 for molmass (g/mol -> kg/mol)
        #- 10^-9 for vmr (ppb = nmol/mol -> mol/mol)
        #- 10^9 for mc (kg/m3 -> ug/m3)
        # mc = 10^-3 vmr (pressure * molmass) / (R * temperature)

        try:
            if molecule in chemistry.EXPRESSED_AS:
                molmass = chemistry.get_molmass(chemistry.EXPRESSED_AS[molecule])
                molecule = (molecule
                            + "_expressed_as_"
                            + chemistry.EXPRESSED_AS[molecule].replace(" ", "_"))
            else:
                molmass = chemistry.get_molmass(molecule)
        except ValueError:
            raise ValueError("cannot convert non-chemical species "
                             + molecule + " to mass concentration")

        if at_stp:
            mc_cube.data = mc_cube.data * (
                1e-3 * (STD_P * molmass) / (R * STD_T))
        else:
            _ensure_p_t_compatible(pressure, temperature, mc_cube)

            mc_cube.data = mc_cube.data * (
                1e-3 * (pressure.data * molmass) / (R * temperature.data))

        # set the new units
        mc_cube.units = 'ug/m3'

        # rename the cube
        mc_cube.rename("mass_concentration_of_"+molecule+"_in_air")

    # if we get here then we can't convert to mass conc
    else:
        raise ValueError(
            'Error in cube_to_mass_conc_in_air. Cannot transform: '
            + cubename)

    # set units to those requested by ounits - default is ug/m3
    try:
        mc_cube.convert_units(ounits)
    except:
        raise ValueError('cannot convert to '+ounits)

    # now return new cube
    return mc_cube

def cube_to_vmr_in_air(cube, ounits='ppb',
                       at_stp=False, pressure=None, temperature=None):
    '''
    Convert a cube to volume mixing ratio. By default calculates ppb,
    but can also accept any mole fraction, e.g. umole/mole (= ppm)

    If the input cube is in mass mixing ratio, this is easy to convert
    and we don't need pressure or temperature. For concentrations, however,
    we do. Standard temperature and pressure can be used for convenience.

    :param cube: a cube containing the data to be converted.
    :param ounits: output units, defaulting to ppb.
    :param at_stp: whether or not to use standard temperature and pressure.
                   Only applies to conversions from mass concentration.
                   Default False.
    :param pressure: cube containing pressure data, needed for conversions from
                     mass concentration when at_stp is set to False (otherwise
                     ignored).
    :param temperature: cube containing temperature data, needed for conversions
                        from mass concentration when at_stp is set to False
                        (otherwise ignored).
    :returns: vmr_cube, containing the converted data.

    Examples:

    Load some data:

    >>> import config
    >>> import pp_data
    >>> sample_data_path = config.SAMPLE_DATADIR+'aqum_output/oper/'
    >>> ppfiles = [sample_data_path+'aqum_20140325_T+000_QM18.pp']
    >>> pp = pp_data.PPData()
    >>> gcl = pp.readdata(ppfiles, short_name_list=['O3'])
    >>> o3_mmr = gcl[0]
    >>> print("{} {} {:.3e}".format(o3_mmr.name(), o3_mmr.units,
    ...                             o3_mmr.data.max()))
    mass_fraction_of_ozone_in_air kg kg-1 9.434e-08

    Default unit for output is ppb.

    >>> o3_ppb = cube_to_vmr_in_air(o3_mmr, at_stp=True)
    >>> print("{} {} {:.3f}".format(o3_ppb.name(), o3_ppb.units,
    ...                             o3_ppb.data.max()))
    mole_fraction_of_ozone_in_air ppb 56.937

    Also supports conversion from a mass concentration.
    This example confirms that while it can invert cube_to_mass_conc_in_air,
    such use naturally introduces imprecision and should be avoided.

    >>> o3 = cube_to_mass_conc_in_air(o3_mmr, at_stp=True)
    >>> o3_ppb_2 = cube_to_vmr_in_air(o3, at_stp=True)
    >>> print("{} {} {:.3f}".format(o3.name(), o3.units,
    ...                             o3.data.max()))
    mass_concentration_of_ozone_in_air ug/m3 113.642
    >>> print("{:.3f} {:.3f} {:.3e}".format(
    ...     o3_ppb.data.max(),
    ...     o3_ppb_2.data.max(),
    ...     abs(o3_ppb.data-o3_ppb_2.data).max()))
    56.937 56.937 1.526e-05
    '''

    # There are a limited number of phenomena types we are
    # able to transform to VMR. Test each of these in turn

    cubename = cube.name()

    cubetype, molecule = cube_type(cubename)

    if cubetype == 'mole_fraction':
        # just copy cube, no transformation needed
        vmr_cube = cube.copy()

    # transform from mass mixing ratio
    elif cubetype == 'mass_fraction':
        vmr_cube = cube.copy()

        # ensure we know what units we're starting with
        vmr_cube.convert_units('kg/kg')

        # now sort out the data.
        # vmr = mmr * (molmass_air / molmass)
        try:
            molmass = chemistry.get_molmass(molecule)
        except ValueError:
            raise ValueError("cannot convert non-chemical species "
                             + molecule + " to volume mixing ratio")
        vmr_cube.data = vmr_cube.data*MMASS_AIR/molmass

        # set the new units
        vmr_cube.units = "mole/mole"

        # rename the cube
        vmr_cube.rename("mole_fraction_of_"+molecule+"_in_air")

    # convert from mass concentration
    elif cubetype == 'mass_concentration':
        vmr_cube = cube.copy()

        # ensure we know what units we're starting with
        vmr_cube.convert_units('kg/m3')

        # vmr = mc (R * temperature) / (pressure * molmass)
        # So we need molecular mass of the cube (and in kg/mol not g/mol)
        try:
            molmass = chemistry.get_molmass(molecule)/1000
        except ValueError:
            raise ValueError("cannot convert non-chemical species "
                             + molecule + " to volume mixing ratio")

        # Two options here - use standard temperature and pressure or
        # use cubes of pressure and temperature
        if at_stp:
            vmr_cube.data = vmr_cube.data * ((R * STD_T) / (STD_P * molmass))
        else:
            _ensure_p_t_compatible(pressure, temperature, vmr_cube)

            vmr_cube.data = vmr_cube.data * (
                (R * temperature.data) / (pressure.data * molmass))

        # set the new units
        vmr_cube.units = "mole/mole"

        # rename the cube
        vmr_cube.rename("mole_fraction_of_"+molecule+"_in_air")

    # if we get here then we can't convert to VMR
    else:
        raise ValueError("cannot convert "
                         + cubename + " to volume mixing ratio")

    # set units to those requested by ounits - default is ppb
    try:
        vmr_cube.convert_units(ounits)
    except:
        raise ValueError('cannot convert to '+ounits)

    return vmr_cube

def cube_to_molecules_conc_in_air(cube, ounits='cm-3',
                                  at_stp=False, pressure=None, temperature=None):
    """
    Convert to number of molecules concentration ('nmc') in air.
    By default converts to (molecules) cm-3, but can convert to any number
    concentration.

    :param cube: a cube containing the data to be converted.
    :param ounits: output units, defaulting to (molecules) cm-3
    :param at_stp: whether or not to use standard temperature and pressure.
                   Only applies to conversions from mmr and vmr.
                   Default False.
    :param pressure: cube containing pressure data, needed for conversions from
                     mmr and vmr when at_stp is set to False (otherwise
                     ignored).
    :param temperature: cube containing temperature data, needed for conversions
                        from mmr and vmr when at_stp is set to False
                        (otherwise ignored).
    :returns: nmc_cube, containing the converted data.


    Read in some test data:

    >>> import config
    >>> import pp_data
    >>> sample_data_path = config.SAMPLE_DATADIR+'aqum_output/oper/'
    >>> ppfiles = [sample_data_path+'aqum_20140325_T+000_QM18.pp']
    >>> pp = pp_data.PPData()
    >>> gcl = pp.readdata(ppfiles, short_name_list=['O3'])
    >>> o3_mmr = gcl[0]

    1. Test conversion from mmr to nmc at stp

    >>> print("{} {} {:.3e}".format(o3_mmr.name(), o3_mmr.units,
    ... o3_mmr.data.max()))
    mass_fraction_of_ozone_in_air kg kg-1 9.434e-08
    >>> o3_nmc = cube_to_molecules_conc_in_air(o3_mmr, ounits='cm-3',
    ... at_stp=True)
    >>> print("{} {} {:.3e}".format(o3_nmc.name(), o3_nmc.units,
    ... o3_nmc.data.max()))
    number_concentration_of_ozone_molecules_in_air cm-3 1.426e+12

    2. Test conversion from vmr to nmc at stp

    >>> o3_vmr = cube_to_vmr_in_air(o3_mmr, at_stp=True)
    >>> print("{} {} {:.3e}".format(o3_vmr.name(), o3_vmr.units,
    ... o3_vmr.data.max()))
    mole_fraction_of_ozone_in_air ppb 5.694e+01
    >>> o3_nmc = cube_to_molecules_conc_in_air(o3_vmr, ounits='cm-3',
    ... at_stp=True)
    >>> print("{} {} {:.3e}".format(o3_nmc.name(), o3_nmc.units,
    ... o3_nmc.data.max()))
    number_concentration_of_ozone_molecules_in_air cm-3 1.426e+12

    3. Test conversion from mc to nmc at stp

    >>> o3_mc = cube_to_mass_conc_in_air(o3_mmr, at_stp=True)
    >>> print("{} {} {:.3e}".format(o3_mc.name(), o3_mc.units,
    ... o3_mc.data.max()))
    mass_concentration_of_ozone_in_air ug/m3 1.136e+02
    >>> o3_nmc = cube_to_molecules_conc_in_air(o3_mc, ounits='cm-3',
    ... at_stp=True)
    >>> print("{} {} {:.3e}".format(o3_nmc.name(), o3_nmc.units,
    ... o3_nmc.data.max()))
    number_concentration_of_ozone_molecules_in_air cm-3 1.426e+12

    4. Test conversion from mmr to nmc using input temperature and pressure

    >>> import datetime
    >>> start_dt = datetime.datetime(2014, 8, 31, 0)
    >>> end_dt = datetime.datetime(2014, 9, 4, 0)
    >>> sample_data_path = config.SAMPLE_DATADIR+'aqum_output/oper/'
    >>> ppfiles = [sample_data_path+'prodm_op_aqum_20140901_18.000.pp']
    >>> pp = pp_data.PPData ()
    >>> gcl = pp.readdata(ppfiles, start_datetime=start_dt,
    ...     end_datetime=end_dt, short_name_list=['p','T','O3'])
    >>> pressure = gcl[0]
    >>> temperature = gcl[1]
    >>> o3_mmr = gcl[2]
    >>> print("{} {} {:.3e}".format(o3_mmr.name(), o3_mmr.units,
    ... o3_mmr.data.max()))
    mass_fraction_of_ozone_in_air kg kg-1 1.272e-05
    >>> o3_nmc = cube_to_molecules_conc_in_air(o3_mmr, ounits='cm-3',
    ... pressure=pressure, temperature=temperature)
    >>> print("{} {} {:.3e}".format(o3_nmc.name(), o3_nmc.units,
    ... o3_nmc.data.max()))
    number_concentration_of_ozone_molecules_in_air cm-3 5.548e+12

    """

    # There are a limited number of phenomena types we are
    # able to transform to number concentration in air. Test
    # each of these in turn.

    cubename = cube.name()
    cubetype, molecule = cube_type(cubename)

    # first check if it is already mass concentration in air
    if cubetype == 'number_concentration':
        nmc_cube = cube.copy()

    elif cubetype == 'mass_fraction':
        nmc_cube = cube.copy()

        # convert mass mixing ratio to kg/kg
        nmc_cube.convert_units('kg/kg')

        # Get molar mass
        try:
            molmass = chemistry.get_molmass(molecule)
        except ValueError:
            raise ValueError("cannot convert non-chemical species "
                             + molecule + " to molecules concentration")

        # mmr to molecules/cm3:
        # N (cm^-3) = (10^-6) * mmr * (m_air/m_species) * (p/kT)
        if at_stp:
            nmc_cube.data = (1e-6 * nmc_cube.data * (MMASS_AIR/molmass) *
                             STD_P / (Boltzmann * STD_T))
        else:
            _ensure_p_t_compatible(pressure, temperature, nmc_cube)

            nmc_cube.data = (1e-6 * nmc_cube.data * (MMASS_AIR/molmass) *
                             pressure.data / (Boltzmann * temperature.data))
        #Fix units
        nmc_cube.units = 'cm-3'

        #Set new name for cube
        nmc_cube.rename('number_concentration_of_'+molecule+'_molecules_in_air')

    elif cubetype == 'mole_fraction':
        nmc_cube = cube.copy()

        # Convert volume mixing ratio to ppb
        nmc_cube.convert_units('ppb')

        # ppb to molecules/cm3
        # N(cm^-3) = (10^-15) * ppb * (p/kT)
        if at_stp:
            nmc_cube.data = (1e-15 * nmc_cube.data *
                             STD_P / (Boltzmann * STD_T))

        else:
            _ensure_p_t_compatible(pressure, temperature, nmc_cube)

            nmc_cube.data = (1e-15 * nmc_cube.data *
                             pressure.data / (Boltzmann * temperature.data))

        #Fix units
        nmc_cube.units = 'cm-3'

        #Set new name for cube
        nmc_cube.rename('number_concentration_of_'+molecule+
                        '_molecules_in_air')

    elif cubetype == 'mass_concentration':
        nmc_cube = cube.copy()

        # Convert volume mixing ratio to ug/m3
        nmc_cube.convert_units('ug/m3')

        # Get molar mass
        try:
            molmass = chemistry.get_molmass(molecule)
        except ValueError:
            raise ValueError("cannot convert non-chemical species "
                             + molecule + " to molecules concentration")

        # ug/m3 to molecules/cm3
        # N(cm^-3) = (10^-12).mc.Avogadro/(M_species )
        nmc_cube.data = (1e-12 * nmc_cube.data *
                         Avogadro / molmass)

        #Fix units
        nmc_cube.units = 'cm-3'

        #Set new name for cube
        nmc_cube.rename('number_concentration_of_'+molecule+
                        '_molecules_in_air')


    # if we get here then we can't convert to number molecule conc
    else:
        raise ValueError(
            'Error in cube_to_molecule_conc_in_air. Cannot transform: '
            + cubename)

    # set units to those requested by ounits - default is cm-3
    try:
        nmc_cube.convert_units(ounits)
    except:
        raise ValueError('cannot convert to '+ounits)

    return nmc_cube


def cube_type(cubename):
    '''
    Test the phenomenon in the cube based on the name

    :param cubename: name of cube

    :returns: (cubetype, molecule), where molecule is the molecule name
              and cubetype has the options of:

              * 'mass_concentration' (eg typically units of ug/m3)
              * 'mass_fraction' (mass mixing ratio, mmr, eg typically kg/kg)
              * 'mole_fraction' (volume mixing ration, vmr, eg typically ppb)
              * 'number_concentration' (eg typically units of molecules/cm3)
              * None - if no matches (molecule is also set to None in this case)

    >>> cubetype, species = cube_type('mass_concentration_of_ozone_in_air')
    >>> print(cubetype, species)
    mass_concentration ozone

    >>> cubetype, species = cube_type('mass_fraction_of_nitrogen_dioxide_in_air')
    >>> print(cubetype, species)
    mass_fraction nitrogen_dioxide

    >>> cubetype, species = cube_type('mole_fraction_of_pm2p5_dry_aerosol_in_air')
    >>> print(cubetype, species)
    mole_fraction pm2p5_dry_aerosol

    >>> cubetype, species = cube_type('number_concentration_of_ozone_molecules_in_air')
    >>> print(cubetype, species)
    number_concentration ozone

    If there are no matches, then returns None, None:

    >>> cubetype, species = cube_type('air_temperature')
    >>> print(cubetype, species)
    None None

    '''

    # first set up the matches for all the types of
    # phenomena we can transform so that we can
    # use them in the if tests
    match = re.search(
        r"^(mass_concentration|mass_fraction|mole_fraction)_of_(.+)_in_air$",
        cubename)

    if not match:
        #Try number of molecules concentation
        match = re.search(r"^(number_concentration)_of_(.+)_molecules_in_air$",
                          cubename)

    if match:
        cubetype = match.group(1)
        molecule = match.group(2)
    else:
        cubetype = None
        molecule = None

    return cubetype, molecule

def _ensure_p_t_compatible(pressure, temperature, cube=None):
    """
    Helper function to verify that the provided pressure and temperature cubes
    are as expected. That is:

    - not None
    - appropriate names of air_temperature and air_pressure
    - standard units Pa and K, converting if necessary
    - same dimensions as each other
    - (optionally) same dimensions as a given data cube

    :param pressure: a cube to check for pressure data
    :param temperature: a cube to check for temperature data
    :param cube: a cube to compare dimensions with the pressure and temperature
                 cubes
    :returns: None, although pressure and temperature may have different units

    Load some example pressure and temperature data:

    >>> import adaq_functions
    >>> ini_dict, sites_data, od, md_list = adaq_functions.get_exampledata(
    ...     exampletype="3d") # doctest: +ELLIPSIS
    Reading inifile .../example_data_3d.ini
    Number of sites:  5
    >>> pressure = md_list[0].extract(short_name="p", singlecube=True)
    >>> temperature = md_list[0].extract(short_name="T", singlecube=True)
    >>> o3 = md_list[0].extract(short_name="O3", singlecube=True)

    Confirm that everything is compatible:

    >>> _ensure_p_t_compatible(pressure, temperature, o3)

    More importantly, check unexpected inputs are handled correctly:

    >>> _ensure_p_t_compatible(None, None)
    Traceback (most recent call last):
      File ".../cube_chemistry.py", line 556, in _ensure_p_t_compatible
        raise ValueError('did not pass in a pressure '
    ValueError: did not pass in a pressure cube named air_pressure

    >>> _ensure_p_t_compatible(pressure=temperature, temperature=pressure)
    Traceback (most recent call last):
      File ".../cube_chemistry.py", line 556, in _ensure_p_t_compatible
        raise ValueError('did not pass in a pressure '
    ValueError: did not pass in a pressure cube named air_pressure

    >>> temp = temperature.extract(iris.Constraint(model_level_number=[1,10]))
    >>> _ensure_p_t_compatible(pressure, temp)
    Traceback (most recent call last):
      File ".../cube_chemistry.py", line 568 in _ensure_p_t_compatible
        raise ValueError('Temperature and pressure not co-located')
    ValueError: Temperature and pressure not co-located

    >>> pressure.convert_units("atm")
    >>> print(pressure.units)
    atm
    >>> _ensure_p_t_compatible(pressure, temperature)
    >>> print(pressure.units)
    Pa
    """

    # check that pressure is really pressure
    if pressure is None or pressure.name() != 'air_pressure':
        raise ValueError('did not pass in a pressure '
                         'cube named air_pressure')

    # check that temperature is really temperature
    if temperature is None or temperature.name() != 'air_temperature':
        raise ValueError('did not pass in a temperature cube '
                         'named air_temperature')

    # check if cubes co-located
    coord_comparison = iris.analysis.coord_comparison(temperature, pressure)
    if (coord_comparison['not_equal'] or
            coord_comparison['non_equal_data_dimension']):
        raise ValueError('Temperature and pressure not co-located')

    # pressure in Pa and T in K using convert_units
    try:
        pressure.convert_units('Pa')
    except:
        raise ValueError('could not convert pressure units to Pa')
    try:
        temperature.convert_units('K')
    except:
        raise ValueError('could not convert temperature units to K')

    # check that p and T co-located with input cube
    # only need to test one as p and T already co-located
    if cube is not None:
        coord_comparison = iris.analysis.coord_comparison(temperature, cube)
        if (coord_comparison['not_equal'] or
                coord_comparison['non_equal_data_dimension']):
            raise ValueError('T/p and data not co-located')


def calculate_nox(cubelist):
    """
    Calculate NOx = NO + NO2, as a volume mixing ratio in ppb.
    Prefers input data in ppb, but will attempt to convert if necessary.

    Note that when converting from a mass concentration, standard temperature
    and pressure is assumed. If this assumption is not appropriate, ensure
    the data are converted before calling this function.

    :param cubelist: a Cubelist containing NO and NO2 data
    :returns: the same cubelist, which will contain the calculated NOx cube

    >>> import config
    >>> import pp_data
    >>> sample_data_path = config.SAMPLE_DATADIR+'aqum_output/oper/'
    >>> pp = pp_data.PPData()
    >>> gcl = pp.readdata([sample_data_path+'aqum_20140330_T+000_QM18.pp'],
    ...     short_name_list=['NO2', 'NO'])
    >>> no = pp.extract(short_name='NO', gridded=True, singlecube=True)
    >>> no2 = pp.extract(short_name='NO2', gridded=True, singlecube=True)
    >>> print('{:.3e} {} {:.3e} {}'.format(no.data.max(), no.units,
    ...                                    no2.data.max(), no2.units))
    1.604e-07 kg kg-1 9.028e-08 kg kg-1
    >>> gcl = calculate_nox(gcl)
    >>> nox = pp.extract(short_name='NOx', gridded=True, singlecube=True)
    >>> print('{:.3f} {}'.format(nox.data.max(), nox.units))
    199.253 ppb
    """

    #Extract required components
    try:
        no = cubelist.extract(iris.AttributeConstraint(short_name='NO'),
                              strict=True)
    except iris.exceptions.ConstraintMismatchError:
        warnings.warn('NOx not calculated: NO missing')
        return cubelist
    try:
        no2 = cubelist.extract(iris.AttributeConstraint(short_name='NO2'),
                               strict=True)
    except iris.exceptions.ConstraintMismatchError:
        warnings.warn('NOx not calculated: NO2 missing')
        return cubelist

    #Convert to correct units
    no = cube_to_vmr_in_air(no, "ppb", at_stp=True)
    no2 = cube_to_vmr_in_air(no2, "ppb", at_stp=True)

    #Sum to get NOx and tidy up name and units.
    nox = no + no2
    nox.rename('mole_fraction_of_nox_in_air')
    nox.units = 'ppb'

    #attempt to intersect the attributes of the constituent cubes
    for key, val in no.attributes.items():
        if key in no2.attributes and no2.attributes[key] == val:
            nox.attributes[key] = val
    nox.attributes['short_name'] = 'NOx'

    cubelist.append(nox)
    return cubelist

def calculate_ox(cubelist):
    """
    Calculate Ox = O3 + NO2, as a volume mixing ratio in ppb.
    Prefers input data in ppb, but will attempt to convert if necessary.

    Note that when converting from a mass concentration, standard temperature
    and pressure is assumed. If this assumption is not appropriate, ensure
    the data are converted before calling this function.

    :param cubelist: a Cubelist containing O3 and NO2 data
    :returns: the same cubelist, which will contain the calculated Ox cube

    >>> import config
    >>> import pp_data
    >>> sample_data_path = config.SAMPLE_DATADIR+'aqum_output/oper/'
    >>> pp = pp_data.PPData()
    >>> gcl = pp.readdata([sample_data_path+'aqum_20140330_T+000_QM18.pp'],
    ...     short_name_list=['NO2','O3'])
    >>> o3 = pp.extract(short_name='O3', gridded=True, singlecube=True)
    >>> no2 = pp.extract(short_name='NO2', gridded=True, singlecube=True)
    >>> print('{:.3e} {} {:.3e} {}'.format(o3.data.max(), o3.units,
    ...                                    no2.data.max(), no2.units))
    1.127e-07 kg kg-1 9.028e-08 kg kg-1
    >>> gcl = calculate_ox(gcl)
    >>> ox = pp.extract(short_name='Ox', gridded=True, singlecube=True)
    >>> print('{:.3f} {}'.format(ox.data.max(), ox.units))
    68.413 ppb
    """

    #Extract required components
    try:
        o3 = cubelist.extract(iris.AttributeConstraint(short_name='O3'),
                              strict=True)
    except iris.exceptions.ConstraintMismatchError:
        warnings.warn('Ox not calculated: O3 missing')
        return cubelist
    try:
        no2 = cubelist.extract(iris.AttributeConstraint(short_name='NO2'),
                               strict=True)
    except iris.exceptions.ConstraintMismatchError:
        warnings.warn('Ox not calculated: NO2 missing')
        return cubelist

    #Convert to correct units
    o3 = cube_to_vmr_in_air(o3, "ppb", at_stp=True)
    no2 = cube_to_vmr_in_air(no2, "ppb", at_stp=True)

    #Sum to get Ox and tidy up name and units.
    ox = o3 + no2
    ox.rename('mole_fraction_of_ox_in_air')
    ox.units = 'ppb'

    #attempt to intersect the attributes of the constituent cubes
    for key, val in o3.attributes.items():
        if key in no2.attributes and no2.attributes[key] == val:
            ox.attributes[key] = val
    ox.attributes['short_name'] = 'Ox'

    cubelist.append(ox)
    return cubelist

def calculate_nh4_aerosol(cubelist, pm='PM10'):
    """
    Calculate mass concentration of NH4 from NH4NO3 and NH42SO4
    Assumes input in concentration units

    :param cubelist: Iris cubelist, containing cube of mass concentration of
                     NH4NO3 (ammonium nitrate) and cube of mass concentration of
                     (NH4)2SO4 (ammonium sulfate)
    :param pm: String, size of PM, eg 'PM10' if using PM10_NH4NO3 and PM10_NH42SO4

    :returns: Iris cubelist as input, but with extra iris cube of mass
              concentration of NH4 (ammonium) aerosol

    >>> import config
    >>> import pp_data
    >>> sample_data_path = config.SAMPLE_DATADIR+'aqum_output/oper/'
    >>> pp = pp_data.PPData()
    >>> gcl = pp.readdata([sample_data_path+'aqum_20140330_T+000_QM18.pp'],
    ... short_name_list=['PM10_NH4NO3','PM10_NH42SO4'])
    >>> nh4no3 = pp.extract(short_name='PM10_NH4NO3', singlecube=True,
    ... gridded=True)
    >>> nh42so4 = pp.extract(short_name='PM10_NH42SO4', singlecube=True,
    ... gridded=True)
    >>> print('{:.3f} {} {:.3f} {}'.format(nh4no3.data.max(), nh4no3.units, \
nh42so4.data.max(), nh42so4.units))
    47.922 ug m-3 31.480 ug m-3
    >>> gcl = calculate_nh4_aerosol(gcl,pm='PM10')
    >>> nh4 = gcl.extract(iris.AttributeConstraint(short_name='PM10_NH4'),
    ... strict=True)
    >>> print('{:.3f} {}'.format(nh4.data.max(), nh4.units))
    12.227 ug m-3
    """

    nh4no3 = cubelist.extract(iris.AttributeConstraint(short_name=pm+'_NH4NO3'),
                              strict=True)
    nh42so4 = cubelist.extract(iris.AttributeConstraint(short_name=pm+'_NH42SO4'),
                               strict=True)

    #Set up NH4 components of NH4NO3 and NH42SO4
    nh4_nh4no3 = None
    nh4_nh42so4 = None

    mmass_nh4 = chemistry.get_molmass('ammonium')


    if nh4no3 is not None:
        pattern = r"^(mass_concentration_of_ammonium_nitrate_in_pm)(.+)(_dry_aerosol_in_air)$"
        match = re.search(pattern, nh4no3.name())
        if not match:
            raise ValueError('Ammonium nitrate should be in mass concentration units')
        mmass_nh4no3 = chemistry.get_molmass('ammonium nitrate')
        nh4_nh4no3 = nh4no3 * mmass_nh4 / mmass_nh4no3

    if nh42so4 is not None:
        pattern = r"^(mass_concentration_of_ammonium_sulfate_in_pm)(.+)(_dry_aerosol_in_air)$"
        match = re.search(pattern, nh42so4.name())
        if not match:
            raise ValueError('Ammonium sulphate should be in mass concentration units')
        mmass_nh42so4 = chemistry.get_molmass('ammonium sulfate')
        nh4_nh42so4 = nh42so4 * 2. * mmass_nh4 / mmass_nh42so4

    if nh4no3 is not None:
        nh4 = nh4_nh4no3
        if nh42so4 is not None:
            if nh4no3.units != nh42so4.units:
                raise ValueError('Units are not the same for NH4NO3 and NH42SO4')
            nh4 += nh4_nh42so4
        nh4.rename(nh4no3.name().replace('ammonium_nitrate', 'ammonium'))
        nh4.convert_units(nh4no3.units) #Ensure in the same unit format as input data
        nh4.attributes['label'] = nh4no3.attributes['label']
    elif nh42so4 is not None:
        nh4 = nh4_nh42so4
        nh4.rename(nh42so4.name().replace('ammonium_sulfate', 'ammonium'))
        nh4.convert_units(nh42so4.units) #Ensure in the same unit format as input data
        nh4.attributes['label'] = nh42so4.attributes['label']
    else:
        warnings.warn('Both input data of NH4NO3 and NH42SO4 are None: NH4 cannot be calculated')
        nh4 = None


    if nh4 is not None:
        nh4.attributes['short_name'] = pm+'_NH4'
        cubelist.append(nh4)

    return cubelist


def calculate_no3_aerosol(cubelist, pm='PM10'):
    '''
    Calculate mass concentration of NO3 aerosol from NH4NO3 aerosol

    :param cubelist: Iris cubelist containing cube of mass concentration
                     of NH4NO3 (ammonium nitrate)
    :param pm: String, size of PM, eg 'PM10' if using PM10_NH4NO3, or
               'PM2p5' if using 'PM2p5_NH4NO3'
    :returns: Iris cubelist as input, but also containing cube of mass
              concentration of NO3 (nitrate) aerosol

    >>> import config
    >>> import pp_data
    >>> sample_data_path = config.SAMPLE_DATADIR+'aqum_output/oper/'
    >>> pp = pp_data.PPData()
    >>> gcl = pp.readdata([sample_data_path+'aqum_20140330_T+000_QM18.pp'],
    ... short_name_list=['PM10_NH4NO3'])
    >>> nh4no3 = pp.extract(short_name='PM10_NH4NO3', singlecube=True,
    ... gridded=True)
    >>> print('{:.3f} {}'.format(nh4no3.data.max(), nh4no3.units))
    47.922 ug m-3
    >>> gcl = calculate_no3_aerosol(gcl, pm='PM10')
    >>> no3 = gcl.extract(iris.AttributeConstraint(short_name='PM10_NO3'),
    ... strict=True)
    >>> print('{:.3f} {}'.format(no3.data.max(), no3.units))
    37.122 ug m-3
    '''

    nh4no3 = cubelist.extract(iris.AttributeConstraint(short_name=pm+'_NH4NO3'),
                              strict=True)

    if not nh4no3:
        warnings.warn(pm+'_NO3 cannot be calculated as '+pm+
                      '_NH4NO3 not in input cubelist')
        return cubelist

    pattern = r"^(mass_concentration_of_ammonium_nitrate_in_pm)(.+)(_dry_aerosol_in_air)$"
    match = re.search(pattern, nh4no3.name())
    if not match:
        raise ValueError('Ammonium nitrate should be in mass concentration units')

    mmass_no3 = chemistry.get_molmass('nitrate')
    mmass_nh4no3 = chemistry.get_molmass('ammonium nitrate')

    #Calculate NO3
    no3 = nh4no3 * mmass_no3 / mmass_nh4no3
    no3.rename(nh4no3.name().replace('ammonium_nitrate', 'nitrate'))
    no3.convert_units(nh4no3.units)
    no3.attributes['short_name'] = pm+'_NO3'
    no3.attributes['label'] = nh4no3.attributes['label']

    cubelist.append(no3)

    return cubelist


def calculate_so4_aerosol(cubelist, pm='PM10'):
    '''
    Calculate mass concentration of SO4 aerosol from (NH4)2SO4 aerosol

    :param cubelist: Iris cubelist containing cube of mass concentration
                     of (NH4)2SO4 (ammonium sulfate)
    :param pm: String, size of PM, eg 'PM10' if using PM10_NH42SO4, or
               'PM2p5' if using 'PM2p5_NH42SO4'
    :returns: Iris cubelist as input, but also containing cube of mass
              concentration of SO4 (sulfate) aerosol

    >>> import config
    >>> import pp_data
    >>> sample_data_path = config.SAMPLE_DATADIR+'aqum_output/oper/'
    >>> pp = pp_data.PPData()
    >>> gcl = pp.readdata([sample_data_path+'aqum_20140330_T+000_QM18.pp'],
    ... short_name_list=['PM10_NH42SO4'])
    >>> nh42so4 = pp.extract(short_name='PM10_NH42SO4', singlecube=True,
    ... gridded=True)
    >>> print('{:.3f} {}'.format(nh42so4.data.max(), nh42so4.units))
    31.480 ug m-3
    >>> gcl = calculate_so4_aerosol(gcl, pm='PM10')
    >>> so4 = gcl.extract(iris.AttributeConstraint(short_name='PM10_SO4'),
    ... strict=True)
    >>> print('{:.3f} {}'.format(so4.data.max(), so4.units))
    22.885 ug m-3
    '''

    nh42so4 = cubelist.extract(iris.AttributeConstraint(short_name=pm+'_NH42SO4'),
                               strict=True)
    if not nh42so4:
        warnings.warn(pm+'_SO4 cannot be calculated as '+pm+
                      '_NH42SO4 not in input cubelist')
        return cubelist

    pattern = r"^(mass_concentration_of_ammonium_sulfate_in_pm)(.+)(_dry_aerosol_in_air)$"
    match = re.search(pattern, nh42so4.name())
    if not match:
        raise ValueError('Ammonium sulfate should be in mass concentration units')

    mmass_so4 = chemistry.get_molmass('sulfate anion')
    mmass_nh42so4 = chemistry.get_molmass('ammonium sulfate')

    #Calculate SO4
    so4 = nh42so4 * mmass_so4 / mmass_nh42so4
    so4.rename(nh42so4.name().replace('ammonium_sulfate', 'sulfate'))
    so4.convert_units(nh42so4.units)
    so4.attributes['short_name'] = pm+'_SO4'
    so4.attributes['label'] = nh42so4.attributes['label']

    cubelist.append(so4)

    return cubelist

def calculate_total_voc(cubelist):
    """
    Calculate total VOC concentration from VOCs used in model

    :param cubelist: Iris cubelist containing cube of mass concentration
                     of VOCs
    :returns: Iris cubelist as input, but also containing cube of total VOC
              in units of ppb

    >>> import config
    >>> import pp_data
    >>> sample_data_path = config.SAMPLE_DATADIR+'aqum_output/'
    >>> pp = pp_data.PPData()
    >>> gcl = pp.readdata([sample_data_path+'prods_op_aqum_20191001_18.000.pp'],
    ...    short_name_list=['HCHO','C2H6','MeCHO','C3H8','Me2CO','C5H8','C3H6',
    ...                     'TOLUENE','CH3OH','C4H10','XYLENE','C2H4'])
    >>> HCHO = pp.extract(short_name='HCHO', gridded=True, singlecube=True)
    >>> C2H4 = pp.extract(short_name='C2H4', gridded=True, singlecube=True)
    >>> print('{:3e} {} {:3e} {}'.format(HCHO.data.mean(), HCHO.units,
    ...                                  C2H4.data.mean(), C2H4.units))
    3.427712e-10 kg kg-1 1.845376e-10 kg kg-1
    >>> gcl = calculate_total_voc(gcl)
    >>> total_voc = pp.extract(short_name='Total_VOC', gridded=True, singlecube=True)
    >>> print('{:.3f} {}'.format(total_voc.data.mean(), total_voc.units))
    2.192 ppb
    """

    voc_species = ['HCHO', 'C2H6', 'MeCHO', 'C3H8',
                   'Me2CO', 'C5H8', 'C3H6', 'TOLUENE', 'CH3OH',
                   'C4H10', 'XYLENE', 'C2H4']

    total_voc = 0
    for species in voc_species: 
        voc = cubelist.extract(iris.AttributeConstraint(short_name=species),
                              strict=True)
        voc_vmr = cube_to_vmr_in_air(voc, 'ppb', at_stp=True)
        total_voc = total_voc + voc_vmr


    total_voc.rename('mass_fraction_of_total_voc_in_air')
    total_voc.attributes['short_name'] = 'Total_VOC'
    total_voc.attributes['label'] = cubelist[0].attributes['label']
    total_voc.attributes['cubetype'] = 'mole_fraction'
    total_voc.units = 'ppb'
    cubelist.append(total_voc)
    return cubelist

def calc_dobson_units(md_list, p_target_top, short_name):

    '''
    This function takes a model data list containing 3D data cubes for
    temperature, pressure and species concentration. Using this input,
    it determines the Dobson units within a user specified vertical
    range that is defined by p_target_top (upper limit). The p_target_top
    value corresponds to the satellite's upper retrieval limit. The lower
    vertical limit currently defaults to the surface. A p_target_bottom
    variable should be added at some point so that the satellite's
    lower retrieval limit can be captured. This would require additional
    code in the __calc_boundary_box function to determine the dimensions
    of the 'new' bottom grid cell.

    This function returns a cube containing gridded Dobson unit data and
    gridded concentration data for each species. These arrays can be
    used for diagnostic plotting in the calling routine.
    (e.g. mmr2dobsonunits.py).

    :param md_list: the model data for temperature, pressure and species \
    concentration. Note that if Dobson units are being calculated for more \
    than one species, the pressure and temperature columns are partially \
    overwritten with NaN values on the first iteration of the code. The \
    species concentration column is also partially overwritten. Required.
    :param p_target_top: the upper pressure level in Pa. Dobson units are \
    calculated below this level. This value corresponds to the upper limit \
    in the atmosphere that the satellite can retrieve to. Required.
    :param short_name: Short name of the chemical species to calculate Dobson \
    units for. Required.

    :returns: du_cube, sp_cube

    * du_cube is a cube containing gridded Dobson units data for the required species.
    * sp_cube is a cube containing gridded concentration data for required species.

    Note that if aqcops data is being used calc_dobson_units can fail in
    __calc_boundary_box due to a memory error if more than one satellite
    pass time has been loaded and a long time series (approx. one month
    or greater) is being processed.

    Examples:

    Load some data. 3D model data for temperature, pressure and species
    concentration is loaded from the prodm file. A prods file for
    orography is also loaded so that the 3D altitude data can be determined.

    >>> import config
    >>> import pp_data
    >>> import datetime
    >>> start_dt = datetime.datetime(2017, 7, 2, 6)
    >>> end_dt = datetime.datetime(2017, 7, 2, 23)
    >>> p_target_top = 45000.0
    >>> sample_data_path = config.SAMPLE_DATADIR + '/aqum_output/oper/3d'
    >>> ppfiles = [sample_data_path+'/prodm*2017070[1-2]*.012.pp',
    ... sample_data_path+'/prods_op_aqum_20170701_18.000.pp']
    >>> md_list = pp_data.PPData()
    >>> md_list.readdata(ppfiles, start_datetime=start_dt,
    ...     end_datetime=end_dt, short_name_list=['p','T','O3'], forecast_day=1)
    [<iris 'Cube' of air_pressure / \
(Pa) (time: 2; model_level_number: 63; grid_latitude: 182; grid_longitude: 146)>,
    <iris 'Cube' of air_temperature / \
(K) (time: 2; model_level_number: 63; grid_latitude: 182; grid_longitude: 146)>,
    <iris 'Cube' of mass_fraction_of_ozone_in_air / \
(kg kg-1) (time: 2; model_level_number: 63; grid_latitude: 182; grid_longitude: 146)>]

    >>> du_cube, sp_cube = calc_dobson_units(md_list, p_target_top, short_name='O3')
    Extracting temperature cube
    Extracting pressure cube
    Extracting species cube
    Calculating Dobson units for species

    >>> print(du_cube) #doctest: +ELLIPSIS
    Tropospheric_ozone_column_abundance_in_Dobson_Units / \
(DU) (time: 2; grid_latitude: 182; grid_longitude: 146)
         Dimension coordinates:
              time                                                  \
x                 -                    -
              grid_latitude                                         \
-                 x                    -
              grid_longitude                                        \
-                 -                    x
         Auxiliary coordinates:
              surface_altitude                                      \
-                 x                    x
         Derived coordinates:
              altitude                                              \
-                 x                    x
         Scalar coordinates:
              forecast_day: 1.0 Days
              level_height: 21324.0... m, bound=(0.0, 42648...) m
              sigma: 0.5, bound=(0.0, 1.0)
         Attributes:
              STASH: m01s34i001
              label: pp
              short_name: O3_DU
              source: Data from Met Office Unified Model
         Cell methods:
              mean: time (1 hour)
              sum: model_level_number

    >>> print(du_cube.attributes['short_name'], du_cube.units)
    O3_DU DU

    >>> print('{:.2f} {:.2f}'.format(np.max(du_cube.data), np.min(du_cube.data)))
    28.28 10.94
    '''


    ## Check that altitude coordinate has been created and loaded
    try:
        md_list.gridded_cube_list[0].coord('altitude')
    except:
        raise ValueError('Altitude coordinate not available!')

    ## Extract the altitude data and determine model level box
    altitude_coord = md_list.gridded_cube_list[0].coord('altitude')
    orography = md_list.gridded_cube_list[0].coord('surface_altitude').points
    orography = iris.util.broadcast_to_shape(
        orography, altitude_coord.bounds.shape, (1, 2))

    box_hgt_bounds = altitude_coord.bounds

    altitudes_arr = altitude_coord.points

    ## Calculate height of each vertical box
    heights_arr = np.diff(box_hgt_bounds, n=1, axis=3)
    heights_arr = np.squeeze(heights_arr)


    ## Extract temperature cube

    print('Extracting temperature cube')

    tconst = iris.Constraint(name='air_temperature')
    t_cubelist = md_list.gridded_cube_list.extract(tconst)
    t_cube = t_cubelist[0]
    t_arr = t_cube.data

    ## Check that temperature data has loaded
    if t_cube.name() != 'air_temperature':
        raise ValueError('did not pass in a cube named air_temperature')

    ## Check temperature units
    assert t_cube.units == 'K', 'Temperature units not in Kelvin!'

    ## Check that dimensions for t_cube are ordered as expected.
    nx, ny, nz, nt = __check_dim_coords(t_cube)


    ## Extract pressure cube

    print('Extracting pressure cube')

    pconst = iris.Constraint(name='air_pressure')
    p_cubelist = md_list.gridded_cube_list.extract(pconst)
    p_cube = p_cubelist[0]
    p_arr = p_cube.data

    ## Check that pressure data has loaded
    if p_cube.name() != 'air_pressure':
        raise ValueError('did not pass in a cube named air_pressure')

    ## Check temperature units
    assert p_cube.units == 'Pa', 'Pressure units not in Pa!'

    ## Check that dimensions for p_cube are ordered as expected.
    nx, ny, nz, nt = __check_dim_coords(p_cube)


    ## Check if temperature and pressure cubes are co-located
    coord_comparison = iris.analysis.coord_comparison(t_cube,
                                                      p_cube)

    if (coord_comparison['not_equal'] or coord_comparison['non_equal_data_dimension']):
        raise ValueError('Temperature and pressure not co-located')


    ## Extract species cubes

    print('Extracting species cube')
    sp_cubelist = md_list.gridded_cube_list.extract(iris.AttributeConstraint(
        short_name=short_name))
    sp_cube = sp_cubelist[0]
    sp_arr = sp_cube.data

    ## Check that sp and t co-located with input cube
    ## Only need to test one as p and t already co-located
    coord_comparison = iris.analysis.coord_comparison(t_cube,
                                                      sp_cube)

    if (coord_comparison['not_equal'] or
            coord_comparison['non_equal_data_dimension']):
        raise ValueError('T/p and concentration not co-located')


    ## Check if species concentration unit have mass mixing ratio in kg kg-1
    ## If not, raise error (NB: will need to add conversions later...)
    if not sp_cube.units == 'kg kg-1':
        raise ValueError(
            sp_cube.name + 'does not have units kg kg-1!')


    ## Check that dimensions for sp_cube are ordered as expected.
    nx, ny, nz, nt = __check_dim_coords(sp_cube)


    ## Set up the heights array with a time dimension.
    if heights_arr.shape == (nz, ny, nx):
        heights_arr = heights_arr.reshape(1, nz, ny, nx).repeat(nt, 0)


    ## Calculate DU for surface -> 450 hPa and for 'top box'
    print('Calculating Dobson units for species')

    ## Loop through every column and set the data in grid
    ## cells above the required p_target_top to nan

    for itime in range(nt):
        for ilat in range(ny):
            for ilon in range(nx):

                ## Firstly check if top box has been identified. This is
                ## determined by checking for NaN values in the pressure
                ## data column. NaN values will be present if the code is
                ## looping over multiple species and it is not the first iteration
                nan_inds = np.argwhere(np.isnan(p_arr[itime, :, ilat, ilon]))

                ## If the code is not on the first iteration, the new top boundary box
                ## does not need to be re-identified, just modify the species column.
                ## If the code is on the first iteration, __calc_boundarybox must run
                if len(nan_inds) > 0:
                    sp_arr[itime, nan_inds, ilat, ilon] = np.nan

                else:
                    __calc_boundarybox(p_arr[itime, :, ilat, ilon],
                                       t_arr[itime, :, ilat, ilon],
                                       sp_arr[itime, :, ilat, ilon],
                                       box_hgt_bounds,
                                       altitudes_arr[:, ilat, ilon],
                                       heights_arr[itime, :, ilat, ilon],
                                       ilon, ilat,
                                       p_target_top)


    ## Calculate Dobson Units - step 1
    col = sp_arr * ((p_arr * heights_arr) / t_arr)
    col_total = np.nansum(col, axis=1)

    ## Get the molar mass of the species
    match = re.search(r"^mass_fraction_of_(.+)_in_air$", sp_cube.name())
    mmass_sp = chemistry.get_molmass(match.group(1))

    ## Calculate Dobson Units - step 2
    du_arr = 1.0E+5 * (273.0/STD_P) * (MMASS_AIR / mmass_sp) * col_total


    ## Set up a cube for Dobson units
    du_cube = sp_cube.copy()
    du_cube.data = du_cube.data * 0.0
    du_cube = du_cube.collapsed('model_level_number', iris.analysis.SUM)
    du_cube.remove_coord('model_level_number')
    du_cube.remove_coord('forecast_period')

    ## Set up new units and names in Dobson units cube
    dobson_units = cf_units.Unit('DU')
    du_cube.rename('Tropospheric_'+match.group(1)+'_column_abundance_in_Dobson_Units')
    du_short_name = du_cube.attributes['short_name'] + '_DU'
    du_cube.units = dobson_units
    du_cube.attributes['short_name'] = du_short_name


    ## Save Dobson unit data in Dobson unit cube
    du_cube.data = du_arr


    return du_cube, sp_cube


def __check_dim_coords(cube):

    """
    This function takes cubes of temperature, pressure, and species
    concentration from calc_dobson_units and checks that they have
    the correct dimensions in the expected order:
    (time, model level number, latitude, and longitude). It returns
    the length of each dimension. This function can only be called
    from calc_dobson_units.

    :param cube: a cube containing 4D model data. Required

    :returns nx, ny, nz, nt: where nx is the length of the X (longitude) \
    coordinate; ny is the length of the Y (latitude) coordinate; nz is \
    the length of the Z (model level number) coordinate and nt is the \
    length of the T (time) coordinate.
    """

    ## Set nx, ny, nz and nt to zero initially
    nx = 0
    ny = 0
    nz = 0
    nt = 0

    ## Check that more than two time points have been loaded.
    ## Code fails with less than two time points, possibly
    ## due to a bug in Iris.

    assert cube.coord_dims('time') != (), (
        'Time dimension currently needs to have at least 2 times!')

    ## Loop over the coordinates
    for coord in cube.dim_coords:

        if iris.util.guess_coord_axis(coord) == 'X':
            assert cube.coord_dims(coord) == (3,), (
                'Longitude is not correctly ordered!')
            nx = len(coord.points)

        if iris.util.guess_coord_axis(coord) == 'Y':
            assert cube.coord_dims(coord) == (2,), (
                'Latitude is not correctly ordered!')
            ny = len(coord.points)

        if iris.util.guess_coord_axis(coord) == 'Z':
            assert cube.coord_dims(coord) == (1,), (
                'Model level number is not correctly ordered!')
            nz = len(coord.points)

        if iris.util.guess_coord_axis(coord) == 'T':
            assert cube.coord_dims(coord) == (0,), (
                'Time is not correctly ordered!')
            nt = len(coord.points)


    ## Check that all dimensions are non-zero
    assert nx > 0, 'Longitude dimension has zero length!'
    assert ny > 0, 'Latitude dimension has zero length!'
    assert nz > 0, 'Model level number dimension has zero length!'
    assert nt > 0, 'Time dimension has zero length!'


    return nx, ny, nz, nt



def __calc_boundarybox(pvals, tvals, spvals,
                       box_hgt_bounds,
                       zvals, hvals,
                       ilon, ilat,
                       p_target_top):

    """
    This function locates the grid box in the vertical dimension that
    corresponds most closely to the user specified upper pressure value
    (p_target_top). It then calculates a new top box based on the
    location of the nearest vertical box and p_target_top. Data in
    grid boxes above the new top box is set to NaN.
    Note that functionality should be added so that the grid box
    closest to a lower pressure value (which would be determined by
    the satellite's lower retrieval limit) can also be determined.
    Currently the lower vertical limit is assumed to be the surface.
    This function returns 3D arrays of pressure, temperature, species
    concentration and grid box height modified to only include data
    from the surface to the new top grid box.
    This function can only be called from calc_dobson_units.

    :param pvals: the column of air pressure data at ilon, ilat, itime. Required.
    :param tvals: the column of air temperature data at ilon, ilat, itime. Required.
    :param spvals: the column of species concentration data at ilon, ilat, itime. Required.
    :param box_hgt_bounds: the original 4D array of box height bounds data. Required.
    :param zvals: the column of altitude data at ilon, ilat. Required.
    :param hvals: the column of grid box height data at ilon, ilat. Required.
    :param ilon: the counter for longitude values. Required.
    :param ilat: the counter for latitude values. Required.
    :param p_target_top: the user specified upper pressure value. Required.

    Please see the diagram below to help explain the variables that
    are used in __calc_boundary_box to calculate the new top box:

    'Boxes' are defined by rho points.

    -- rho points   |-------------|
                    |             |
                    |             |
    .. theta points |.............|
       (p,T,z,sp)   |             |
                    |             |
                    |-------------|   -> p_rho_idz_top
                    |             |
    ++ p target top |+++++++++++++| } -> upper bound => z_top_box_ub = idz
                    |.............| }
                    |      o      | }
                    |             | } -> top box from ++ to --
                    |-------------| } -> lower bound => z_top_box_lb
                    |             |
                    |             |
                    |.............|
                    |             |
                    |             |
                    |-------------|


    * idz = theta point nearest to p_target_top
    * top box is the highest 'box' (between rho levels) that includes p target top.
    * New top box is defined by p target top at the top point and the lower bound
      defined by the rho point below it.
    * o = centre of new top box => New 'theta' point, which through interpolation
      gives p_centre_new, and T for the new top box.
    """

    ## Locate the array index of pressure value closest to p_min
    idz = (np.abs(pvals - p_target_top)).argmin()

    ## Check that p vals and t vals are decreasing with altitude and
    ## that z vals are increasing with altitude.
    ## Mask any 'spurious' values.
    for val in range(len(pvals)-1):
        if pvals[val] < pvals[val+1]:
            pvals[val] = np.nan

        if tvals[val] < tvals[val+1]:
            tvals[val] = np.nan

        if zvals[val] > zvals[val+1]:
            zvals[val] = np.nan


    ## Interpolate p at the upper bound of box idz
    p_rho_idz_top = np.interp(box_hgt_bounds[idz, ilat, ilon][1], zvals, pvals)

    ## Check position of p_rho_idz_top
    if all(np.isfinite(pval) for pval in [pvals[idz+1], p_rho_idz_top, pvals[idz-1]]):
        assert pvals[idz+1] < p_rho_idz_top < pvals[idz-1], (
            "p_rho_idz_top is not between p at adjacent theta levels!")


    ## Determine relative position of p_rho_upper and p_target_top
    ## If p_target_top is less than p_rho_idz_top, the top box is a
    ## fraction of box idz + 1.
    ## If p_target_top is greater than p_rho_upper the top box is a
    ## fraction of box idz.
    if p_target_top < p_rho_idz_top:
        top_box_idz = idz+1

    if p_target_top > p_rho_idz_top:
        top_box_idz = idz

    if p_target_top == p_rho_idz_top:
        top_box_idz = idz


    ## Interpolate the altitude at p_target_top and calculate upper and
    ## lower bounds for the altitude of the new top box
    z_top_box_ub = np.interp(p_target_top, pvals[::-1], zvals[::-1])
    z_top_box_lb = box_hgt_bounds[top_box_idz, ilat, ilon][0]


    ## Check that the interpolated upper altitude bound is sensible
    if (p_target_top < p_rho_idz_top) or p_target_top > p_rho_idz_top:
        assert box_hgt_bounds[top_box_idz, ilat, ilon][0] < \
            z_top_box_ub < box_hgt_bounds[top_box_idz, ilat, ilon][1], (
                'z_top_box_ub is not between the correct model levels!')

    if p_target_top == p_rho_idz_top:
        assert box_hgt_bounds[top_box_idz-1, ilat, ilon][0] < \
            z_top_box_ub < box_hgt_bounds[top_box_idz+1, ilat, ilon][1], (
                'z_top_box_ub is not between the correct model levels!')


    ## Calculate the height of the new top box
    h_new_top_box = z_top_box_ub - z_top_box_lb

    ## Check that the height of the new top box is less than
    ## the height of top_box_idz
    hgt_top_box_idz = (box_hgt_bounds[top_box_idz, ilat, ilon][1]
                       - box_hgt_bounds[top_box_idz, ilat, ilon][0])
    assert h_new_top_box <= hgt_top_box_idz, (
        "h_new_top_box is less than than hgt_top_box_idz")

    ## Interpolate the pressure at the new top box lower bound
    p_top_box_lb = np.interp(z_top_box_lb, zvals, pvals)

    ## Estimate pressure new top box centre
    p_centre_new = p_top_box_lb - ((p_top_box_lb - p_target_top) / 2.0)

    ## Check the pressure at the new top box bottom and the new top box centre
    assert p_top_box_lb > p_target_top, (
        "p_at_newtop_bottom is not greater than p_target_top")

    assert p_top_box_lb > p_centre_new > p_target_top, (
        "p_centre_new is not between p_top_box_lb and p_target_top")


    ## Interpolate temperature at the new top box centre
    t_centre_new = np.interp(p_centre_new, pvals[::-1], tvals[::-1])


    ## Set pressure, temperature and species mixing ratio to NaN
    ## above the new top box and insert new pressure, temperature
    ## and box height for the new top box. Species concentration
    ## remains the same in the new top box.
    pvals[top_box_idz:] = np.nan
    pvals[top_box_idz] = p_centre_new

    tvals[top_box_idz:] = np.nan
    tvals[top_box_idz] = t_centre_new

    spvals[top_box_idz:] = np.nan

    hvals[top_box_idz:] = np.nan
    hvals[top_box_idz] = h_new_top_box



if __name__ == '__main__':

    import doctest
    doctest.testmod()
