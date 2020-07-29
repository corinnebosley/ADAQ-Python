# Crown copyright 2014
"""
Functions relating to chemical names and symbols
"""
from __future__ import division
from __future__ import absolute_import
from __future__ import print_function

import re

# ATOMIC_WEIGHTS contains standard atomic weights
# based on the Commission on Isotopic
# Abundances and Atomic Weights's definitive table of the standard
# atomic weights, published in 2014:
# http://www.ciaaw.org/atomic-weights.htm
# Note that unstable isoptopes do not have a standard atomic weight -
# represent this as None and raise an exception for this.
# Alternative spellings of aluminium, caesium and sulfur are included
ATOMIC_WEIGHTS = {
    'actinium'      : None,
    'aluminium'     :  26.9815385,
    'aluminum'      :  26.9815385,
    'americium'     : None,
    'antimony'      : 121.76,
    'argon'         :  39.948,
    'arsenic'       :  74.921595,
    'astatine'      : None,
    'barium'        : 137.327,
    'berkelium'     : None,
    'beryllium'     :   9.0121831,
    'bismuth'       : 208.9804,
    'bohrium'       : None,
    'boron'         :  10.8135,
    'bromine'       :  79.904,
    'cadmium'       : 112.414,
    'caesium'       : 132.905452,
    'cesium'        : 132.905452,
    'calcium'       :  40.078,
    'californium'   : None,
    'carbon'        :  12.0106,
    'cerium'        : 140.116,
    'chlorine'      :  35.4515,
    'chromium'      :  51.9961,
    'cobalt'        :  58.933194,
    'copernicium'   : None,
    'copper'        :  63.546,
    'curium'        : None,
    'darmstadtium'  : None,
    'dubnium'       : None,
    'dysprosium'    : 162.5,
    'einsteinium'   : None,
    'erbium'        : 167.259,
    'europium'      : 151.964,
    'fermium '      : None,
    'flerovium'     : None,
    'fluorine'      :   8.99840316,
    'francium'      : None,
    'gadolinium'    : 157.25,
    'gallium'       :  69.723,
    'germanium'     :  72.63,
    'gold'          : 196.966569,
    'hafnium'       : 178.49,
    'hassium'       : None,
    'helium'        :   4.002602,
    'holmium'       : 164.93033,
    'hydrogen'      :   1.00811,
    'indium'        : 114.818,
    'iodine'        : 126.90447,
    'iridium'       : 192.217,
    'iron'          :  55.845,
    'krypton'       :  83.798,
    'lanthanum'     : 138.90547,
    'lawrencium'    : None,
    'lead'          : 207.2,
    'lithium'       :   6.9675,
    'livermorium'   : None,
    'lutetium'      : 174.9668,
    'magnesium'     :  24.3055,
    'manganese'     :  54.938044,
    'meitnerium'    : None,
    'mendelevium'   : None,
    'mercury'       : 200.592,
    'molybdenum'    :  95.95,
    'neodymium'     : 144.242,
    'neon'          :  20.1797,
    'neptunium'     : None,
    'nickel'        :  58.6934,
    'niobium'       :  92.90637,
    'nitrogen'      :  14.00728,
    'nobelium'      : None,
    'osmium'        : 190.23,
    'oxygen'        :  15.99977,
    'palladium'     : 106.42,
    'phosphorus'    :  30.973762,
    'platinum'      : 195.084,
    'plutonium'     : None,
    'polonium'      : None,
    'potassium'     :  39.0983,
    'praseodymium'  : 140.90766,
    'promethium'    : None,
    'protactinium'  : 231.03588,
    'radium'        : None,
    'radon'         : None,
    'rhenium'       : 186.207,
    'rhodium'       : 102.9055,
    'roentgenium'   : None,
    'rubidium'      :  85.4678,
    'ruthenium'     : 101.07,
    'rutherfordium' : None,
    'samarium'      : 150.36,
    'scandium'      :  44.955908,
    'seaborgium'    : None,
    'selenium'      :  78.971,
    'silicon'       :  28.085,
    'silver'        : 107.8682,
    'sodium'        :   2.98976928,
    'strontium'     :  87.62,
    'sulfur'        :  32.0675,
    'sulphur'       :  32.0675,
    'tantalum'      : 180.94788,
    'technetium'    : None,
    'tellurium'     : 127.6,
    'terbium'       : 158.92535,
    'thallium'      : 204.3835,
    'thorium'       : 232.0377,
    'thulium'       : 168.93422,
    'tin'           : 118.71,
    'titanium'      :  47.867,
    'tungsten'      : 183.84,
    'ununoctium'    : None,
    'ununpentium'   : None,
    'ununseptium'   : None,
    'ununtrium'     : None,
    'uranium'       : 238.02891,
    'vanadium'      :  50.9415,
    'xenon'         : 131.293,
    'ytterbium'     : 173.054,
    'yttrium'       :  88.90584,
    'zinc'          :  65.38,
    'zirconium'     :  91.224}

# ATOMIC_NAME is a dictionary which translates from
# symbols of elements to names.
ATOMIC_NAMES = {
    'Ac' : 'actinium',
    'Ag' : 'silver',
    'Al' : 'aluminium',
    'Am' : 'americium',
    'Ar' : 'argon',
    'As' : 'arsenic',
    'At' : 'astatine',
    'Au' : 'gold',
    'B'  : 'boron',
    'Ba' : 'barium',
    'Be' : 'beryllium',
    'Bh' : 'bohrium',
    'Bi' : 'bismuth',
    'Bk' : 'berkelium',
    'Br' : 'bromine',
    'C'  : 'carbon',
    'Ca' : 'calcium',
    'Cd' : 'cadmium',
    'Ce' : 'cerium',
    'Cf' : 'californium',
    'Cl' : 'chlorine',
    'Cm' : 'curium',
    'Co' : 'cobalt',
    'Cr' : 'chromium',
    'Cs' : 'caesium',
    'Cu' : 'copper',
    'Db' : 'dubnium',
    'Ds' : 'darmstadtium',
    'Dy' : 'dysprosium',
    'Er' : 'erbium',
    'Es' : 'einsteinium',
    'Eu' : 'europium',
    'F'  : 'fluorine',
    'Fe' : 'iron',
    'Fm' : 'fermium',
    'Fr' : 'francium',
    'Ga' : 'gallium',
    'Gd' : 'gadolinium',
    'Ge' : 'germanium',
    'H'  : 'hydrogen',
    'He' : 'helium',
    'Hf' : 'hafnium',
    'Hg' : 'mercury',
    'Ho' : 'holmium',
    'Hs' : 'hassium',
    'I'  : 'iodine',
    'In' : 'indium',
    'Ir' : 'iridium',
    'K'  : 'potassium',
    'Kr' : 'krypton',
    'La' : 'lanthanum',
    'Li' : 'lithium',
    'Lr' : 'lawrencium',
    'Lu' : 'lutetium',
    'Md' : 'mendelevium',
    'Mg' : 'magnesium',
    'Mn' : 'manganese',
    'Mo' : 'molybdenum',
    'Mt' : 'meitnerium',
    'N'  : 'nitrogen',
    'Na' : 'sodium',
    'Nb' : 'niobium',
    'Nd' : 'neodymium',
    'Ne' : 'neon',
    'Ni' : 'nickel',
    'No' : 'nobelium',
    'Np' : 'neptunium',
    'O'  : 'oxygen',
    'Os' : 'osmium',
    'P'  : 'phosphorus',
    'Pa' : 'protactinium',
    'Pb' : 'lead',
    'Pd' : 'palladium',
    'Pm' : 'promethium',
    'Po' : 'polonium',
    'Pr' : 'praseodymium',
    'Pt' : 'platinum',
    'Pu' : 'plutonium',
    'Ra' : 'radium',
    'Rb' : 'rubidium',
    'Re' : 'rhenium',
    'Rf' : 'rutherfordium',
    'Rg' : 'roentgenium',
    'Rh' : 'rhodium',
    'Rn' : 'radon',
    'Ru' : 'ruthenium',
    'S'  : 'sulfur',
    'Sb' : 'antimony',
    'Sc' : 'scandium',
    'Se' : 'selenium',
    'Sg' : 'seaborgium',
    'Si' : 'silicon',
    'Sm' : 'samarium',
    'Sn' : 'tin',
    'Sr' : 'strontium',
    'Ta' : 'tantalum',
    'Tb' : 'terbium',
    'Tc' : 'technetium',
    'Te' : 'tellurium',
    'Th' : 'thorium',
    'Ti' : 'titanium',
    'Tl' : 'thallium',
    'Tm' : 'thulium',
    'U'  : 'uranium',
    'Uub' : 'ununbium',
    'Uuh' : 'ununhexium',
    'Uuo' : 'ununoctium',
    'Uup' : 'ununpentium',
    'Uuq' : 'ununquadium',
    'Uus' : 'ununseptium',
    'Uut' : 'ununtrium',
    'V'  : 'vanadium',
    'W'  : 'tungsten',
    'Xe' : 'xenon',
    'Y'  : 'yttrium',
    'Yb' : 'ytterbium',
    'Zn' : 'zinc',
    'Zr' : 'zirconium'
}

# ATOMIC_SYM is a dictionary which translates from
# names of elements to their symbols. We
# create the dictionary by inverting ATOMIC_NAMES
ATOMIC_SYM = {v: k for k, v in ATOMIC_NAMES.items()}

# CHEM_FORM is a dictionary of molecular formulae
# to be used to calculate molar masses and also
# when converting masses expressed as mass of an element
CHEM_FORM = {
    'CFC11' : 'CFCl3',
    'CFC113' : 'CCl2FCClF2',
    'CFC113a' : 'CCl3CF3',
    'CFC114' : 'CClF2CClF2',
    'CFC115' : 'CClF2CF3',
    'CFC12' : 'CF2Cl2',
    'HCC140a' : 'CH3CCl3',
    'HCFC141b' : 'CH3CCl2F',
    'HCFC142b' : 'CH3CClF2',
    'HCFC22' : 'CHClF2',
    'acetic acid' : 'CH3COOH',
    'aceto-nitrile' : 'CH3CN',
    'alpha hexachlorocyclohexane' : 'C6H6Cl6',
    'alpha pinene' : 'C10H16',
    'ammonia' : 'NH3',
    'ammonium' : 'NH4',
    'ammonium nitrate' : 'NH4NO3',
    'ammonium sulfate' : 'NH4NH4SO4',
    'aragonite' : 'CaCO3',
    'benzene' : 'C6H6',
    'beta pinene' : 'C10H16',
    'bromine chloride' : 'BrCl',
    'bromine monoxide' : 'BrO',
    'bromine nitrate' : 'BrONO2',
    'butane' : 'C4H10',
    'calcite' : 'CaCO3',
    'carbon dioxide' : 'CO2',
    'carbon monoxide' : 'CO',
    'carbon tetrachloride' : 'CCl4',
    'carbonate anion' : 'CO3',
    'chlorine dioxide' : 'OClO',
    'chlorine monoxide' : 'ClO',
    'chlorine nitrate' : 'ClONO2',
    'chlorophyll-a' : 'C55H72O5N4Mg',
    'dichlorine peroxide' : 'Cl2O2',
    'dichlorineperoxide' : 'Cl2O2',
    'dimethyl sulfide' : 'CH32SCH3',
    'dinitrogen pentoxide' : 'N2O5',
    'dinitrogenpentoxide' : 'N2O5',
    'ethane' : 'C2H6',
    'ethanol' : 'C2H5OH',
    'ethene' : 'C2H4',
    'ethyne' : 'HC2H',
    'formaldehyde' : 'CH2O',
    'formic acid' : 'HCOOH',
    'glyoxal' : 'CHOCHO',
    'halo2402' : 'C2Br2F4',
    'halon1202' : 'CBr2F2',
    'halon1211' : 'CBrClF2',
    'halon1301' : 'CBrF3',
    'halon2402' : 'C2Br2F4',
    'hcc140a' : 'CH3CCl3',
    'hexachlorobiphenyl' : 'C12H4Cl6',
    'hydrogen bromide' : 'HBr',
    'hydrogen chloride' : 'HCl',
    'hydrogen cyanide' : 'HCN',
    'hydrogen peroxide' : 'H202',
    'hydrogen sulfide' : 'H2S',
    'hydroperoxyl radical' : 'HO2',
    'hydroxyl radical' : 'OH',
    'hypobromous acid' : 'HOBr',
    'hypochlorous acid' : 'HOCl',
    'iodine monoxide' : 'IO',
    'isoprene' : 'C5H8',
    'limonene' : 'C10H16',
    'methane' : 'CH4',
    'methanol' : 'CH3OH',
    'methyl bromide' : 'CH3Br',
    'methyl chloride' : 'CH3Cl',
    'methyl hydroperoxide' : 'CH3OOH',
    'methyl peroxy radical' : 'CH3O2',
    'molecular hydrogen' : 'H2',
    'nitrate anion' : 'NO3',
    'nitrate' : 'NO3',
    'nitric acid' : 'HNO3',
    'nitrite anion' : 'NO2',
    'nitrogen dioxide' : 'NO2',
    'nitrogen monoxide' : 'NO',
    'nitrous acid' : 'HNO2',
    'nitrous oxide' : 'N2O',
    'ozone' : 'O3',
    'pentane' : 'C5H12',
    'peroxyacetyl nitrate' : 'CH3COO2NO2',
    'peroxynitric acid' : 'HNO4',
    'phosphate anion' : 'PO4',
    'propane' : 'C3H8',
    'propene' : 'C3H6',
    'sulfate anion' : 'SO4',
    'sulfur dioxide' : 'SO2',
    'toluene' : 'C6H5CH3',
    'trimethylbenzene' : 'C9H12',
    'water vapor' : 'H2O',
    'xylene' : 'C6H4C2H6',
    'acetaldehyde' : 'CH3CHO',
    'acetone' : 'CH3COCH3',
    'acetonyl peroxy radical' : 'CH3COCH2O2',
    'acetonylhydroperoxide' : 'CH3COCH2OOH',
    'atomic bromine' : 'Br',
    'atomic chlorine' : 'Cl',
    'atomic hydrogen' : 'H',
    'atomic nitrogen' : 'H',
    'carbon disulfide' : 'CS2',
    'carbonyl sulfide' : 'OCS',
    'dimethyl sulfoxide' : 'CH3SOCH3',
    'ethyl hydroperoxide' : ' C2H6O2',
    'ethyl peroxy radical' : 'C2H6O2',
    'ground state atomic oxygen' : 'O',
    'hydroxyacetone' : 'C3H6O2',
    'i-propyl hydroperoxide' : 'C3H8O2',
    'methacrolein' : 'C4H6O',
    'methanesulfonic acid' : 'CH3SO3H',
    'methlyglyoxal' : 'C3H4O2',
    'methyl ethyl ketone' : 'C4H8O',
    'methyl nitrate' : 'CH3NO3',
    'nitrate radical' : 'NO3',
    'peracetic acid' : 'CH3CO3H',
    'peroxyacetyl radical' : 'CH3COO2',
    'peroxypropanoyl radical' : 'C3H5O2',
    'peroxypropionyl nitrate' : 'C3H5NO5',
    'propanal' : 'CH3CH2CHO',
    'sulfuric acid' : 'H2SO4'
    }

#Dictionary of species (given as a shortname) that may be expressed
# as a component (full name that is a CHEM_SPECIES key)
EXPRESSED_AS = {
    "nox": "nitrogen dioxide",
    "ox": "ground state atomic oxygen",
    "nmvoc": "carbon",
    }

def get_atomic_symbol(name):
    """
    Given the symbol of an element, return its name.

    >>> get_atomic_symbol('Sulfur')
    'S'

    Code is case insenstitive

    >>> get_atomic_symbol('oXygeN')
    'O'

    Code now supports alternate element spellings
    >>> get_atomic_symbol('Sulphur') # doctest: +ELLIPSIS
    'S'

    >>> get_atomic_symbol('cesium') # doctest: +ELLIPSIS
    'Cs'

    >>> get_atomic_symbol('unobtainium') # doctest: +ELLIPSIS
    Traceback (most recent call last):
      File ".../doctest.py", line 1315, in __run
        compileflags, 1) in test.globs
      File "<doctest __main__.get_atomic_symbol[4]>", line 1, in <module>
        get_atomic_symbol('unobtainium') # doctest: +ELLIPSIS
      File ".../chemistry.py", line 427, in get_atomic_symbol
        raise ValueError, "Unknown element name: " + name
    ValueError: Unknown element name: unobtainium

    """
    # Convert to lower case
    name = name.lower()
    if name == 'sulphur':
        name = 'sulfur'
    elif name == 'aluminum':
        name = 'aluminium'
    elif name == 'cesium':
        name = 'caesium'
    try:
        symbol = ATOMIC_SYM[name]
    except:
        raise ValueError("Unknown element name: " + name)

    return symbol

def get_atomic_name(symbol):
    """
    Given an atomic symbol, return IUPAC name of element

    >>> get_atomic_name('S')
    'sulfur'

    >>> get_atomic_name('W')
    'tungsten'

    >>> get_atomic_name('Zn')
    'zinc'

    """
    try:
        name = ATOMIC_NAMES[symbol]
    except:
        raise ValueError("Unknown element symbol: " + symbol)

    return name

def get_atomic_weight(name):
    """
    Return the atomic weight of an element given its name

    >>> get_atomic_weight('aluminium')
    26.9815385

    The names are case insensitive

    >>> get_atomic_weight('Carbon')
    12.0106

    >>> get_atomic_weight('HELIUM')
    4.002602

    Alternate spellings are also supported

    >>> get_atomic_weight('aluminum')
    26.9815385

    >>> get_atomic_weight('sulfur')
    32.0675

    >>> get_atomic_weight('sulphur')
    32.0675

    >>> get_atomic_weight('caesium')
    132.905452

    >>> get_atomic_weight('cesium')
    132.905452

    If the element is unknown, returns a ValueError

    >>> get_atomic_weight('bogon') # doctest: +ELLIPSIS
    Traceback (most recent call last):
      File ".../doctest.py", line 1289, in __run
        compileflags, 1) in test.globs
      File "<doctest __main__.get_atomic_weight[3]>", line 1, in <module>
        get_atomic_weight('bogon')
      File "standard_atomic_weights.py", line 166, in get_atomic_weight
        raise ValueError, "Unknown element: " + name
    ValueError: Unknown element: bogon

    If the element is unstable there is no standard atomic
    weight defined so also returns an exception

    >>> get_atomic_weight('actinium') # doctest: +ELLIPSIS
    Traceback (most recent call last):
      File ".../doctest.py", line 1289, in __run
        compileflags, 1) in test.globs
      File "<doctest __main__.get_atomic_weight[4]>", line 1, in <module>
        get_atomic_weight('actinium')
      File "standard_atomic_weights.py", line 182, in get_atomic_weight
        raise ValueError, "Unstable element: " + name
    ValueError: Unstable element: actinium

    """

    # Now use the values from the dictionary to
    # return the required value
    try:
        weight = ATOMIC_WEIGHTS[name.lower()]

    # If there is a key error in the dictionary,
    # then the element is unknown. Raise a ValueError
    except KeyError:
        raise ValueError("Unknown element: " + name)
    # If the weight returned is None then
    # it is unstable and the weight is not defined
    if weight is None:
        raise ValueError("Unstable element: " + name)
    return weight


def get_atomic_weight_from_symbol(symbol):
    """
    Return the atomic weight of an element given
    its chemical symbol

    >>> get_atomic_weight_from_symbol('C')
    12.0106

    >>> get_atomic_weight_from_symbol('He')
    4.002602

    If an unknown element is requested, gives a value error

    >>> get_atomic_weight_from_symbol('X')  # doctest: +ELLIPSIS
    Traceback (most recent call last):
      File ".../doctest.py", line 1289, in __run
        compileflags, 1) in test.globs
      File "<doctest __main__.get_atomic_weight_from_symbol[1]>", line 1, \
in <module>
        get_atomic_weight_from_symbol('X')
      File ".../r111_cube_transforms/adaqcode/standard_atomic_weights.py", \
line 348, in get_atomic_weight_from_symbol
        raise ValueError, "Unknown element symbol: " + symbol
    ValueError: Unknown element symbol: X

    As with requesting by name, unstable elements give an error

    >>> get_atomic_weight_from_symbol('Uut')  # doctest: +ELLIPSIS
    Traceback (most recent call last):
      File ".../doctest.py", line 1289, in __run
        compileflags, 1) in test.globs
      File "<doctest __main__.get_atomic_weight_from_symbol[2]>", line 1, \
in <module>
        get_atomic_weight_from_symbol('Uut')
      File ".../standard_atomic_weights.py", line 361, in get_atomic_weight_\
from_symbol
        return get_atomic_weight(name)
      File ".../standard_atomic_weights.py", line 329, in get_atomic_weight
        raise ValueError, "Unstable element: " + name
    ValueError: Unstable element: ununtrium

    """
    return get_atomic_weight(get_atomic_name(symbol))

def count_elements(formula):
    """
    Take a chemical formula as input
    and return a dictionary with a count of the
    names and number of each element in it

    >>> print(count_elements('O3'))
    {'O': 3}

    >>> count_elements('C55H72O5N4Mg') == {'C': 55, 'H': 72,
    ... 'O': 5, 'N': 4, 'Mg': 1}
    True

    If the string cannot be made sense of as
    a chemical formula an error will be raised

    >>> print(count_elements('H3PoooO5'))  # doctest: +ELLIPSIS
    Traceback (most recent call last):
       File ".../doctest.py", line 1289, in __run
         compileflags, 1) in test.globs
       File "<doctest __main__.count_elements[2]>", line 1, in <module>
         print count_elements('H3PoooO5')
       File ".../chemistry.py", line 457, in count_elements
         raise ValueError, "malformed formula: " + formula
    ValueError: malformed formula: H3PoooO5

    """

    # using a regular expression, split up
    # the name into symbols with an Element
    # name followed by a number e.g. O2, Fe, Ca3
    # in alternate slots. If things don't match
    # the pattern, they will go to the
    # rubbish array

    pattern = r'([A-Z][a-z]{0,1}\d*)'

    pieces = re.split(pattern, formula)
    data = pieces[1::2]
    rubbish = [item for item in pieces[0::2] if item]

    # if we get anything in the rubbish list, then this
    # is an error
    if rubbish:
        raise ValueError("malformed formula: " + formula)

    # now split the valid data into an integer and
    # a chemical symbol and put a count of the
    # elements into a dictionary

    # create empty dictionary to hold result
    counted = {}

    # This pattern captures the element and
    # the number seperately

    little_pattern = r'([A-Z][a-z]{0,1})(\d*)'

    # Now loop over these element - number strings,
    # and split into an element name and an int
    # If there is no number this mean the number is 1
    # e.g. H2O means 2 x H and 1 x O

    for atoms in data:
        match = re.search(little_pattern, atoms)

       # this should always match but just in case...
        if match:

            # The element symbol is the first match
            element = match.group(1)

            # the count of that element is the second match
            # if the string is zero length => 1
            if not match.group(2):
                count = 1
            else:
                count = int(match.group(2))

            # now increment dictionary with these values
            if element in counted:
                counted[element] += count
            else:
                counted[element] = count
    return counted

def get_molecule_name(name):
    """
    Takes a name and checks if it matches a pattern for something that
    might be a molecule.
    If it does match one of these return the potential molecule name.

    >>> name = 'mass_fraction_of_nitrogen_dioxide_in_air'
    >>> print(name, get_molecule_name(name))
    mass_fraction_of_nitrogen_dioxide_in_air nitrogen_dioxide

    >>> name = 'atmosphere_mole_content_of_ozone'
    >>> print(name, get_molecule_name(name))
    atmosphere_mole_content_of_ozone ozone

    """

    # Set up a list of the patterns which look like
    # they might be chemical compounds
    pattern_pre = [
        r"^atmosphere_moles_of_(.+)",
        r"^atmosphere_mass_content_of_(.+)",
        r"^atmosphere_mole_content_of_(.+)",
        r"^mass_fraction_of_(.+)_in_.+",
        # note that we need the more restrictive
        # match first, so that if there are additional
        # parts at the end these match and we only
        # capture the compound name
        r"^mole_content_of_(.+)_in_.+",
        r"^mole_content_of_(.+)",
        r"^mole_concentration_of_(.+)_in_.+",
        r"^mass_concentration_of_(.+)_in_.+",
        r"^mole_fraction_of_(.+)_in_.+",
        r"^mole_flux_of_(.+)",
        r"^moles_of_(.+)_in_.+",
        r"^moles_of_(.+)"
        ]
    # loop over all the potential patterns
    for pattern in pattern_pre:
        # test the input using each pattern in turn
        match = re.search(pattern, name)
        # if there is a match then we have captured the name we want
        if match:
            molecule = match.group(1)
            # sometime we might have a tendency which ends in us "due_to_"...
            # check for this and strip off the end
            match = re.search(r'(.+)_due_to_.+', molecule)
            if match:
                molecule = match.group(1)
            return molecule
    # if we get here we have failed to match
    # raise a suitable exception...
    raise ValueError('Error in get_molecule_name. Name does not match: '
                     + name)

def get_molmass(molecule):
    """
    function which returns the molecular mass of a chemical
    species based on its name

    >>> molecule='nitrogen_dioxide'
    >>> print(molecule, get_molmass(molecule))
    nitrogen_dioxide 46.00682


    If attempt to get the molecular mass of a cube where the molecule
    name is not recognised, a value error is raised

    To extend this function to include an additional molecule
    add the molecular mass to the molmass dictionary

    >>> molecule='spam'
    >>> print(molecule, get_molmass(molecule)) # doctest: +ELLIPSIS
    Traceback (most recent call last):
      File ".../doctest.py", line 1289, in __run
        compileflags, 1) in test.globs
      File "<doctest __main__.get_molmass[5]>", line 1, in <module>
        print molecule, get_molmass(molecule)
      File ".../chemistry.py", line 704, in get_molmass
        raise ValueError('Error in get_molmass. Unknown molecule: '+molecule)
    ValueError: Error in get_molmass. Unknown molecule: spam

    """

    # note that we need to replace underscores with spaces
    molecule = molecule.replace('_', ' ')

    # Now try and look up the molecule passed in
    # if it is not in the dictionary raise an exception
    # with a helpful message
    try:
        formula = CHEM_FORM[molecule]
    except:
        raise ValueError('Error in get_molmass. Unknown molecule: '+molecule)

    # count the elements
    try:
        count = count_elements(formula)
    except:
        raise ValueError('Error in get_molmass' +
                         ' Cant count elements in  ' + molecule +
                         ' with formula ' + formula)

    mol_mass = 0.
    for key in count:
        try:
            mol_mass += get_atomic_weight_from_symbol(key) * count[key]
        except:
            raise ValueError('Error in get_molmass. '
                             'Unknown element: ' + key)
    return mol_mass

def elemental_to_molecular(name):
    """
    Function to help with conversion from mass fractions expressed as
    mass of a single elemental in a compound to mass of the molecule

    such as from:

    mass_fraction_of_sulfur_dioxide_in_air_expressed_as_sulfur

    to:

    mass_fraction_of_sulfur_dioxide_in_air

    Input:

    name - the name of the phenomenon to be converted

    Returns:

    convfac - the scaling factor to use for conversions

    newname - the name to convert to

    Method:

    First check if the name matches the correct pattern
    If so, calculate the mass of the element, the of that element in
    the molecule and the mass of the molecule. Then scale
    the mass by the ratio mass_molecule/mass_element and rename the
    cube

    >>> convfac, newname = elemental_to_molecular(
    ... 'mass_fraction_of_sulfur_dioxide_expressed_as_sulfur_in_air')
    >>> print("{:10,.8f} {}".format(convfac, newname))
    1.99788072 mass_fraction_of_sulfur_dioxide_in_air

    >>> convfac, newname = elemental_to_molecular(
    ... 'mass_fraction_of_nitrogen_dioxide_expressed_as_nitrogen_in_air')
    >>> print("{:10,.8f} {}".format(convfac, newname))
    3.28449349 mass_fraction_of_nitrogen_dioxide_in_air

    If the name does not match any of the patterns this is an error:

    >>> convfac, newname = elemental_to_molecular(
    ... 'mass_fraction_of_sulfur_dioxide_in_air') # doctest: +ELLIPSIS
    Traceback (most recent call last):
      File ".../doctest.py", line 1289, in __run
        compileflags, 1) in test.globs
      File "<doctest __main__.elemental_to_molecular[2]>", line 1, in <module>
        mol_mass, elemental_mass, newname = elemental_to_molecular(\
'mass_fraction_of_sulfur_dioxide_in_air')
      File ".../chemistry.py", line 863, in elemental_to_molecular
        + 'Unable to convert: ' + name)
    ValueError: Error in elemental_to_molecular.Unable to convert: \
mass_fraction_of_sulfur_dioxide_in_air

    Various other problems with the name cause errors:

    >>> convfac, newname = elemental_to_molecular(
    ... 'mass_fraction_of_sulfur_dioxide_expressed_as_nitrogen_in_air')
    ... # doctest: +ELLIPSIS
    Traceback (most recent call last):
      File ".../doctest.py", line 1289, in __run
        compileflags, 1) in test.globs
      File "<doctest __main__.elemental_to_molecular[5]>", line 1, in <module>
        convfac, newname = elemental_to_molecular(\
'mass_fraction_of_sulfur_dioxide_expressed_as_nitrogen_in_air')
      File ".../adaqcode/chemistry.py", line 886, in elemental_to_molecular
        ' for ' + name)
    ValueError: Error in elemental_to_molecular. Element nitrogen not in \
formula: SO2 for mass_fraction_of_sulfur_dioxide_expressed_as_nitrogen_in_air


    """
    # Set up a list of the patterns which we will try to convert
    # this is deliberately very restrictive

    pattern_list = [
        r"^(atmosphere_mass_content_of_)(.+)_expressed_as_(.+)",
        r"^(mass_concentration_of_)(.+)_expressed_as_(.+)(_in_air)$",
        r"^(mass_fraction_of_)(.+)_expressed_as_(.+)(_in_air)$"
        ]

    # loop over all the potential patterns
    for pattern in pattern_list:

        # test the input using each pattern in turn
        match = re.search(pattern, name)

        if match:
            # construct new name from pattern
            name1 = match.group(1)
            molecule = match.group(2)

            element_name = match.group(3)
            newname = name1 + molecule

            # check if this match is for _in_air
            # if so, add that to the end
            if len(match.groups()) == 4:
                newname = newname + match.group(4)

            # get the mass of the molecule
            try:
                mol_mass = get_molmass(molecule)
            except:
                raise ValueError('Error in elemental_to_molecular.'
                                 + ' Unknown molecule: ' + molecule)

            # Now try and look up the formula for the molecule passed in
            # if it is not in the dictionary raise an exception
            # with a helpful message

            try:
                # note that to use dictionary, replace underscores with spaces
                formula = CHEM_FORM[molecule.replace('_', ' ')]
            except:
                raise ValueError('Error in elemental_to_molecular.'
                                 + ' Unknown molecule: ' + molecule)

            # find out how many atoms of the element
            # are in the molecule

            # count all the elements to get a dictionary
            try:
                count = count_elements(formula)
            except:
                raise ValueError('Error in elemental_to_molecular.'
                                 + ' Problem with formula: ' + formula)

            # we have an element name and we need an element
            # symbol to use with the dictionary
            try:
                element_symb = get_atomic_symbol(element_name)
            except:
                raise ValueError('Error in elemental_to_molecular.'
                                 + ' Problem with getting symbol for element: '
                                 + element_name)

            # Now we query the dictionary to find out the number of
            # atoms of the element in our molecule
            try:
                n_element = count[element_symb]
            except:
                raise ValueError('Error in elemental_to_molecular.'
                                 + ' Element ' +  element_name +
                                 ' not in formula: ' + formula +
                                 ' for ' + name)

            # We can now find the elemental mass in the molecule
            elemental_mass = get_atomic_weight_from_symbol(element_symb) \
                             * n_element

            # the conversion factor is the mass of the the molecule
            # divided by the total mass of the element in that
            # molecule
            convfac = (mol_mass / elemental_mass)

            return convfac, newname

    # If we get here none of the patterns match
    raise ValueError('Error in elemental_to_molecular.'
                     + 'Unable to convert: ' + name)


if __name__ == "__main__":

    import doctest
    doctest.testmod()
