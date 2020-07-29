"""
Constants for this library
"""
#: Standard temperature. Units K (equivalent to 20C).
#: Recommended by Defra for conversion from ppb to ug/m3 :
#: http://uk-air.defra.gov.uk/air-pollution/faq?question=16
STD_T = 20. + 273.
#: Standard pressure 1013hPa in units of Pa (same source
#: as standard T).
STD_P = 1013.*100.

#: Molar Mass of Air. Units g/mol. Value from:
#: http://www.engineeringtoolbox.com/molecular-mass-air-d_679.html
#: Note this is specific to Earth in the present day. If we want
#: to look at paleo studies or exoplanets this might not be correct.
MMASS_AIR = 28.97

#: Ratio of molar mass of water vapour to dry air :
#: https://en.wikipedia.org/wiki/Density_of_air
MMASS_RATIO = 0.62201

#: Specific Gas Constant. Units J/K/kg :
#: https://en.wikipedia.org/wiki/Density_of_air
SPC_GAS_CONST = 287.058

#: Universal Gas Constant. Units J/K/mol :
#: https://en.wikipedia.org/wiki/Gas_constant
UNIV_GAS_CONST = 8.3144598
