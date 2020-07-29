#!/usr/bin/env python
"""
Script to test doctests still work in all modules which are known to contain
doctests. Run using ./doctests.py at the command line.
Please ensure you are first using the current release of scitools using:
module load scitools

You can also run these tests on SPICE for faster turn around and lower impact
on your machine. NB to do this the code must be checked out to a networked
file location i.e. NOT /data/local. SPICE CAN'T SEE YOUR /data/local!

sbatch doctests_spice_rhel7.sh


"""
from __future__ import print_function

import doctest
import glob
import warnings
import sys

import six

import adaqcode
import adaqscripts

#Turn off warnings as these clog up output and are not tested in doctests
warnings.filterwarnings("ignore")

#: List of modules known to contain doctests.
#: If a new module containing doctests is added, then it needs to also
#: be included in this list.
DOCTEST_MODULES = [adaqcode.adaq_cmap, adaqcode.adaq_data,
                   adaqcode.adaq_functions, adaqcode.adaq_plotting,
                   adaqcode.adaq_vertical_plotting, adaqcode.aurn_data,
                   adaqcode.aq_indices, adaqcode.array_statistics,
                   adaqcode.camsaqobs_data, adaqcode.chemistry,
                   adaqcode.cube_chemistry, adaqcode.cube_functions,
                   adaqcode.cube_statistics, adaqcode.cube_time,
                   adaqcode.field_layer, adaqcode.field_plot,
                   adaqcode.inifile, adaqcode.line_plot,
                   adaqcode.maccens_data, adaqcode.meteorology,
                   adaqcode.name_data,
                   adaqcode.name_functions,
                   adaqcode.nimrod_data, adaqcode.pollenobs_data,
                   adaqcode.pp_data, adaqcode.plotting_functions,
                   adaqcode.quantile_matching, adaqcode.sites_info,
                   adaqcode.shell_commands,
                   adaqcode.statistics_plotting, adaqcode.timeseries_plot,
                   adaqcode.timeseries_stats, adaqcode.trajectory_data,
                   adaqcode.trajectory_plot,
                   adaqcode.weather_regime, adaqcode.website,
                   adaqscripts.aq_mass_retrieve, adaqscripts.name_field_plot,
                   adaqscripts.aq_plot, adaqscripts.aq_plot_3d,
                   adaqscripts.aq_plot_regimes, adaqscripts.mmr2dobsonunits,
                   adaqscripts.retrieve_aurn, adaqscripts.rolling_stats,
                   adaqscripts.plot_trajectory, adaqscripts.rsmc_plot,
                   adaqscripts.name_postage_stamp_plot,
                   adaqscripts.name_field_plot_diff,
                   adaqscripts.ens_exceedance_plot,
                   adaqscripts.ens_percentile_plot
                  ]
DOCTEST_FILES = glob.glob("adaqdocs/user_guide/*.rst")

# Test ecgrib_data only if iris has the capability to load GRIB files
try:
    import iris_grib
    DOCTEST_MODULES.append(adaqcode.ecgrib_data)
except ImportError:
    GRIB = False


def run_all_doctests():
    """
    Function to run all doctests

    :returns: total_failures, integer, the number of failed doctests.
    """

    total_failures = 0
    total_tests = 0

    for module in DOCTEST_MODULES:
        print('Testing ', module)
        failure_count, test_count = doctest.testmod(module)
        total_failures += failure_count
        total_tests += test_count
        #Ensure output is returned to output file if running under SPICE
        sys.stdout.flush()
        sys.stderr.flush()

    for filename in DOCTEST_FILES:
        print('Testing ', filename)
        failure_count, test_count = doctest.testfile(filename)
        total_failures += failure_count
        total_tests += test_count
        #Ensure output is returned to output file if running under SPICE
        sys.stdout.flush()
        sys.stderr.flush()

    msg = 'total failures: {0} from {1} tests'.format(
        total_failures, total_tests)
    print(msg)

    return total_failures

if __name__ == '__main__':
    if six.PY2:
        print('Using python2')
    else:
        print('Using python3')

    nfailures = run_all_doctests()
