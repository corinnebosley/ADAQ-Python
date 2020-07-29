"""
ADAQ python code
"""

__version__ = 1.0

#Firstly, insure all code will run under cron etc
import os
import sys
import six
import matplotlib as mpl
#if isinstance(sys.stdout, file):
if not sys.stdout.isatty(): #isttay checks if output is to a terminal
    try:
        # Not in a CASE tool
        if os.getpgrp() != os.tcgetpgrp(sys.stdout.fileno()):
            #Called in background
            mpl.use('agg')
    except OSError:
        #Called via nohup
        mpl.use('agg')

#Ensure that modules are available on the path.
PACKAGE_DIR = os.path.dirname(os.path.realpath(__file__))
sys.path.append(PACKAGE_DIR)

#Make key variables also available at top level
import config
from config import SAMPLE_DATADIR, CODE_DIR, TEST_DIR

#Import all submodules
import adaq_cmap
import adaq_data
import adaq_functions
import adaq_plotting
import adaq_vertical_plotting
import array_statistics
import aurn_data
import aq_indices
import camsaqobs_data
import chemistry
import constants
import cube_chemistry
import cube_functions
import cube_statistics
import cube_time
import ecgrib_data
import field_layer
import field_plot
import inifile
import line_plot
import maccens_data
import meteorology
import name_functions
import name_data
import nimrod_data
import plotting_functions
import pollenobs_data
import pp_data
import quantile_matching
import shell_commands
import sites_info
import statistics_plotting
import timeseries_plot
import timeseries_stats
import trajectory_data
import trajectory_plot
import weather_regime
import website
