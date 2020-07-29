"""
Contains site-specific configuration settings.
"""
import os

#: Data directory containing sample data (site-specific)
SAMPLE_DATADIR = '/data/users/apdg/python_sample_data/'

#: Location of the top level of this version of the code,
#: based on location of this config.py code.
CODE_DIR = os.path.dirname(os.path.realpath(__file__))+'/../'

#: Location of directory that should be writable to by current user
TEST_DIR = os.path.expandvars("$DATADIR")+'/ADAQ_PythonCode/test'
