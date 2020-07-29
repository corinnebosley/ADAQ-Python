#!/bin/bash -l
#SBATCH --mem=50000
#SBATCH --ntasks=1
#SBATCH --output=doctests-rhel7-%j.out 
#SBATCH --error=doctests-rhel7-%j.out 
#SBATCH --time=50
#SBATCH --partition=rhel7
#SBATCH --export=NONE

# Runs the doctests on SPICE instead of desktop
# to use this script submit with: 
# sbatch doctests_spice_rhel7.sh
# NB to do this the code must be checked out to a networked
# file location i.e. NOT /data/local. SPICE CAN'T SEE YOUR 
# /data/local!

# Test python 2
echo 'Testing doctests using scitools/default_legacy-current (python2)'
module unload scitools
module load scitools/default_legacy-current
echo $SSS_TAG_DIR

# run doctests
time ./doctests.py

# Test python 3
echo '*********************************************************'
echo 'Testing doctests using scitools/default-current (python3)'
module unload scitools
module load scitools/default-current
echo $SSS_TAG_DIR

# run doctests
time ./doctests.py

echo 'Finished doctests_spice_rhel7.sh'
