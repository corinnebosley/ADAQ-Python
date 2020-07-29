.. _adaqscripts:

Introduction to ADAQ scripts
============================
There are a few scripts available which combine some/all of the ADAQ routines
for particular purposes. These are designed to be run with inifiles.
For example :mod:`name_field_plot`.py does gridded field plotting of NAME data, 
while :mod:`aq_plot`.py does standard verification of Air Quality metrics and
gridded field plots. These scripts are all in the adaqscripts subdirectory 
of ADAQ python. To run these scripts, simply take a copy of the 
inifile with the same name, eg aq_plot.ini (to go with aq_plot.py), 
modify this to point to your required directories and change some options, 
then run, treating it as a unix script, by using:

.. code-block:: ksh

  $  /location/of/adaqcodetrunk/adaqscripts/aq_plot.py /location/of/my/inifile/aq_plot.ini
  
Within the Met Office, there is a copy of the trunk 
kept in ~apdg/PythonCode/ADAQ_Python/trunk. This should be updated whenever 
any code changes are made in the trunk. However this could mean that the
code is changed while you are working with it. It is therefore recommended
that you firstly checkout a copy of the trunk for your own use (just be 
careful not to make any changes which you then commit directly back to the 
trunk!).
To check out a version of the trunk:

.. code-block:: ksh
   
   cd $DATADIR #Or another suitable directory
   fcm checkout fcm:adaq_python_tr
  


.. autosummary::
   :toctree:
   
   aq_mass_retrieve
   aq_plot
   aq_plot_3d
   aq_plot_regimes
   ens_exceedance_plot
   ens_percentile_plot
   mmr2dobsonunits
   name_field_plot
   name_field_plot_diff
   name_postage_stamp_plot
   retrieve_aurn
   rolling_stats
   plot_trajectory
   rsmc_plot
