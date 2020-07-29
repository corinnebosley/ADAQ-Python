Installation Guide
==================
ADAQ Python requires very little setting up. However, in order to be able to
use ADAQ Python it is necessary to have installed the MOSciTools suite of
python packages. To load these at the Met Office type:

.. code-block:: ksh

    module load scitools

In addition, to be able to run the examples and do
the exercises in the user guide you need to tell ADAQ Python where the 
sample data is stored. To do this open :mod:`config` (see adaqcode/config.py)
and modify the SAMPLE_DATADIR to provide the location of the sample data.

.. note:: On all computers the forward slash '/' should be used to separate
          directories.


.. _cartopy_setup_ref:

Setting Up Mapping in Cartopy
-----------------------------

.. note:: Note that Cartopy is already configured to find the map data 
          on the Met Office Scientific desktop so users of those computers can
          ignore the following instructions.

ADAQ python uses mapping data from Natural Earth (http://www.naturalearthdata.com/).
The python package Cartopy, which is part of the MOSciTools suite of packages, will
download map data for plotting as required. However, if your computer is not
connected to internet permanently or you are likely to plot a large number of maps
the maps can be downloaded from Natural Earth and saved on your computer. There are
two ways to tell Cartopy how to find your downloaded maps.

  1. If you can access the cartopy installation then you can create a siteconfig.py
     file with the cartopy install which sets the config key 'pre_existing_data_dir'
     to the location of the map data.

  2. If you can't access the cartopy installation then it is possible to make use 
     of site.getusersitepackages 
     (see ​https://docs.python.org/2/library/site.html#site.getusersitepackages) 
     and add a file called cartopy_userconfig.py there. This file should contain::

        def update_config(config_dict):
            config_dict[‘pre_existing_data_dir’] = ‘/path/to/my/NaturalEarthFiles’


