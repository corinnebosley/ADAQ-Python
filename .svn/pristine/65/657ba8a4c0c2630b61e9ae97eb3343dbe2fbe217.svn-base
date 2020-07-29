
name_data.py
============

.. toctree::
      :maxdepth: 2

.. note::
   The aim of this data class is to read all gridded forms of NAME data with combinations of X, Y, Z and T dimensions. However, due to limitations in iris currently only horizontally gridded and time series data can be loaded into an adaq_data class. Trajectory data will be handled in a separate data class.

NameData Class
--------------     

.. automodule:: name_data
       :members:

.. _field_attributes_ref:

Field Attributes
----------------
In iris all the information found in the column headers in a NAME output file is stored on the iris cube as an attribute. The exception to this is any data which relates to a dimension (currently X, Y, Z and T). To extract a cube with particular attributes it is necessary to supply a dictionary of attributes (in the form of key, value pairs) to the NameData class. For example if we wish to load only air concentrations our dictionary key, value pair would look like this:

>>> field_attributes = {'Quantity':'Air Concentration'}

Possible keys in NAME II gridded output files are:
   * Quantity
   * Species
   * Species Category
   * Time Av or Int

Possible keys in NAME II time series output files are: 
   * Quantity
   * Species
   * Species Category

Possible keys in NAME III output files are:
   * Ensemble Av
   * Horizontal Av or Int
   * Name
   * Quantity
   * Sources
   * Species
   * Species Category
   * Time Av or Int
   * Vertical Av or Int
   * Prob Perc
   * Prob Perc Ens
   * Prob Perc Time
   * D

Possible keys in NAME III output format 2 output files are
   * Species category
   * Field Name
   * Quantity
   * Species
   * Source/source group
   * Particle size range
   * Ensemble av info
   * Time av/int info
   * Horizontal av/int info
   * Vertical av/int info
   * Prob & %-ile info

.. note ::
   Not all cubes will contain all the attributes listed above. For example if no probability output was requested "D" will not be a cube attribute. i.e. The empty lines in the column header in the NAME file will not appear in the cube attributes list.

If you are unsure about the values associated with these keys load your data into an iris cube and print the result. For example:

>>> import iris
>>> cubes = iris.load('Fields_grid4_C1_T1_201208011300.txt')
>>> print cubes[0]

Yields the following output::

  CAESIUM-137_AIR_CONCENTRATION / (Bq / m^3) (latitude: 200; longitude: 200)
       Dimension coordinates:
            latitude                                  x               -
            longitude                                 -               x
       Scalar coordinates:
            height: 50.0 m
            time: 2012-08-01 13:00:00, bound=(2012-08-01 10:00:00, 2012-08-01 13:00:00)
       Attributes:
            End of release: 01/08/2012 13:00 UTC
            Ensemble Av: No ensemble averaging
            Horizontal Av or Int: No horizontal averaging
            Met data: Single Site Flow.SiteFlow
            NAME Version: NAME III (version 5.4.4)
            Name: Unnamed Field Req 8
            Number of field cols: 5
            Number of preliminary cols: 4
            Quantity: Air Concentration
            Release height: 0.000 to 1.000m agl
            Release location: 0.9655E   50.9127N
            Run duration: 3hr 0min
            Run name: DUNGENESS_01082012_0952Z
            Run time: 01/08/2012 11:19:20.802 UTC+01:00
            Source strength: 1.000000 Bq / s
            Sources: All sources
            Species: CAESIUM-137
            Species Category: RADIONUCLIDE
            Start of release: 01/08/2012 10:00 UTC
            Time Av or Int: 3hr 0min average
            Vertical Av or Int: No vertical averaging
       Cell methods:
            mean: time

And we can see that the value associated with the "Species" key is "CAESIUM-137" and the value associated with the "Name" key is "Unnamed Field Req 8".

.. note::
    The keys and values in the dictionary are case sensitive
