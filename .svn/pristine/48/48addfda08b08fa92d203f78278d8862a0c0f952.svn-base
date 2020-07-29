Further Information
===================

.. _directory_structure:

Directory structure
^^^^^^^^^^^^^^^^^^^

**adaqcode**
  * This should contain all code that may want to be called from elsewhere.
  * Any code in here should be referenced in the __init__.py file
  * Should generally be based around the ADAQData class where possible.
  
**adaqscripts**
  * These are stand-alone scripts, that make use of adaqcode
  * These should generally be called using filename.py inifilename.ini
  
**adaqdocs**
  * This contains extra documentation
  * Note documentation specific to individual modules should live with the module itself.
  * This is split into several subdirectories:
     * user_guide - this contains the user guide, which is written as testable rst files
     * developers_guide - instructions for developers
     * figures - image files for use in the gallery or other documentation.

**adaqsandpit**
  * This contains all code which is useful or likely to be useful to several members of ADAQ
    but has not yet been integrated into adaqcode/adaqscripts
  * Examples of data used by this code can be placed in the sandpit_data directory
  * Documentation should be added at the top of each file

**top-level**
  * Contains pylintrc for checking coding standards
  * Contains top-level documentation-index and code to build documentation.
     

.. _adding_dataclass:

Creating a new class for reading in data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Introduction
------------

All new data formats that are read in, should be read either directly into an :class:`adaq_data.ADAQData` object, 
or more likely, a subclass of ADAQData. As a subclass it will still have all the same methods and
functionality of :class:`adaq_data.ADAQData`, but with extra details specific to this file format, such as how 
read in this format of data and to convert to short_name. If you are new to classes, have a look at this good
`basic introduction <http://www.pythoncentral.io/series/python-classes-tutorial/>`_ about classes and subclasses 
(aka dervied classes).

There are two compulsory methods for a new data format class, which are :ref:`dataclass_init`, to declare
the class as a subclass of :class:`adaq_data.ADAQData`, plus a :ref:`dataclass_readdata` method to contain specifics on 
reading this data format. 
If reading gridded data, the :ref:`dataclass_extract_sites` method is highly recommended.

newformat_data.py
-----------------

Any new data formats should go in a new sensibly named file, for example newformat_data.py (all lowercase).
This module for 'newformat' should start with some basic descriptive documentation and importing of adaq_data:

.. code-block:: python

   """
   Class for reading from newformat data files into an ADAQData class.
   """ 

   import adaq_data

Some dictionaries should then be included to allow conversion from the cube name or other identifier, 
that will allow a CF compliant standard name and a short_name to be set. For example:

.. code-block:: python
   
   #: Conversion between newformat names and CF standard names
   NEWFORMAT_2_STDNAME = {
                          'newformat name for ozone' : 'mass_concentration_of_ozone_in_air',
                          'newformat name for PM10' : 'mass_concentration_of_pm10_ambient_aerosol_in_air'
                         } 
   
   #: Conversion between newformat names and short names
   NEWFORMAT_2_SHORTNAME = {
                            'newformat name for ozone' : 'O3',
                            'newformat name for PM10' : 'PM10'
                           } 

Starting data class
-------------------
                         
The new class should then be opened as a subclass of :class:`adaq_data.ADAQData`, also with some 
documentation, which ideally will eventually be used as an example of how to use this class. The class
name should be given in CamelCase.

.. code-block:: python
   
   class NewFormatData(adaq_data.ADAQData):
       """
       Subclass of ADAQData, which contains extra functionality specific to newformat data.
       
       ** Example: **
       Add example of using this code later, using doctests.
       """
    
.. _dataclass_init:

__init__()
----------

The first method should be an __init__ method which uses the __init__ method from :class:`adaq_data.ADAQData`,
as well as adding any newformat specific attributes. For example:

.. code-block:: python
   
    def __init__(self, label='newformat'):
        """
        Initiates a class from newformat data
        as a subclass of :class:`adaq_data.ADAQData`.  
        """

        adaq_data.ADAQData.__init__(self)
        
        #Now some optional extra attributes specific to this type of data:

        self.label           = label #Label
        self.short_name_list = None #List of short names
        self.start_datetime  = None #Starting datetime
        self.end_datetime    = None #End datetime
        self.filenames       = None #List of raw data filenames
        self.sites_data      = None #sites_data object

.. _dataclass_readdata:

readdata()
----------

There should then be a method for reading data in. This should have keywords which can be used
for setting constraints on the loading of cubes. Note all methods within a class should take
self as its first argument. This method should return a gridded_cube_list, or possibly 
a sites_cube_list depending on the data type. Here we assume we are reading in gridded data.

.. code-block:: python
   
    def readdata(self, filenames=None, short_name_list=None, 
                 start_datetime=None, end_datetime=None):
        """
        Create gridded_cube_list from filenames specified.
        
        :param filenames: list of files to read, can include wildcards etc
        :param short_name_list: list of short names (None gets everything)
        :param start_datetime: Start time of data in datetime format.           
        :param start_datetime: End time of data in datetime format.     
        """      
        
        #Set class attributes on basis of keywords if set
        if filenames is not None:
            self.filenames = filenames
        if short_name_list is not None:    
            self.short_name_list = short_name_list
        if start_datetime is not None:
            self.start_datetime = start_datetime
        if end_datetime is not None:
            self.end_datetime = end_datetime
                 
        #Raise an error if no filenames are set 
        #(may have been set not as part of readdata keywords)
        if self.filenames is None:
            raise ValueError('No filenames given')
        
        #Now set up iris constraints:
        constraints = None

        if self.forecast_day is not None:
            fcst_day_constraint = iris.Constraint(forecast_day
                                                  = self.forecast_day)
            constraints = constraints & fcst_day_constraint

        if self.short_name_list is not None:
            sname_constraint = iris.AttributeConstraint(short_name=lambda c:
                                                      c in self.short_name_list)
            constraints = constraints & sname_constraint

        if self.start_datetime is not None:
            time_constraint = iris.Constraint(time=lambda c:
                                              c.bound[0]>= self.start_datetime)
            constraints = constraints & time_constraint
            
        if self.end_datetime is not None:
            time_constraint = iris.Constraint(time=lambda c:
                                              c.bound[1]<= self.end_datetime)
            constraints = constraints & time_constraint
        
        #Finally read in data using iris load:
        
        self.gridded_cube_list = iris.load(self.filenames,
                                           constraints,
                                           callback=self.__callback )

        return self.gridded_cube_list
            
__callback()
------------

To read data into an iris cube, often a callback is required. Callbacks can also
be useful for setting attributes or cube names from which constraints can be based
on (for example the short_name_list) - this makes the reading much more efficient.
Callbacks can also be used to set attributes or coordinates based on filenames - 
for example ensemble member numbers. Note the callback is not expected to be used
outside of this module, so prefix it with __ to make it private.

.. code-block:: python
   
    def __callback(self, cube, field, filename):
    	"""
        Private method to provide callback to iris.load.
        Adds standard_name and short_name if known.
        Adds label attribute.
        """
        
        #Use conversion dictionaries to setup cube name
        # and short_name.
        cubename = cube.name()
        if cubename in NEWFORMAT_2_STDNAME:
            cube.rename(NEWFORMAT_2_STDNAME[cubename])
        if cubename in NEWFORMAT_2_SHORTNAME:
            cube.attributes['short_name'] = NEWFORMAT_2_SHORTNAME[cubename]
        else:
            cube.attributes['short_name'] = cubename  
        
        #Give label to cube
        cube.attributes['label'] = self.label      

.. _dataclass_extract_sites:

extract_sites()
---------------

Now gridded data should be being read in successfully. It is now a simple 
process to extract site-specific data, assuming sites data is given in 
:class:`sites_info.SitesInfo` format. To do this all that needs to be known is the X and Y
coordinate axis names. This would generally be 'grid_longitude' and 'grid_latitude'
or 'longitude' and 'latitude', or something similar. Here we show an example using
longitude and latitude. We make use of :meth:`adaq_data.ADAQData.extract_scl_from_gridded`
which is defined in the parent class :class:`adaq_data.ADAQData`.

.. code-block:: python

    def extract_sites(self, sites_data=None):
        """
        Extract site-specific data from gridded_cube_list into
        sites_cube_list, given
        
        :param sites_data: site information dictionary from :class:`sites_info.SitesInfo`
        """

        if self.sites_data is None:
            self.sites_data = sites_data
        if sites_data is None:
            if self.sites_data is None:
                raise ValueError("No sites requested")
            else:
                sites_data = self.sites_data     

        self.sites_cube_list = self.extract_scl_from_gridded(sites_data,
                                                             'longitude',
                                                             'latitude')
        
        return self.sites_cube_list                                                     

doctests
--------

Finally, write some doctests (:ref:`code_testing_doctests`), especially in the example section at the beginning
of the class, and ensure the doctests get tested:

.. code-block:: python

    if __name__ == '__main__':

        import doctest
        doctest.testmod()


Final notes
-----------
   
This is essentially all that is required in the class. There may be other functions 
specific to this data type, for example to convert unusual time formats or for unit
conversion. There may also be, either as part of this class, or as another class in the
same file, methods such as get_filenames (eg :meth:`nimrod_data.NimrodData.get_filenames`)
and an accompanying get_fileinfo (eg :meth:`nimrod_data.NimrodData.get_fileinfo`). These
are useful when the filenames are of a common format such that they can be used to 
determine the dates within the file etc.
These methods can then be used to limit the filenames down on the basis of a start and end 
date, or forecast period. Iris is much more efficient at readding in data if it has a 
shorter list of filenames to read from rather than just giving it an entire directory.
               
          
   
