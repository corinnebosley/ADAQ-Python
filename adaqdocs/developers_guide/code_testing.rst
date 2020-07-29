Testing Code
============
.. _code_testing_doctests:

Doctests
^^^^^^^^

It is important to ensure that all new code that is written is testable and will continue
to be tested as part of the ADAQPython code base. For this, we make use of doctests.
These supply a simple python test within the documentation string, which is 
then immediately followed by the expected answer. When running the doctests, 
the code is run and checked against the expected answers.

Within the documentation strings (enclosed within triple quotes), code can be 
included as examples, by typing the code following >>>. For example:

.. code-block:: python

  """
  >>> a, b = 3, 4
  """

In particular, expected results can be printed, by including them on line immediately
following the code line, on a plain line, eg:

.. code-block:: python

  """
  >>> print a, b, a+b
  3 4 7
  """
  
For dictionaries, as the printing order can not be guaranteed, either:

.. code-block:: python

   """
   >>> mydict == {'a':1, 'c':2, 'b':3}
   True
   """
   
OR:

.. code-block:: python

   """
   >>> for key in sorted(mydict.keys()):
   ...    print key, mydict[key]
   a 1
   b 3
   c 2
   """
   

These examples can then be tested using

.. code-block:: python

 if __name__ == "__main__":
     import doctest
     doctest.testmod()     

This will check all code segments give expected answers. It will raise errors
if it finds any answers do not match the expected answer.

Alternatively, while developing, to just test a single function/method within the module:

.. code-block:: python

 if __name__ == "__main__":
     import doctest
     doctest.run_docstring_examples(mymethod, globals()) 


If any of the code outputs are pointing eg at a particular directory which
may change depending on user, then can use ellipsis (three dots) to represent 
this section (note similar to ls * it may match other incorrect answers unless
used sensibly). However doctest needs to be informed of its usage:

.. code-block:: python

  """
  import config
  print config.SAMPLE_DATADIR+'AURN_obs' # doctest: +ELLIPSIS
  /.../AURN_obs
  """
  
If answers get split over several lines, then to keep doc strings tidy, 
may be useful to use NORMALIZE_WHITESPACE which treats all sequences 
of whitespace (blanks and newlines) as equal: 

.. code-block:: python

  """
  >>> print x # doctest: +NORMALIZE_WHITESPACE
  [  0.00000000e+00   1.74929180e-07   6.45126085e-05   2.01771603e-03
  3.48425356e-02   2.69860042e-01   6.44658106e-01   9.22545804e-01
  9.93557700e-01   9.99683948e-01   9.99841974e-01   1.00000000e+00]
   """
   
If code needs to be split over multiple lines, for example using for loops:

.. code-block:: python

  """
  >>> for i in range(5):
  ...   if i > 2:
  ...      print i   
  """
  
doctests.py
^^^^^^^^^^^  

A script to run all known doctests can be used: :mod:`doctests`.py - note that 
to use this at the Met Office, you first need to type 

.. code-block:: ksh

    module load scitools
    
The doctests should be run under both python2 and python3 versions of scitools
to check compatability with both.

You can also run these tests on SPICE for faster turn around and lower
impact on your machine. This will also automatically check both python2 and python3
versions of scitools. NB to do this the code must be checked out to a networked
file location i.e. NOT /data/local. SPICE CAN'T SEE YOUR /data/local!

.. code-block:: ksh

    sbatch doctests_spice_rhel7.sh

This runs the code on spice using the SLURM queuing system. It will produce
a file doctests-xxxxxx.out where xxxxxx is the ID of the job run in 
batch on SPICE.

You can see if you job is still running with the command:

.. code-block:: ksh

    squeue -u $USER

You can also set up the batch script to email you when it is finished. See
the SPICE user guide for more information.

.. _code_testing_new_sample_data:


Creating new sample data
^^^^^^^^^^^^^^^^^^^^^^^^
Some new developments may require the generation of new sample data files, e.g. to 
be used in doctests. The suggested approach is to create directory under your own 
account which links to the existing sample data (under apdg) but also contains
the new sample data you require.

Create the directory for your sample data:

* .. code-block:: ksh

    mkdir /data/users/my_user_id/python_sample_data

and then make a symbolic link to the data files at '/data/users/apdg/python_sample_data':

* .. code-block:: ksh

    cd /data/users/my_user_id/python_sample_data
    ln -s /data/users/apdg/python_sample_data/* .


For large air quality files it may be helpful to create new sample data files using :mod:`adaq_functions`.py. 
You can run this module as 'main', calling the method 
:meth:`adaq_functions.create_example_data`. Examples are given at the bottom of 
:mod:`adaq_functions`.py. Modify the code or the files created in the way that you need for your 
sample data. The method reads model pp output (the location of which is specified in 
the relevant ini file, e.g. 'example_data_3d.ini'). New netcdf files are created and these 
are the sample data files which can be used in doctests (the netcdf files are much faster 
to load than the pp files).

The sample data files created are placed in the following directories:

   ..adaqcode/gridded_cube_list

   ..adaqcode/sites_cube_list
   
These files can then be manually copied to the desired location 
(under /data/users/my_user_id/python_sample_data) and the ADAQ code modified to point to the new 
sample data. The 'SAMPLE_DATADIR' is specified in 
:mod:`config`.py. and this should be modified to point to '/data/users/my_user_id/python_sample_data'
This gives the top level directory for sample data. The actual 
directory where netcdf sample data are loaded from is specified 
in the method :meth:`adaq_functions.get_exampledata`.py. 

.. _code_testing_sample_data_to_mass:

Saving the sample data to MASS
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When the code using the sample data is committed to the trunk a copy of the sample data should be
saved to MASS. The sample data is stored in MASS at:

   moose:/adhoc/users/apdg/python_sample_data/

There is a tar file for each of the subfolders in the python_sample_data directory.

To add a new tar file you will first need to make sure that apdg is pointing to the correct
user name on MASS. Login as apdg then:

* .. code-block:: ksh

   cd .moosedir
   cp moose_apdg moose

Then tar the directory in the python_sample_data folder:

* .. code-block:: ksh

   cd /data/users/apdg/python_sample_data
   tar -cvf <my_directory>.tar <my_directory>

Then copy to MASS

* .. code-block:: ksh

   moo put <my_directory>.tar moose:/adhoc/users/apdg/python_sample_data

Finally delete the tar file

* .. code-block:: ksh

   rm <my_directory>.tar

.. note::

   If you are adding to an existing directory in python_sample_data then the tar
   file for the directory may already exist on MASS. It will be necessary to delete
   the tar file on MASS first using *moo rm <tarfile>* before following the instructions
   above.
   
   
