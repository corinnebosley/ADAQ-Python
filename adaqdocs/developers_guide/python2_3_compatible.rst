Writing python2 and python3 compatible code
===========================================

.. _python2and3_compatible:

ADAQpython is now able to be run under both python2 and python3. Any new code
therefore needs to keep this ability. When writing new code, some hints on how
to achieve this are given below. There is also a 
`cheat sheet <http://python-future.org/compatible_idioms.html>`_ which has extra tips.

Print function
--------------
When printing, use print as a function -  note this requires an extra import
statement first to ensure it works under python2 if printing more than a
single statement:

Single statement example:

.. code-block:: python

  print('hello')
  
Multiple statements:

.. code-block:: python

  from __future__ import print_function
  print('hello', 'world')
  
In the case of strings/floats, the print format can be variable between different
versions (for example whitespace and precision). 
To overcome this, use `format strings <https://docs.python.org/3.4/library/string.html#formatspec>`_ .
For example:

.. code-block:: python

  fmtstr = '{:>6},{:>4},{:8.3f},{:7.3f},{:4d},{:>8},{:>16}'
  print(fmtstr.format(*site))
  
.. code-block:: text

  GB0002,  AH,  52.504, -3.034, 370,   RURAL,      Aston_Hill

In the case of NumPy arrays, it is advised to use the function
`np.set_printoptions() <https://docs.scipy.org/doc/numpy-1.13.0/reference/generated/numpy.set_printoptions.html>`_ .
For example:

.. code-block:: python

  import numpy as np
  np.set_printoptions(formatter={'float':lambda x: '{:5.2f}'.format(x)})
  print(cube.data)
  
.. code-block:: text

  [16.00 20.00 37.00 39.00 31.00 40.00 37.00 23.00]


at the end, remember to undo the settings ...

.. code-block:: python

  np.set_printoptions()
  
Division
--------
  
Division in python3 is always as a float. If integer division is required, then
use it explicitly:

.. code-block:: python
  
   2 // 3 #=0
  
Octals
------

When using integers, ensure they do not start with a 0 as these will
be interpreted as an octal instead. This is particularly noticable for 
datetimes:

.. code-block:: python

  datetime.datetime(2018,08,06) #Don't use this
  datetime.datetime(2018,8,6) #Use this instead

np.genfromtxt
-------------
     
numpy.genfromtxt causes lots of issues with strings verses bytes. 
The preferred method for new code is to use 
`pandas.read_table <https://pandas.pydata.org/pandas-docs/stable/generated/pandas.read_table.html>`_

.. code-block:: python

  table = pd.read_table(filename, sep='\s+')
  
If the header names need changing then this can be done easily using eg:

.. code-block:: python

  columns_in = table.columns
  columns_out = {c: c.lower() for c in table.columns}
  for col_in, col_out in columns_out.items():
      if col_out == 'lat':
          col_out = 'latitude'
  columns_out[col_in] = col_out
  table.rename(columns=columns_out, inplace=True)
  
If this data is still required in the same format as would be read in using
np.genfromtxt (a numpy record array), then use:

.. code-block:: python

  data = np.array(table.to_records(index=False))

Dictionaries
------------

Looping over keys:

.. code-block:: python

  for key in mydict:
      print(key)
     
Looping over key and values:

.. code-block:: python

   for key, value in mydict.items():
       print(key, value)
      
Dictionaries are not guaranteed to be printed in any particular order. 
This can be resolved by sorting the loop over the dictionary keys alphabetically: 

.. code-block:: python
    
   for key in sorted(mydict):
       print key, mydict[key]

      
Range
-----

range() is now an object in python3. This is fine for using in loops, but
if required as a list, then need to specifically convert to a list:

.. code-block:: python

  for i in range(3):
    print(i)
    
.. code-block:: python

  myrange = list(range(3))
