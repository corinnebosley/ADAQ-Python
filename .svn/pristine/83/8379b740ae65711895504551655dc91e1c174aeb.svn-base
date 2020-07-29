Coding Standards
================

.. _coding_standards:

Standard Rules
^^^^^^^^^^^^^^

 * Indenation. Use 4 space indentation. Tabs must not be used.
 * Limit all lines to a maximum of 79 characters. The preferred way of wrapping long lines is by using Python's implied line continuation inside parentheses, brackets and braces. Long lines can be broken over multiple lines by wrapping expressions in parentheses. These should be used in preference to using a backslash for line continuation.
 * Separate top-level function and class definitions with two blank lines.
 * Use blank lines in functions, sparingly, to indicate logical sections.
 * Imports should usually be on separate lines.
 * Avoid extraneous whitespace.
 * Always surround these binary operators with a single space on either side: assignment (=), augmented assignment (+=, -= etc.), comparisons (==, <, >, !=, <>, <=, >=, in, not in, is, is not), Booleans (and, or, not).
 * Don't use spaces around the = sign when used to indicate a keyword argument.
 * Comments that contradict the code are worse than no comments. Always make a priority of keeping the comments up-to-date when the code changes!
 * Comments should be complete sentences.
 * Use inline comments sparingly.
 * Write docstrings for all public modules, functions, classes, and methods.
   
   *  `PEP 257 <http://legacy.python.org/dev/peps/pep-0257/>`_  describes good docstring conventions.
 
 * Never use the characters 'l' (lowercase letter el), 'O' (uppercase letter oh), or 'I' (uppercase letter eye) as single character variable names.
 * Modules should have short, all-lowercase names. Underscores can be used in the module name if it improves readability. Python packages should also have short, all-lowercase names, although the use of underscores is discouraged.
 * Class names should normally use the CapWords convention.
 * Function names should be lowercase.
 * Constants are usually defined on a module level and written in all capital letters with underscores separating words. 

Additional rules from Science IT
--------------------------------

 * The use of the `Doctest module <https://docs.python.org/2/library/doctest.html>`_ to incorporate small unit tests into your code is strongly encouraged.
 
   * See :ref:`code_testing_doctests` for more information on usage in the ADAQ Python code repository.
   
 * Docstrings should be as per PEP 8, with the addition that restructured text (reST) directives recognised by the 
   `Sphinx utility <http://sphinx-doc.org/>`_ may be used to apply additional formatting to module, function, class and method docstrings and to create a set of webpages. 

   * See :ref:`sphinx_markup` for more information on usage in the ADAQ Python code repository.

Recommendations for File Organisation
-------------------------------------

In order to make information, constants, functions etc. easier to find in files it is recommended that developers adhere to the following rules.

 * Constants and dictionaries should be placed at the top of each file immediately after any import statements.
 * Functions and methods should be listed in alphabetical order unless there is a strong logical reason for using a different order. For example if methods are always used in a set order then they could be placed in that order in the file.

Import orders
-------------

Modules should be imported into code in the following order groups - within
each group, they should be imported alphabetically with 'from' at the end:

 1. Any module related purely to python2 and python3 compatibility:
 
  .. code-block:: python
      
      from __future__ import division
      from __future__ import print_function
      from six.moves.builtins import str
      
 2. Standard python library imports
 
  .. code-block:: python
     
      import os
      import sys
      import warnings
      from datetime import datetime, timedelta
 
 3. Scientific Software Stack imports
 
  .. code-block:: python
  
      import cf_units
      import iris
      import numpy as np
      import matplotlib.pyplot os plt
      
 4. Local imports from adaqcode/adaqscripts:
  
  .. code-block:: python

      import config
      import sites_info
      
Note running these under python3's pylint should generally give no errors, however running under 
python2 will raise errors - for this reason pylint should always be run under python3.

.. _pylint:

Pylint
^^^^^^

Pylint is a useful tool to check if your python code corresponds to the required coding standards.
For use in the ADAQ python repository, there is a pylintrc file in the top level directory which sets 
up the required standards. This will check your code against pep8 and give it a 
score from minus infinity to 10. For code to be admitted to the trunk it needs to
have a score of at least 5, but ideally with a much higher score - all errors that
it complains about must be justifiable.

Single file
-----------

To test a single python file using pylint:

This should be done using python3 - the python2 version complains about issues
that are caused by using python3 coding standards (although some of these can be useful). 
So first load the python3 version of the Scientific Software Stack 
(currently default-current):

.. code-block:: ksh

   module unload scitools
   module load scitools

pylint can then be run as follows:

.. code-block:: ksh

    pylint pythonfile.py


pylint_branch.py
^^^^^^^^^^^^^^^^

To run pylint on all python code that has been modified as part of a branch, use
:mod:`pylint_branch`.py and run on the command line. Note this should be done using the python3
version of scitools (currently default-current) - the python2 version complains about issues
that are caused by using python3 coding standards.

.. code-block:: ksh

   module unload scitools
   module load scitools
   ./pylint_branch.py

This will give a pylint report for all modified python code, plus a summary report of the scores, which
can then be copied and pasted into the ticket.
