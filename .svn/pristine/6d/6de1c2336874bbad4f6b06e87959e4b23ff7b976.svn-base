Documentation
=============
   
Documentation for ADAQ Python Code is built using `Sphinx <http://sphinx-doc.org/>`_
which creates the documentation directly from docstrings (comments at the 
beginning of code enclosed within triple quotes) in the code. This ensures
that documentation can be easily kept up-to-date. The use of :ref:`code_testing_doctests` 
within the docstrings also allows for code to be easily tested and 
for examples to always be working. Further documentation is included in .rst files for
extra information not directly linked to code.

.. _building_docs:

Building html documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Simple solution:

.. code-block:: ksh
  
   ./build_html_docs.scr [$htmldirectory]

where $htmldirectory is optional - this is the location of where
you want the html to appear (for example ~/public_html/python_docs)
If this is not set, then you can immediately open the top level page by
doing:
.. code-block:: ksh
   
   cd ../_build/html  #Or $htmldirectory if set
   firefox index.html

Alternatively:

.. code-block:: ksh

  CODEDIR=/path/to/adaqcode
  export PYTHONPATH=$PYTHONPATH:$CODEDIR
  make html
  copy entire _build directory into public_html directory

To add a new module (eg newmodule.py) into the documentation:

 * Add a newmodule.rst file - see module_data.rst for example - change references to the new module
 * Note this also needs to be added to fcm
 * Modify reference.rst to include the new module
 * Also need to add newmodule into __init__.py
 * If doctests are included, also need to add module into :data:`doctests.DOCTEST_MODULES`

.. _sphinx_markup:
 
Sphinx markup
^^^^^^^^^^^^^

Given below are some useful examples of the markup used for sphinx. For further instructions, 
see the `Sphinx documentation <http://sphinx-doc.org/contents.html>`_  and the 
`reStructuredText Primer <http://sphinx-doc.org/rest.html>`_ .

General
-------

For bullet points:

.. code-block:: python

  """
    *  For bullet points 
  """

For numbered lists, eg.

#. Item A
#. Item B

.. code-block:: python

  """
    #. Item A
    #. Item B
  """

**Bold text:**

.. code-block:: python

 """
 **Bold**
 """

*Italic text:*

.. code-block:: python

 """
 *Italic*
 """

For code segments:

.. code-block:: python

  """
  >>> This is my code
  output_of_code
  """

To turn text into a heading:

.. code-block:: python

  """
  My heading
  ==========

  My subheading
  -------------
  """

For parameters/keywords:
    
.. code-block:: python

  """
  :param filenames: list of files to read
  """

Formulae, eg inline: :math:`x_{amp} = inversef_{oc}(f_{mc}(x_{mp}))` 

.. code-block:: python

   """
   :math:`x_{amp} = inversef_{oc}(f_{mc}(x_{mp}))` 
   """
   

Example code:

* eg python:

.. code-block:: python

  A=[1,2,3,4]
  for a in A:
    print A

.. code-block:: python

  """
  .. code-block:: python

    A=[1,2,3,4]
    for a in A:
      print a
  """

* eg. ksh:

.. code-block:: ksh

  for a in 1 2 3 4 
  do
    echo a
  done  
    
.. code-block:: python

  """
  .. code-block:: ksh

    A=1,2,3,4
    for a in A
    do
      echo a
    done  
  """      

Referencing
-----------

To give a link to an external url, eg `Sphinx utility <http://sphinx-doc.org/>`_ :

.. code-block:: python
  
  """
  `Sphinx utility <http://sphinx-doc.org/>`_ 
  """

For referencing other modules/classes, eg
site information data from a class, eg :class:`sites_info.SitesInfo`:

.. code-block:: python
  
  """
  site information from :class:`sites_info.SitesInfo`
  """

From a method, eg :meth:`sites_info.SitesInfo.read_from_file`:

.. code-block:: python

  """
  get site information data from :meth:`sites_info.SitesInfo.read_from_file`
  """

From a module, eg :mod:`sites_info`:

.. code-block:: python

  """
  site information data from :mod:`sites_info`
  """


To reference a section of documentation, eg to My Section:

First set up a label for it:

.. code-block:: text
   
   .. _my_section:
   
   My Section
   ==========
   
Can then refer to this elsewhere:
   
.. code-block:: text

   :ref:`my_section`   

.. _editingdoc_modulerst:

Modifying module .rst files
---------------------------  

To document the entire module, add the following to its rst file:

.. code-block:: python

  """
  .. automodule:: modulename
    :members:
  """

To just document a class:

.. code-block:: python

  """
  .. autoclass:: modulename.ClassName
    :members:
  """

To just document a function:

.. code-block:: python

  """
  .. autoclass:: modulename.function
  """

To document a global variable, firstly add comments for the variable 
using #: Comments in the line above the variable, then:

.. code-block:: python

  #: Set variable to 5.
  VARIABLE=5

  """
  .. autodata:: modulename.VARIABLE
  """






