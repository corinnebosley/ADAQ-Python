ADAQ Sandpit 
============

Introduction
^^^^^^^^^^^^

The ADAQ sandpit is a space to put code which is useful to or may be useful to more than
one member of ADAQ but which has not yet been integrated into the full ADAQ Python Library.

Scripts placed in this directory do not need to be code reviewed but users should follow 
the instructions given below before placing code in this directory.

This does mean that code in this directory is not fully tested, may be accidently broken
and should **not** be used operationally.

Instructions for Use
^^^^^^^^^^^^^^^^^^^^

#. Create a ticket (optional):

   * In the Trac pages, click on 'new ticket'.
   * Give a brief summary of the ticket and a more detailed description
   * Set Type and Component and an optional Milestone and Keywords
   * Generally set Owner to your own username
   * Can also cc other people so they get emails when the ticket is updated
   * 'Create ticket'
 
#. Create a branch 

   * *If setting up for the first time, firstly add the following line into your ~/.fcm file:*
     
     *set::url::adaq_python svn://fcm8/ADAQ_PythonCode_svn/ADAQ_PythonCode*
     
   * .. code-block:: ksh
   
       fcm branch-create $Branch_name --rev-flag NONE fcm:adaq_python
   
   * Where $Branch_name is a name for your branch 

#. Checkout branch:

   * .. code-block:: ksh
   
       cd $DATADIR #Or another suitable directory
       fcm checkout fcm:adaq_python_br/dev/$userid/$Branch_name 
   
   * This will put the code into $DATADIR/$Branch_name    

#. Modify or add some code

   * Add the code to the *adaqsandpit* directory 
   * Add some documentation at the beginning of the code so other users know what
     the code does and how to run it    

#. If a new file is being added, add the code to fcm:

   * .. code-block:: ksh
      
      fcm add newcode.py

#. Update documentation for this code:

   * Open the file *adaqsandpit/contents.rst* and add some text which looks like:
   * .. code-block:: text

      My new code
      -----------

      .. automodule:: newcode
             :members:

#. Check documentation is OK - this can be done by :ref:`building_docs`. The documentation can
   be seen by opening _build/html/index.html in a web browser.

#. Commit code to branch:

   * .. code-block:: ksh
      
      fcm commit

#. Merge in latest version of trunk into branch:
   
   * .. code-block:: ksh
      
      cd $branch_name
      fcm merge fcm:adaq_python_tr
   
   * Resolve any conflicts:
     
     .. code-block:: ksh
      
      fcm conflicts
   
   * When happy, re-commit to branch:

     .. code-block:: ksh
      
      fcm commit

#. Now can commit to trunk:
      
   *  .. code-block:: ksh
       
       fcm checkout fcm:adaq_python_tr
       cd trunk
       fcm merge fcm:adaq_python_br/dev/$user/$branch_name
       fcm commit

