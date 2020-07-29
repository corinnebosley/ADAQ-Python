.. _code_dev_guide:

Code Development Guide
======================

Initial Code Development
^^^^^^^^^^^^^^^^^^^^^^^^

1. Create a ticket:

 * In the Trac pages, click on 'new ticket'.
 * Give a brief summary of the ticket and a more detailed description
 * Set Type and Component and an optional Milestone and Keywords
 * Generally set Owner to your own username
 * Can also cc other people so they get emails when the ticket is updated
 * Click 'Create ticket'

2. Create a branch 

 * *If setting up for the first time, firstly add the following line into your ~/.fcm file:*

   *set::url::adaq_python svn://fcm8/ADAQ_PythonCode_svn/ADAQ_PythonCode*

 * .. code-block:: ksh

     fcm branch-create $Branch_name fcm:adaq_python --ticket=$N

 * Where $Branch_name is a name for your branch and $N is the ticket number you created above.  
 * This will create a branch rXXX_$Branch_name where rXXX is the revision number XXX of the current 
   version of the trunk that this has been branched from. 
   
3. Comment on ticket to give a link to the branch:

 * .. code-block:: trac-wiki

    Branch created: [source:ADAQ_PythonCode/branches/dev/userid/rXXX_Branch_name rXXX_Branch_name]

 * Where userid is your linux user id.
 * Also change status of ticket to 'In Progress' (first click on 'Modify Ticket').

4. Checkout branch:

 * .. code-block:: ksh

     cd $DATADIR #Or another suitable directory
     fcm checkout fcm:adaq_python_br/dev/$userid/rXXX_$Branch_name 

 * This will put the code into $DATADIR/rXXX_$Branch_name    

5. Modify or add some code. Remember to follow our :ref:`coding_standards` while developing code. 
   It is much easier to follow the coding standards from the start of the work, than
   to try and make the changes later. If adding code to read a new format of data, then
   see :ref:`adding_dataclass`.
   
 * For some developments you may need to prepare some new sample data to work with your new 
   doctests. See :ref:`code_testing_new_sample_data`. 

6. Note that when developing and using  ADAQ Python  at the Met Office you should use the latest version 
   of the scitools environment. To do this type: 

 * .. code-block:: ksh

     module load scitools

Preparing for review
^^^^^^^^^^^^^^^^^^^^   

6. Once happy with the code, add some doctests to provide examples of how to use the code.
   This should also provide tests to ensure nobody can break the code in the future.

 * For more information, see :ref:`code_testing_doctests`.
 * Also need to test you haven't broken any other doctests in the code by running the 
   whole module.
 * If adding doctests into a new module, the module needs to be added to 
   :data:`doctests.DOCTEST_MODULES` list.       

7. Code should be added or deleted from fcm as appropriate: 

 * If a new module (file) is being added, add the code to fcm:

   .. code-block:: ksh

         fcm add newcode.py

 * If a module (file) is being deleted, delete the code (and its .rst file) from fcm:

   .. code-block:: ksh

        fcm del oldcode.py
        fcm del oldcode.rst

 * If deleting files, also manually delete the .pyc file which may have been created:

   .. code-block:: ksh

        rm oldcode.pyc   
   
8. Update documentation for this code:

 * The documentation for each function is generated from the docstrings in the
   code. This means that the code and the webpages have exactly the same information.
   The system that creates the documentation (Sphinx) uses .rst files which have
   the same name as the .py file for the module to add each module into
   the functions being documented. 
 * If you are adding new functions to an existing module the .rst file should
   already exist. Ensure that your docstrings contain sufficient explanation of the purpose
   and arguments of the code. These docstrings can contain 
   formated text. For more details on the mark up used see :ref:`sphinx_markup`. 
   Also check the .rst file in case functions etc are defined separately 
   and new code needs to be added individually. 
   (If it uses the automodule approach below this should not be necessary).
 * If your changes alter the arguments and/or functionality of an existing
   function make sure that the docstrings are changed accordingly.
 * If it is a new module, first create a .rst file for this code 
   (same name as .py file, just with .rst suffix instead). This file should 
   contain the following (replace filename with the name of the python file
   to be documented).
 * .. code-block:: text

    filename.py
    ============

    .. toctree::
          :maxdepth: 2

    New Class
    ---------------

    .. automodule:: filename
           :members:

 * For more information, see :ref:`editingdoc_modulerst`.

 * Any new .rst file needs to be added to reference.rst (or if it
   is in adaqscripts, then add to adaqscripts/adaqscripts.rst). It also needs
   to be added to fcm:

   .. code-block:: ksh

    fcm add newcode.rst
     
 * An example docstring is below. This uses rst markup to highlight the input parameter
   (shell command) and the return value. It also has two examples of usage which
   also act as doctests.
 
 * .. code-block:: text

     """
    Generic function to call shell using input command and print
    standard out. Waits for command to finish running.

    :param command: shell command

    :returns: returncode - command's return code.

    >>> returncode = call_shell('echo hello')
    hello
    <BLANKLINE>
    >>> print returncode
    0

    """
    
9. Check documentation is OK 

 * This can be done by :ref:`building_docs`. The documentation can
   be seen by opening _build/html/index.html in a web browser:
   Note this should be done using the python3
   version of scitools (currently default-current).
   
   .. code-block:: ksh
     
        cd /location/of/branch
        cd adaqdocs
        module unload scitools
        module load scitools      
        ./build_html_docs.scr
        firefox ../_build/html/index.html

 * Check that any new functions are included in the documentation, 
   and that all formatting works properly.

10. Run pylint to check against coding standards:
   
 * See :ref:`pylint` for instructions and coding standards.
 * For code to be admitted to the trunk it needs to
   have a score of at least 5, but ideally with a much higher score - all errors that
   it complains about must be justifiable.
 * The scores for all modified modules must then be included on the ticket. This is done
   easily by running :mod:`pylint_branch`.py. Note this should be done using the python3
   version of scitools (currently default-current).

   .. code-block:: ksh

       cd /location/of/branch
       module unload scitools
       module load scitools
       ./pylint_branch.py     
 
 * The summary report from this should be copied & pasted onto the ticket.      

11. Once happy with code changes, run ./:mod:`doctests`.py to check that doctests in all 
    modules still work and nothing has been broken. Note: the doctests need to be run under
    both python2 and python3 versions of the Scientific Software Stack. The easiest way to
    do this is to run doctests_spice_rhel7.sh which is set up to run both of these automatically.
    To do this the code must be checked out or copied to a networked file location 
    i.e. NOT /data/local.

 * .. code-block:: ksh
 
     cd /location/of/branch
     sbatch doctests_spice_rhel7.sh
 
 * The output from this will be in the latest created file doctests-xxxxx.out - check the 
   results from this both from python3 (top of output file) and python2 (bottom of file).
 
 * Update the ticket to say that the doctests have all been run successfully.      
   
12. Commit code to branch:

 * .. code-block:: ksh

    fcm commit

 * State revision number on ticket:

 * .. code-block:: trac-wiki

     Changes committed at r150

 * r150 will become a hyperlink to the changes at that revision.  

13. Consider whether these changes are significant and may have 
    an impact on any of the emergency response systems. 
    If so, test your branch on these systems.

14. If it has been a while since you started development, (ie the trunk has changed 
    significantly) then update your branch to the latest version of the trunk before review. 
    (See above for guidance on these items)
   
 * Create a new branch at the latest revision of the trunk and check this out
   eg rYYY_$Branch_name and update ticket.
 * Merge in original branch:

   .. code-block:: ksh

     cd rYYY_$Branch_name
     fcm merge fcm:adaq_python_br/dev/userid/rXXX_$Branch_name

 * Again check doctests, pylint scores and documentation are all still OK and haven't been broken.
 * Commit code to branch and update ticket.      

15. Put a link to the branch in the ticket description to make it easier for the reviewer to find.  

  .. code-block:: trac-wiki

    Branch: [source:ADAQ_PythonCode/branches/dev/userid/rXXX_Branch_name rXXX_Branch_name]

  * If you undertook step 14 this will be

  .. code-block:: trac-wiki

    Branch: [source:ADAQ_PythonCode/branches/dev/userid/rYYY_Branch_name rYYY_Branch_name]

   
16. Pass code to someone else to review:
  
 * Comment on ticket that passing code to a particular person to review.
 * Under "Action" click "Assign for review to:" and put the user ID of the reviewer in 
   the box. The status of the ticket will change from  "in_progress" to "reviewing"
   and the reviewer will now own the ticket.

.. _code_dev_guide_reviewer:

Reviewer
^^^^^^^^
 
17. *Code Reviewer* - Review code:

 * Note that if you click on a revision number on the ticket this allows you
   to see the changes from the last revision highlighted. You may want to
   click on 'Revision Log' (top right hand corner) and then 'View Changes'
   to see all changes that have gone in.    
 * Add comment in ticket description:
 * .. code-block:: trac-wiki

     [wiki:ticket/N/CodeReview CodeReview]

 * Where N is the ticket number 
 * This will create a link to a new page - click it to get to page.
 * Change 'Using the template' to 'CodeReview', then click 'Create this page'
 * Checkout the latest version of the branch:

   .. code-block:: ksh

     cd $DATADIR #Or another suitable directory
     fcm checkout fcm:adaq_python_br/dev/$userid/rXXX_$Branch_name 
     module load scitools

 * Answer the questions on the code review template and give any further comments.
 * Consider whether these changes are significant and may have an impact on the emergency response systems.
   If so, check any testing that has been done or request it if you feel it's necessary.
 * Run ./:mod:`doctests`.py - this can pick up any 
   hard-coded directory locations. (For more information about doctests, see
   :ref:`code_testing_doctests`): Note: the doctests need to be run under
   both python2 and python3 versions of the Scientific Software Stack. The easiest way to
   do this is to run doctests_spice_rhel7.sh which is set up to run both of these automatically. 
   To do this the code must be checked out or copied to a networked file location 
   i.e. NOT /data/local.

 * .. code-block:: ksh

        cd /location/of/branch
        sbatch doctests_spice_rhel7.sh
 
 * The output from this will be in the latest created file doctests-xxxxxx.out - check the 
   results from this both from python3 (top of output file) and python2 (bottom of file).

 * Run ./:mod:`pylint_branch`.py - to check any errors flagged.
   (For more information about this, see :ref:`pylint`). Note this should be done using the python3
   version of scitools (currently default-current).

   .. code-block:: ksh

      cd /location/of/branch
      module unload scitools
      module load scitools
      ./pylint_branch.py

 * Build documentation (see :ref:`building_docs`) and read in web browser to check 
   new documenation makes sense. Note this should be done using the python3
   version of scitools (currently default-current).

   .. code-block:: ksh

      cd /location/of/branch
      cd adaqdocs
      module unload scitools
      module load scitools
      ./build_html_docs.scr
      firefox ../_build/html/index.html

 * If there are no problems, then can give approval to be committed to trunk. 
   To do this, click 'Modify Ticket' and then "approve & assign to" and put the 
   authors name in the box. The status will change to "approved".
 * If not, pass the ticket back to the code developed to correct problems -
   once these have been sorted and recommitted to branch, approval can be given. 
 * NB if the trunk has been merged in the branch, it may be easier to use the following
   command, run from within the checked out branch, to see the differences:

   ..  code-block:: ksh   

     fcm diff -g ^/ADAQ_PythonCode/trunk@xxx ^/ADAQ_PythonCode/branches/dev/userid/Branch_name@yyy

  where xxx is the version of the trunk when last merged into the branch 
  and yyy is the latest modified version of the branch.   

Finishing off
^^^^^^^^^^^^^
   
18. Merge in latest version of trunk into branch:
   
 * .. code-block:: ksh

    cd $branch_name
    fcm merge fcm:adaq_python_tr

 * Resolve any conflicts:

   .. code-block:: ksh

    fcm conflicts

 * Rerun doctests to ensure nothing has broken due to changes in the trunk:
 
   .. code-block:: ksh

        cd /location/of/branch
        sbatch doctests_spice_rhel7.sh 
 
 * If you created new sample data to test your changes on, you may have modified SAMPLE_DATADIR in 
   :mod:`config`.py to point to one of your own directories containing your new sample data. 
   Remember to change this back again to the usual location 
   '/data/users/apdg/python_sample_data/' 

 * When happy, re-commit to branch:

   .. code-block:: ksh

    fcm commit

19. Now can commit to trunk:
      
 *  .. code-block:: ksh

     cd $DATADIR #Or another suitable directory
     fcm checkout fcm:adaq_python_tr
     cd trunk
     fcm merge fcm:adaq_python_br/dev/$user/$branch_name
     fcm commit

 * Start the commit message with #N: where N is the ticket number. Give a brief overview of changes
   - probably just the ticket summary.    

 * State revision number of commit to trunk on ticket:

   .. code-block:: trac-wiki

     Changes committed to trunk at r150


20. Update the checked out copy of the code stored under apdg

 *  .. code-block:: ksh

     ssh -Y apdg@vld???
     cd /home/h03/apdg/PythonCode/ADAQ_Python/trunk
     fcm update

 * Replace vld??? with the name of a desktop computer

21. If documentation has been modified, this should be updated in the main html pages.
    Note this should be done using the python3 version of scitools (currently default-current).
    Documentation can be updated:
    
 *  .. code-block:: ksh

     ssh -Y apdg@vld???
     cd /home/h03/apdg/PythonCode/ADAQ_Python/trunk/adaqdocs
     module unload scitools
     module load scitools
     ./build_html_docs.scr ~apdg/public_html/adaq_pythoncode

 * Replace vld??? with the name of a desktop computer

22. If you modified or created new sample data, ensure this is copied to the 
    usual location '/data/users/apdg/python_sample_data/' specified in :mod:`config`.py.
    You should also copy the new sample data to MASS. See :ref:`code_testing_sample_data_to_mass`.

23. Mark ticket as resolved - fixed. This will also close the ticket.

       
