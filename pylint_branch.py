#!/usr/bin/env python
"""
Script to check for all python code modified as part of this branch
and give pylint report for this code and a summary report of scores

usage: ./pylint_branch.py
"""
from __future__ import print_function

import os.path
import subprocess

def run_pylint_branch(command=None):
    """
    Run pylint on all code modified as part of this branch.
    If command keyword is set then this should provide a unix command
    to return a list of files to be tested, eg
    command = 'ls adaqcode/*.py adaqscripts/*.py' to check all python code
    in these two directories.
    """
    summary = []
    rcfile = os.path.dirname(os.path.realpath(__file__)) + '/pylintrc'

    #Set up command to determine which .py files have been modified in branch
    if command is None:
        command = 'fcm diff -b | grep "^Index: " | cut -c8- | grep ".py"'

    p1 = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)

    for pyfile in p1.stdout:
        pyfile = pyfile.decode().strip()
        print('*******************************************')
        print('Pylint report for ', pyfile)

        #Extract directory name - pylint should be run from same directory
        #as code to ensure imported modules can be found.
        dirname = os.path.dirname(pyfile)
        if not dirname:
            dirname = './'
            
        filename = os.path.basename(pyfile)

        command = ('cd ' + dirname + '; pylint --rcfile=' +
                   rcfile + ' ' + filename)
        p2 = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)

        for out in p2.stdout:
            outstr = out.decode()
            print(outstr)
            #Also get score for summary report
            if outstr[:27] == 'Your code has been rated at':
                summary.append(pyfile+': '+outstr.split()[6].split('/')[0])

    #Print out summary
    print('*******************************************')
    print('*******************************************')
    print('Summary:')
    for line in summary:
        print(line)


if __name__ == '__main__':
    run_pylint_branch()
