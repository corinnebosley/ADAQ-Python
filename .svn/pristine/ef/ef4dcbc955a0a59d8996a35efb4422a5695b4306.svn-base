"""
Functions that are used to call and interact with unix shell commands.
"""
from __future__ import print_function

from six.moves.builtins import range
import subprocess
import time

def check_mass_dir(massdir):
    """
    Check if mass directory exists

    :param massdir: string, directory name to check on mass

    :returns: True, False or 'Unknown'.
              If True then directory exists on mass
              If False then directory does not exist on mass
              If 'Unknown' then there has been an error with checking,
              probably due to mass being unavailable.

    Note this cannot sensibly be tested with a doctest as this relies on the
    user having access to mass.
    """

    command = 'moo test -d '+massdir

    process = subprocess.Popen(command, shell=True,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)
    #Print standard output and standard error as command runs
    for line in iter(process.stdout.readline, b''):
        result = line.decode().strip() #Print as ascii
        if result == 'true':
            print('Directory exists')
            return True
        elif result == 'false':
            print('Directory does not exist')
            return False

    #If something else, then check return code:
    returncode = process.wait()
    print('Not known if directory exists')
    print('Return code from moo test:', returncode)
    return 'Unknown'

def call_mass(selectfile, massdir, outputdir, file_patterns=None,
              massretries=0, massretrydelay=60, retrieve=True):
    """
    Set up moo select or moo get shell command and call it using shell.
    If any errors with mass, will keep retrying based on massretries
    and mass retrydelay.

    :param selectfile: filename of file containing moo select query
    :param massdir: directory(s) or tar file on mass to retrieve from
    :param outputdir: directory on local machine to place retrieved files
    :param file_patterns: list of file name patterns for using with moo get
                          if selectfile=None
    :param massretries: number of times to retry mass retrieval
    :param massretrydelay: sleep time in seconds between each retry
    :param retrieve: Logical. If set to True retrieves files from mass.
                     If False, sets up all required moo select files and returns
                     the command to be run manually by user.
    :returns: String, containing command submitted to shell

    Example of running command, but with retrieve set to False so doesn't
    actually do retrieval as this is likely to fail if user does not have
    mass access or mass is currently not working!:

    >>> command = call_mass('selectfilename.txt', 'moose:/devfc/suiteid/',
    ... '/output/file/dir', retrieve=False)
    Moose Retrieval Command:
    moo select -f selectfilename.txt moose:/devfc/suiteid/ /output/file/dir

    >>> print(command)
    moo select -f selectfilename.txt moose:/devfc/suiteid/ /output/file/dir

    """

    if selectfile is None:
        if file_patterns is None:
            cmd = 'moo get -f ' + massdir + ' ' + outputdir
        else:
            moo_list = [massdir + '/' + p for p in file_patterns]
            cmd = 'moo get -f '+ ' '.join(moo_list) + ' ' + outputdir
    else:
        cmd = 'moo select -f ' + selectfile + ' ' + massdir + ' ' + outputdir
    print('Moose Retrieval Command:')
    print(cmd)
    if retrieve:
        for retry in range(-1, int(massretries)):
            if retry > 0:
                print('Retry number ', retry)
            #Set shell command going
            process = subprocess.Popen(cmd, shell=True,
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.STDOUT)
            #Print standard output and standard error as command runs
            for line in iter(process.stdout.readline, b''):
                print(line.decode()) #Print as ascii
            #Wait for command to finish running and get return code to check
            returncode = process.wait()
            print('Return code:', returncode)
            if returncode == 0:
                break
            elif retry < int(massretries)-1:
                #Non-zero return code, so try again after a period of time
                print("Error retrieving from mass, sleeping for ", \
                      massretrydelay, "seconds before retrying...")
                time.sleep(int(massretrydelay))

        if returncode != 0:
            #Fatal error
            raise IOError("Mass retrieval failed")
    return cmd

def call_shell(command):
    """
    Generic function to call shell using input command and print
    standard out. Waits for command to finish running.

    :param command: shell command

    :returns: returncode - command's return code.

    >>> returncode = call_shell('echo hello')
    hello
    <BLANKLINE>
    >>> print(returncode)
    0

    """

    process = subprocess.Popen(command, shell=True,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)
    #Print standard output and standard error as command runs
    for line in iter(process.stdout.readline, b''):
        print(line.decode()) #Print as ascii

    #Wait for command to finish running and get return code to check
    returncode = process.wait()

    return returncode

if __name__ == '__main__':

    import doctest
    doctest.testmod()
