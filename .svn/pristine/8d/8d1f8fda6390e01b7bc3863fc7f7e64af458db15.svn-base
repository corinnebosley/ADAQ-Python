#!/usr/bin/env python
"""
Utility to concatenate a contiguous sequence of AQUM LBC files, discarding time
points at the end of each file except the last. The time at the end of each 
file is expected to match the time at the start of the next.

Usage: aqum_lbc_cat lbc_file_1[,lbc_file_2][...] output_lbc_file

Note: This script is intended for LBC input files from UM vn10.1+.
Before that the integer constant 'num_times' in the header was used so would
need to be updated if concatenating files produced by older versions.
"""
from __future__ import print_function

import os
import sys
#import mule #Not currently in python3 software stack

STASH_OROGRAPHY = 31001


def get_args():

    """
    Check for valid command line arguments and return them.
    """

    if len(sys.argv) != 3:
        sys.exit('[FAIL] Expecting two args: <input file list> <output file>')
    infile_list = sys.argv[1].split(',')
    outfile = sys.argv[2]
    for infile in infile_list:
        if not os.path.isfile(infile):
            sys.exit('[FAIL] Input file {} does not exist'.format(infile))
    if os.path.isfile(outfile):
        print('[WARN] Output file {} already exists - overwriting'.format(
              outfile))
        print()
    return infile_list, outfile


def validity_time(field):

    """
    Return validity time of field as a formatted string.
    """

    year = field.lbyr
    month = field.lbmon
    day = field.lbdat
    hr = field.lbhr
    mins = field.lbmin
    secs = field.lbsec
    string = "{0:04d}-{1:02d}-{2:02d}T{3:02d}:{4:02d}:{5:02d}".format(
             year,month,day,hr,mins,secs)
    return string


def get_lbc_file(infile):

    """
    Load LBC file 'infile'. Return header data and valid times of first
    and last time-dependent fields. The first field in the file is expected to 
    hold the orography. The remaining fields are assumed to be time-dependent. 
    """

    lbc_in = mule.LBCFile.from_file(infile)
    # Correct grid type = should be LAM(no wrap), rotated
    if (lbc_in.fixed_length_header.horiz_grid_type == 0 or 
        lbc_in.fixed_length_header.horiz_grid_type == 3):
        lbc_in.fixed_length_header.horiz_grid_type = 103

    if lbc_in.fields[0].lbuser4 != STASH_OROGRAPHY:
        raise UserWarning(
            'File {} does not have orography in first field'.format(infile))
    t_first = validity_time(lbc_in.fields[1])
    t_last = validity_time(lbc_in.fields[-1])
    return lbc_in, t_first, t_last


def set_file_end_time(lbc_data, field):

    """
    Set last validity time in the fixed length header of 'lbc_data' to the 
    validity time of 'field'.
    """

    lbc_data.fixed_length_header.t2_year = field.lbyr 
    lbc_data.fixed_length_header.t2_month = field.lbmon
    lbc_data.fixed_length_header.t2_day = field.lbdat
    lbc_data.fixed_length_header.t2_hour = field.lbhr
    lbc_data.fixed_length_header.t2_minute = field.lbmin
    lbc_data.fixed_length_header.t2_second = field.lbsec


def main():

    """
    Concatenate LBC files.
    """

    infile_list, outfile = get_args()

    # Copy first file header for output 
  
    infile = infile_list[0]  
    lbc_in, t_first, t_last = get_lbc_file(infile)
    lbc_out = lbc_in.copy()

    # Reference the fields for output: includes orography and time-dependent
    # fields. The last time point will be discarded and replaced from the 
    # next input file which should be based on a more recent forecast.

    nf_tdep = lbc_in.integer_constants.num_field_types - 1
    new_fields = lbc_in.fields[:-nf_tdep]

    # Copy output fields from first file and get time-dependent fields from 
    # the remaining input files, copying all but the last time point from each.
    # The time of the first field should match that of the last field from 
    # the previous file. 

    infile_prev = infile

    for infile in infile_list[1:]:

        print('Copying fields from {}'.format(infile_prev))
        print('  First validity time {}'.format(t_first))
        print('  Last validity time {}, discarded'.format(t_last))
        lbc_out.fields.extend(new_fields)

        t_last_prev = t_last 
        lbc_in, t_first, t_last = get_lbc_file(infile)
        if (t_first != t_last_prev):
            raise UserWarning(
'Unexpected time in file {}: does not match last time in previous file'.format(
                              infile))
        new_fields = lbc_in.fields[1:-nf_tdep]

        infile_prev = infile

    # Copy of last file must include the last time point 

    print('Copying fields from {}'.format(infile_prev))
    print('  First validity time {}'.format(t_first))
    print('  Last validity time {}, retained'.format(t_last))

    lbc_out.fields.extend(new_fields)
    new_fields = lbc_in.fields[-nf_tdep:]
    lbc_out.fields.extend(new_fields)

    # Update header and write output file

    set_file_end_time(lbc_out, lbc_in.fields[-1])
    lbc_out.to_file(outfile)


if __name__ == '__main__':
    import six
    if six.PY3:
        raise ValueError('mule not currently available for python3 SSS')
    else:
        import mule
    main()


