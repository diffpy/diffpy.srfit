#!/usr/bin/env python
##############################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2010 The Trustees of Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Chris Farrow
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
##############################################################################

"""Input utilities."""

__all__ = ["inputToString"]

import os.path

def inputToString(inpt):
    """Convert input from various modes to a string.

    This is useful when you want a method to accept a string, open file object
    or file name.

    inpt    --  An open file-like object, name of a file
                or a string containing the input.

    Returns the input in a string
    Raises IOError if the input is supected to be a file name, but the file
    cannot be found.

    """
    # Get the input into a string
    inptstr = ""
    if hasattr(inpt, "read"):
        inptstr = inpt.read()
    # FIXME check for typos in the file name
    elif os.path.exists(inpt) or (len(inpt) < 80 and inpt.count("\n") == 0):
        with open(inpt, 'r') as infile:
            inptstr = infile.read()
    else:
        inptstr = inpt

    return inptstr

# End of file
