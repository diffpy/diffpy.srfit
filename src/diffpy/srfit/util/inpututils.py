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


def inputToString(input):
    """Convert input from various modes to a string.

    This is useful when you want a method to accept a string, open file object
    or file name.

    input    --  An open file-like object, name of a file
                or a string containing the input.

    Returns the input in a string
    Raises IOError if the input is supected to be a file name, but the file
    cannot be found.
    """
    # Get the input into a string
    inptstr = ""
    if hasattr(input, "read"):
        inptstr = input.read()
    # TODO remove handling of string input accept only file or filename
    # FIXME check for typos in the file name
    elif os.path.exists(input) or (len(input) < 80 and input.count("\n") == 0):
        with open(input, "r") as infile:
            inptstr = infile.read()
    else:
        inptstr = input

    return inptstr


# End of file
