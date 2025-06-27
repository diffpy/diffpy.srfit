#!/usr/bin/env python
##############################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2008 The Trustees of Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Chris Farrow
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
##############################################################################
"""Name utilities."""

__all__ = ["isIdentifier", "validateName"]

import re

reident = re.compile(r"^[a-zA-Z_]\w*$")


def isIdentifier(s):
    """Check to see if a python string is a valid identifier.

    From http://code.activestate.com/recipes/413487/
    """
    if reident.match(s) is None:
        return False
    return True


def validateName(name):
    """Validate that a name is a valid identifier.

    Raises ValueError if the name is invalid.
    """
    # Check that the name is valid
    if not isIdentifier(name):
        raise ValueError(f"Name {name} is not a valid identifier")
    return


# End of file
