########################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2008 Trustees of the Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Chris Farrow
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
########################################################################

"""Utitiles used throughout SrFit.

"""

def hasinterface(obj, cls):
    """Determine whether an object has the interface defined within cls."""
    attrs = [attr for attr in dir(cls) if not attr.startswith("__")]
    for attr in attrs:
        if not hasattr(obj, attr):
            return False
    return True

def verifyinterface(obj, cls):
    """Verify that obj has the interface defined within cls.

    Raises TypeError if the interface is not satisfied.

    """
    if hasinterface(obj, cls):
        return

    m = "'%s' does not have the interface required by '%s'"%(obj, cls.__name__)
    raise TypeError(m)
    return

# End of file
