#!/usr/bin/env python
########################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2009 Trustees of the Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Chris Farrow
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
########################################################################
"""Base structure wrapper class.

The BaseStructure is a ParameterSet with functionality required by all
structure adaptors.

"""

from diffpy.srfit.fitbase.parameterset import ParameterSet

class BaseStructure(ParameterSet):
    """Base class for structure adapters.

    Attributes:
    stru    --  The adapted object

    """

    @classmethod
    def canAdapt(self, stru):
        """Return whether the structure can be adapted by this class."""
        return False

    def getLattice(self):
        """Get a ParameterSet containing the lattice Parameters.

        The returned ParameterSet may contain other Parameters than the lattice
        Parameters. It is assumed that the lattice parameters are named "a",
        "b", "c", "alpha", "beta", "gamma".
        
        """
        raise NotImplementedError("The must be overloaded")

    def getSites(self):
        """Get a list of ParameterSets that represents the sites.

        The site positions must be accessible from the list entries via the
        names "x", "y", and "z".

        """
        raise NotImplementedError("The must be overloaded")

