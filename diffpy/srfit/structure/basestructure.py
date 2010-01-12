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
from diffpy.srfit.util import public

@public
class BaseStructure(ParameterSet):
    """Base class for structure adapters.

    BaseStructure derives from ParameterSet and provides methods that help
    interface the ParameterSet with the space group constraint methods in the
    sgconstraints module and to ProfileGenerators.

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

    def getScatterers(self):
        """Get a list of ParameterSets that represents the scatterers.

        The site positions must be accessible from the list entries via the
        names "x", "y", and "z". The ADPs must be accessible as well, but the
        name and nature of the ADPs (U-factors, B-factors, isotropic,
        anisotropic) depends on the adapted structure.

        """
        raise NotImplementedError("The must be overloaded")

    def getSpaceGroup(self):
        """Get the HM space group symbol for the structure."""
        raise NotImplementedError("The must be overloaded")


