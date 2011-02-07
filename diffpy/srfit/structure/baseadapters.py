#!/usr/bin/env python
########################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2011 Trustees of the Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Chris Farrow
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
########################################################################
"""Base class for adapting structures to a uniform adapter interface.

"""

__all__ = ["BaseScattererAdapter", "BaseStructureAdapter",
    "BaseLatticeAdapter"]

from diffpy.srfit.adapters.adaptersmod import ContainerAdapter

class BaseScattererAdapter(ContainerAdapter):
    """Base class for scatterer adapters.

    BaseScattererAdapter provides a getXYZ method that returns a tuple of the
    adapted fractional x, y and z coordinates.

    """

    def getXYZ(self):
        """Get tuple of adapted fractional coordinates.

        The coordinates must be supplied in order: x, y, z.
        
        """
        raise NotImplementedError

    def getADPs(self):
        """Get tuple of adapted ADPs. 
        
        The ADPs must be supplied in order: 11, 22, 33, 12, 13, 23, iso.
        Whether the ADPs are U- or B-factors is dependent on the adapted
        scatterer. If any of these cannot be supplied, supply None instead. If
        none of these can be supplied, return None.

        """
        raise NotImplementedError

# End class BaseScattererAdapter

class BaseLatticeAdapter(ContainerAdapter):
    """Base class for structure adapters.

    This contains the "angunits" attribute required by the getLattice method of
    BaseStructureAdapter. (This defaults to "deg".)

    """
    angunits = "deg"

    def getLatPars(self):
        """Get a tuple of the adapted lattice parameters.

        Parameters must be supplied in order: a, b, c, alpha, beta, gamma

        """
        raise NotImplementedError

# End class BaseLatticeAdapter

class BaseStructureAdapter(ContainerAdapter):
    """Base class for structure adapters.

    BaseStructureAdapter provides methods that help interface structure
    adapters with space group functions.

    """

    def getLattice(self):
        """Get the adapted lattice.

        The returned adapter must have the interface of BaseLatticeAdapter.

        If the lattice does not exist, return None.

        """
        raise NotImplementedError

    def getScatterers(self):
        """Get a tuple of adapted scatterers.

        The scatterers must conform to the interface defined by
        BaseScattererAdapter.

        """
        raise NotImplementedError

# End class BaseStructureAdapter

__id__ = "$Id$"
