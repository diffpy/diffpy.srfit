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
"""Structure wrapper class for structures compatible with SrReal.
"""

__all__ = ["BaseStructure"]

from .basestructure import BaseStructure
from .bvsrestraint import BVSRestraint

class SrRealStructure(BaseStructure):
    """Base class for SrReal-compatible structure adapters.

    This derives from BaseStructure and provides some extended functionality
    provided by SrReal.

    Attributes:
    stru    --  The adapted object

    """

    def restrainBVS(self, prefactor = 1, scaled = False):
        """Restrain the bond-valence sum to zero.

        This adds a penalty to the cost function equal to
        prefactor * bvrmsdiff
        where bvrmsdiff is the rmsdifference between the calculated and
        expected bond valence sums for the structure. If scaled is true, this
        is also scaled by the current chi^2 value so the restraint is roughly
        equally weighted in the fit.

        prefactor   --  A multiplicative prefactor for the restraint 
                        (default 1).
        scaled      --  A flag indicating if the restraint is scaled
                        (multiplied) by the unrestrained point-average chi^2
                        (chi^2/numpoints) (default False).

        Returns the BVSRestraint object for use with the 'unrestrain' method.

        """

        # Create the Restraint object
        res = BVSRestraint(self.stru, prefactor, scaled)
        # Add it to the _restraints set
        self._restraints.add(res)
        # Our configuration changed. Notify observers.
        self._updateConfiguration()
        # Return the Restraint object
        return res

