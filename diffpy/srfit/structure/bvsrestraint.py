########################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2010 Trustees of the Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Chris Farrow
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
########################################################################
"""Bond-valence sum calculator from SrReal.

This can be used as an addition to a cost function during a structure
refinement.

"""

__all__ = ["BVSCalculator"]

from diffpy.srfit.fitbase.restraint import Restraint

from diffpy.srreal.srreal_ext import BVSCalculator

class BVSRestraint(Restraint):
    """Wrapping of BVSCalculator.bvrmsdiff as a Restraint.

    This creates a restraint the calculates the bond-valence sum of a
    structure.

    Attributes:
    _calc       --  The SrReal BVSCalculator instance.
    _stru       --  Structure object supported by BVSCalculator. This is the
                    structure from which the bond-valence sum is calculated.
    prefactor   --  A multiplicative prefactor for the restraint (default 1).
    scaled  --  A flag indicating if the restraint is scaled (multiplied) by
                the unrestrained point-average chi^2 (chi^2/numpoints) 
                (default False).

    """

    def __init__(self, stru):
        """Initialize the Restraint.

        stru       --  Structure object supported by BVSCalculator

        """
        self._calc = BVSCalculator()
        self._stru = stru
        self.prefactor = 1
        self.scaled = False
        return

    def restrain(self, prefactor = 1, scaled = False):
        """Set up the restraint.
        
        prefactor   --  A multiplicative prefactor for the restraint (default 1).
        scaled  --  A flag indicating if the restraint is scaled (multiplied)
                    by the unrestrained point-average chi^2 (chi^2/numpoints)
                    (default False).

        """
        self.prefactor = prefactor
        self.scaled = False
        pass

    def penalty(self, w = 1.0):
        """Calculate the penalty of the restraint.

        w   --  The point-average chi^2 which is optionally used to scale the
                penalty (default 1.0).
        
        """
        self._calc.eval(self._stru)
        penalty = self._calc.bvrmsdiff()
        penalty *= self.prefactor

        if self.scaled: penalty *= w

        return penalty



# End of file
