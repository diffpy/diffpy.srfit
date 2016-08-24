#!/usr/bin/env python
##############################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2009 The Trustees of Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Chris Farrow
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
##############################################################################

"""Structure wrapper class for structures compatible with SrReal.
"""

__all__ = ["SrRealParSet"]


from diffpy.srfit.structure.basestructureparset import BaseStructureParSet
from diffpy.srfit.structure.bvsrestraint import BVSRestraint

class SrRealParSet(BaseStructureParSet):
    """Base class for SrReal-compatible structure adapters.

    This derives from BaseStructureParSet and provides some extended
    functionality provided by SrReal.

    Attributes:
    stru            --  The adapted object
    _usesymmetry    --  A flag indicating if SrReal calculators that operate on
                        this object should use symmetry. By default this is
                        True.

    """

    def __init__(self, *args, **kw):
        BaseStructureParSet.__init__(self, *args, **kw)
        self._usesymmetry = True
        self.stru = None
        return

    def restrainBVS(self, sig = 1, scaled = False):
        """Restrain the bond-valence sum to zero.

        This adds a penalty to the cost function equal to
        bvmsdiff / sig**2
        where bvmsdiff is the mean-squared difference between the calculated
        and expected bond valence sums for the structure. If scaled is True,
        this is also scaled by the current point-averaged chi^2 value so the
        restraint is roughly equally weighted in the fit.

        sig     --  The uncertainty on the BVS (default 1).
        scaled  --  A flag indicating if the restraint is scaled
                    (multiplied) by the unrestrained point-average chi^2
                    (chi^2/numpoints) (default False).

        Returns the BVSRestraint object for use with the 'unrestrain' method.

        """

        # Create the Restraint object
        res = BVSRestraint(self, sig, scaled)
        # Add it to the _restraints set
        self._restraints.add(res)
        # Our configuration changed. Notify observers.
        self._updateConfiguration()
        # Return the Restraint object
        return res

    def useSymmetry(self, use = True):
        """Set this structure to use symmetry.

        This determines how the structure is treated by SrReal calculators.

        """
        self._usesymmetry = bool(use)
        return

    def usingSymmetry(self):
        """Check if symmetry is being used."""
        return self._usesymmetry

    def _getSrRealStructure(self):
        """Get the structure object for use with SrReal calculators.

        If this is periodic, then return the structure, otherwise, pass it
        inside of a nosymmetry wrapper.

        """
        from diffpy.srreal.structureadapter import nosymmetry
        if self._usesymmetry:
            return self.stru
        return nosymmetry(self.stru)
