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
"""Bond-valence sum calculator from SrReal wrapped as a Restraint.

This can be used as an addition to a cost function during a structure
refinement to keep the bond-valence sum within tolerable limits.
"""

__all__ = ["BVSRestraint"]

from diffpy.srfit.exceptions import SrFitError
from diffpy.srfit.fitbase.restraint import Restraint


class BVSRestraint(Restraint):
    """Wrapping of BVSCalculator.bvmsdiff as a Restraint.

    The restraint penalty is the root-mean-square deviation of the theoretical
    and calculated bond-valence sum of a structure.

    Attributes:
    _calc   --  The SrReal BVSCalculator instance.
    _parset --  The SrRealParSet that created this BVSRestraint.
    sig     --  The uncertainty on the BVS (default 1).
    scaled  --  A flag indicating if the restraint is scaled (multiplied)
                by the unrestrained point-average chi^2 (chi^2/numpoints)
                (default False).
    """

    def __init__(self, parset, sig=1, scaled=False):
        """Initialize the Restraint.

        parset  --  SrRealParSet that creates this BVSRestraint.
        sig     --  The uncertainty on the BVS (default 1).
        scaled  --  A flag indicating if the restraint is scaled
                    (multiplied) by the unrestrained point-average chi^2
                    (chi^2/numpoints) (bool, default False).
        """
        from diffpy.srreal.bvscalculator import BVSCalculator

        self._calc = BVSCalculator()
        self._parset = parset
        self.sig = float(sig)
        self.scaled = bool(scaled)
        return

    def penalty(self, w=1.0):
        """Calculate the penalty of the restraint.

        w   --  The point-average chi^2 which is optionally used to scale the
                penalty (float, default 1.0).
        """
        # Get the bvms from the BVSCalculator
        stru = self._parset._getSrRealStructure()
        self._calc.eval(stru)
        penalty = self._calc.bvmsdiff

        # Scale by the prefactor
        penalty /= self.sig**2

        # Optionally scale by w
        if self.scaled:
            penalty *= w

        return penalty

    def _validate(self):
        """This evaluates the calculator.

        Raises SrFitError if validation fails.
        """
        from numpy import nan

        p = self.penalty()
        if p is None or p is nan:
            raise SrFitError("Cannot evaluate penalty")
        v = self._calc.value
        if len(v) > 1 and not v.any():
            emsg = (
                "Bond valence sums are all zero.  Check atom symbols in "
                "the structure or define custom bond-valence parameters."
            )
            raise SrFitError(emsg)
        return

    # End of class BVSRestraint


# End of file
