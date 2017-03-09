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

"""Restraints class.

Restraints are used by RecipeOrganizers to organize restraint equations.
Restraints store an Equation, bounds on its value, and the form of the penalty
function for breaking a restraint. This penalty is added to the residual
equation calculated by a FitRecipe.
"""

__all__ = ["Restraint"]

from numpy import inf

from diffpy.srfit.fitbase.validatable import Validatable
from diffpy.srfit.exceptions import SrFitError


class Restraint(Validatable):
    """Restraint class.

    Attributes
    eq      --  An equation whose evaluation is compared against the restraint
                bounds.
    lb      --  The lower bound on the restraint evaluation (default -inf).
    ub      --  The lower bound on the restraint evaluation (default inf).
    sig     --  The uncertainty on the bounds (default 1).
    scaled  --  A flag indicating if the restraint is scaled (multiplied) by
                the unrestrained point-average chi^2 (chi^2/numpoints)
                (default False).

    The penalty is calculated as
    (max(0, lb - val, val - ub)/sig)**2
    and val is the value of the calculated equation.  This is multipled by the
    average chi^2 if scaled is True.

    """

    def __init__(self, eq, lb = -inf, ub = inf, sig = 1, scaled = False):
        """Restrain an equation to specified bounds.

        eq      --  An equation whose evaluation is compared against the
                    restraint bounds.
        lb      --  The lower bound on the restraint evaluation (float, default
                    -inf).
        ub      --  The lower bound on the restraint evaluation (float, default
                    inf).
        sig     --  The uncertainty on the bounds (default 1).
        scaled  --  A flag indicating if the restraint is scaled (multiplied)
                    by the unrestrained point-average chi^2 (chi^2/numpoints)
                    (bool, default False).

        """
        self.eq = eq
        self.lb = float(lb)
        self.ub = float(ub)
        self.sig = float(sig)
        self.scaled = bool(scaled)
        return

    def penalty(self, w = 1.0):
        """Calculate the penalty of the restraint.

        w   --  The point-average chi^2 which is optionally used to scale the
                penalty (default 1.0).

        Returns the penalty as a float

        """
        val = self.eq()
        penalty = (max(0, self.lb - val, val - self.ub) / self.sig)**2

        if self.scaled:
            penalty *= w

        return penalty

    def _validate(self):
        """Validate my state.

        This validates eq.

        Raises SrFitError if validation fails.

        """
        if self.eq is None:
            raise SrFitError("eq is None")
        from diffpy.srfit.equation.visitors import validate
        try:
            validate(self.eq)
        except ValueError, e:
            raise SrFitError(e)

        # Try to get the value of eq.
        try:
            val = self.eq()
        except TypeError, e:
            raise SrFitError("eq cannot be evaluated")
        finally:
            if val is None:
                raise SrFitError("eq evaluates to None")

        return

# End class Restraint

# End of file
