#!/usr/bin/env python
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
"""Restraints class. 

Restraints are used by RecipeOrganizers to organize restraint equations.
Restraints store an Equation, bounds on its value, and the form of the penalty
function for breaking a restraint. This penalty is added to the residual
equation calculated by a FitRecipe.

"""
__all__ = ["Restraint"]

from numpy import inf

class Restraint(object):
    """Restraint class.

    Attributes
    eq      --  An equation whose evaluation is compared against the restraint
                bounds.
    lb      --  The lower bound on the restraint evaluation (default -inf).
    ub      --  The lower bound on the restraint evaluation (default inf).
    prefactor   --  A multiplicative prefactor for the restraint (default 1).
    power   --  The power of the penalty (default 2).
    scaled  --  A flag indicating if the restraint is scaled (multiplied) by
                the unrestrained point-average chi^2 (chi^2/numpoints) 
                (default False).

    The penalty is calculated as 
    prefactor * max(0, lb - val, val - ub) ** power
    and val is the value of the calculated equation.  This is multipled by the
    average chi^2 if scaled is True.

    """

    def __init__(self):
        """Initialization. """
        self.eq = None
        self.lb = -inf
        self.up = inf
        self.prefactor = 1
        self.power = 2
        self.scaled = False
        return

    def restrain(self, eq, lb = -inf, ub = inf, prefactor = 1, power = 2,
            scaled = False):
        """Restrain an equation to specified bounds.
        
        eq      --  An equation whose evaluation is compared against the
                    restraint bounds.
        lb      --  The lower bound on the restraint evaluation (default -inf).
        ub      --  The lower bound on the restraint evaluation (default inf).
        prefactor   --  A multiplicative prefactor for the restraint (default 1).
        power   --  The power of the penalty (default 2).
        scaled  --  A flag indicating if the restraint is scaled (multiplied)
                    by the unrestrained point-average chi^2 (chi^2/numpoints)
                    (default False).

        """
        self.eq = eq
        self.lb = lb
        self.ub = ub
        self.prefactor = prefactor
        self.power = power
        self.scaled = scaled
        return

    def penalty(self, w = 1.0):
        """Calculate the penalty of the restraint.

        w   --  The point-average chi^2 which is optionally used to scale the
                penalty (default 1.0).
        
        """
        val = self.eq()
        penalty = self.prefactor *\
                max(0, self.lb - val, val - self.ub) ** self.power

        if self.scaled:
            penalty *= w

        return penalty

# End class Restraint

# version
__id__ = "$Id$"

#
# End of file
