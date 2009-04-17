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

Restraints are used by a FitModel to organize restraint equations. Restraints
store an Equation, bounds on its value, and the form of the penalty function
for breaking a restraint.
"""

from numpy import inf
from diffpy.srfit.equation import Clicker 

class Restraint(object):
    """Constraint class.

    Attributes
    eq      --  An equation whose evaluation is used to set check against the
                restraint.
    lb      --  The lower bound on the restraint evaluation (default -inf).
    ub      --  The lower bound on the restraint evaluation (default inf).
    prefactor   --  A multiplicative prefactor for the restraint (default 1).
    power   --  The power of the penalty (default 2).
    _value  --  The last evaluated value of the restraint penalty
    _clicker    --  Used to compare with the Equation clicker in order to update
                    _value.

    The penalty is calculated as 
    prefactor * max(0, lb - val, val - ub) ** power
    and val is the value of the calculated equation.
    """

    def __init__(self):
        """Initialization. """
        self._clicker = Clicker()
        self.eq = None
        self.lb = -inf
        self.up = inf
        self.prefactor = 1
        self.power = 2
        self._value = 0
        return

    def restrain(self, eq, lb = -inf, ub = inf, prefactor = 1, power = 2):
        """Constrain a Parameter according to an Equation."""
        self.eq = eq
        self.lb = lb
        self.ub = ub
        self.prefactor = prefactor
        self.power = power
        return

    def _penalty(self):
        """Calculate the penalty from scratch."""
        val = self.eq()
        return self.prefactor *\
                max(0, self.lb - val, val - self.ub) ** self.power

    def penalty(self):
        """Calculate the penalty of the restraint."""
        # This only recalculates the penalty if the equation will return a new
        # value. Based on some simple time trials, the worst slowdown this
        # results in is a 21% increase in evaluation time. If the Equation is
        # not going to change that often, then there is considerable speed-up.
        # For a restraint on a single variable.
        ## chance of change     evaluation time vs without clock
        #   0%                  19%
        #   10%                 31%
        #   20%                 43%
        #   30%                 51%
        #   40%                 64%
        #   50%                 73%
        #   60%                 82%
        #   70%                 93%
        #   80%                 101%
        #   90%                 111%
        #   100%                121%
        # The Levenberg-Marquardt algorithm using numerical derivatives will
        # usually change only two values at a time. This means that for a fit
        # problem with 3 variables or more, the clicker check will speed things
        # up.  The Simplex optimization scheme will usually change all the
        # variables at once, which means that this check is literally a waste
        # of time.
        if self._clicker < self.eq.root.clicker:
            self._value = self._penalty()
            self._clicker.click()
        return self._value

# version
__id__ = "$Id$"

#
# End of file
