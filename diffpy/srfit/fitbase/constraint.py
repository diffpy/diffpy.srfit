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
"""Constraint class. 

Constraints are used by a FitModel (and elsewhere) to organize constraint
equations. They store a Parameter object and an Equation object that is used to
compute its value. Constraint contains methods to aid in the creation of this
association.
"""

class Constraint(object):
    """Constraint class.
    
    Attributes
    par     --  A Parameter that is the subject of the constraint.
    eq      --  An equation whose evaluation is used to set the value of the
                constraint.
    """

    def __init__(self):
        """Initialization. """
        self.par = None
        self.eq = None
        return

    def constrain(self, par, eq):
        """Constrain a Parameter according to an Equation."""
        self.par = par
        self.eq = eq
        self.update()
        return

    def update(self):
        """Update the parameter according to the equation."""
        # This will be evaluated quickly thanks to the Equation class.
        val = self.eq()
        # This will only change the Parameter if val is different from the
        # currently stored value.
        self.par.setValue(val)
        return

# version
__id__ = "$Id$"

#
# End of file
