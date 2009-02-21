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

from diffpy.srfit.equation.visitors.Visitor import Visitor

from diffpy.srfit.equation.literals import Clicker

class Evaluator(Visitor):
    """Evaluator visitor."""

    def __init__(self):
        """Initialize."""
        self.value = 0
        self.clicker = Clicker()
        return

    def onArgument(self, arg):
        """Process an Argument node."""
        self.value = arg.value
        return

    def onOperator(self, op):
        """Process an Operator node."""

        # We have to re-evaluate if the operator's clicker is greater than the
        # evaluator's. 
        if op.clicker > self.clicker:

            # Visit each literal in the operator's arg list and update its
            # value.
            for literal in op.args:
                literal.identify(self)

            vals = [arg.value for arg in op.args]
            op.value = op.operation(*vals)

        self.value = op.value
        return


# version
__id__ = "$Id$"

#
# End of file
