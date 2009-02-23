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
"""Validator visitor for validating a tree of Literals.

The Validator walks an equation tree composed of Literals and checks the
validity of each equation as much as possible without evaluating it. It
collects errors in a list.
"""

from diffpy.srfit.equation.visitors.Visitor import Visitor

from diffpy.srfit.equation.literals import Clicker

class Validator(Visitor):
    """Validator error for checking validity of an equation tree.

    Attributes:
    errors  --  List of strings describing errors.
    """

    def __init__(self):
        """Initialize."""
        self.errors = []
        self._nin = 0
        return

    def reset(self):
        """Reset the validator.

        This clears the errors list.
        """
        self.errors = []
        return

    def onArgument(self, arg):
        """Process an Argument node.

        No assumption is made about the argument type.
        """
        self._nin = 1
        return

    def onOperator(self, op):
        """Process an Operator node."""
        # Check name
        if op.name is None:
            m = "'%s' does not have a name"%op
            self.errors.append(m)
        # Check symbol
        if op.symbol is None:
            m = "'%s' does not have a symbol"%op
            self.errors.append(m)
        # Check operation without evaluating it
        if op.operation is None:
            m = "'%s' does not define and operation"%op
            self.errors.append(m)

        localnin = 0
        for literal in op.args:
            literal.identify(self)
            # Count the nout of the arguments.
            localnin += self._nin


        # Check the input/output balance
        if localnin != op.nin:
            m = "'%s' requires %i inputs but receives %i"%(op, op.nin, localnin)
            self.errors.append(m)

        self._nin = op.nout
        return


# version
__id__ = "$Id$"

#
# End of file
