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

from .Visitor import Visitor

from .. import Clicker

class Printer(Visitor):
    """Validator error for checking validity of an equation tree.

    Attributes:
    errors  --  List of strings describing errors.
    """

    def __init__(self):
        """Initialize."""
        self.output = ""
        return

    def onArgument(self, arg):
        """Process an Argument node.

        No assumption is made about the argument type.
        """
        if arg.name is None:
            self.output += str(arg.value)
        else:
            self.output += str(arg.name)
        return

    def onOperator(self, op):
        """Process an Operator node."""
        # We have to deal with infix operators
        if op.name != op.symbol and op.nin == 2:
            self._onInfix(op)
            return

        self.output += str(op.name) + "("

        for idx, literal in enumerate(op.args):
            if idx != 0: self.output += ", "
            literal.identify(self)


        self.output += ")"
        return

    def _onInfix(self, op):
        """Process infix operators."""

        self.output += "("
        op.args[0].identify(self)
        self.output += " %s "%op.symbol
        op.args[1].identify(self)
        self.output += ")"
        return

# version
__id__ = "$Id$"

#
# End of file
