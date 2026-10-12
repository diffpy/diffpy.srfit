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
"""Validator visitor for validating a tree of Literals.

The Validator walks an equation tree composed of Literals and checks the
validity of each equation as much as possible without evaluating it. It
collects errors in a list.

The Validator checks that the input count of each Operator is equal to
the output count of its arguments. It also checks that each object has
the proper interface.
"""

__all__ = ["Validator"]

from diffpy.srfit.equation.literals.abcs import ArgumentABC, OperatorABC
from diffpy.srfit.equation.visitors.visitor import Visitor

msg = "'%s' does not have the interface required by '%s'"


class Validator(Visitor):
    """Validator error for checking validity of an equation tree.

    Attributes
    ----------
    errors : list of str
        The descriptions of the errors found.
    _nin : int
        The number of outputs of the last node visited, used to count
        the inputs of its parent operator.
    """

    def __init__(self):
        self.reset()
        return

    def reset(self):
        """Clear the errors and the input count."""
        self.errors = []
        self._nin = 0
        return

    def on_argument(self, arg):
        """Check an Argument node.

        The Argument must be an instance of ArgumentABC from
        diffpy.srfit.equation.literals.abcs.

        Parameters
        ----------
        arg : Argument
            The Argument to check.

        Returns
        -------
        list of str
            The errors found so far.
        """
        if not isinstance(arg, ArgumentABC):
            m = msg % (arg, ArgumentABC.__name__)
            self.errors.append(m)
        self._nin = 1
        return self.errors

    def on_operator(self, op):
        """Check an Operator node and its arguments.

        The Operator must be an instance of OperatorABC from
        diffpy.srfit.equation.literals.abcs, be single-valued, have a
        name, symbol and operation, and receive as many inputs as it
        takes.

        Parameters
        ----------
        op : Operator
            The Operator to check.

        Returns
        -------
        list of str
            The errors found so far.
        """
        if not isinstance(op, OperatorABC):
            m = msg % (op, OperatorABC.__name__)
            self.errors.append(m)

        # Can only process single-valued functions
        if op.nout != 1:
            m = "'%s' is not single-valued (nout != 1)" % op
            self.errors.append(m)
        # Check name
        if op.name is None:
            m = "'%s' does not have a name" % op
            self.errors.append(m)
        # Check symbol
        if op.symbol is None:
            m = "'%s' does not have a symbol" % op
            self.errors.append(m)
        # Check operation without evaluating it
        if op.operation is None:
            m = "'%s' does not define and operation" % op
            self.errors.append(m)

        localnin = 0
        for literal in op.args:
            literal.identify(self)
            # Count the nout of the arguments.
            localnin += self._nin

        # Check the input/output balance
        if op.nin >= 0 and localnin != op.nin:
            m = "'%s' requires %i inputs but receives %i" % (
                op,
                op.nin,
                localnin,
            )
            self.errors.append(m)

        self._nin = op.nout
        return self.errors


# End of file
