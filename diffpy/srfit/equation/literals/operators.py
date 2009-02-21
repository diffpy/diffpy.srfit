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
"""Operator classes. 

Operators are the joining nodes of the equation tree...
"""

from diffpy.srfit.equation.literals.Literal import Literal

class Operator(Literal):
    """Class for holding a general operator.

    This holds a general operator and records its arguments and name. This
    should be the only Operator that visitors know how to process.
    """

    def __init__(self):
        """Initialization."""
        Literal.__init__(self)
        self.name = None
        self.symbol = None
        self.numargs = 0
        self.args = []
        self.operation = None
        return

    def identify(self, visitor):
        """Identify self to a visitor."""
        visitor.onOperator(self)
        return

    def addLiteral(self, literal):
        """Add a literal to this operator.

        Raises AttributeError if the number of literals would exceed the number
        of arguments allowed by the operator.
        """
        if len(self.args) >= self.numargs:
            raise AttributeError("Cannot accept another literal.")
        self.args.append(literal)
        self.clicker.addSubject(literal.clicker)
        literal.clicker.update()
        return

# Some specified operators

import numpy

class AdditionOperator(Operator):
    """Addition operator."""

    def __init__(self):
        """Initialization."""
        Operator.__init__(self)
        self.name = "add"
        self.symbol = "+"
        self.numargs = 2
        self.operation = numpy.add
        return

class SubtractionOperator(Operator):
    """Subtraction operator."""

    def __init__(self):
        """Initialization."""
        Operator.__init__(self)
        self.name = "subtract"
        self.symbol = "-"
        self.numargs = 2
        self.operation = numpy.subtract
        return

class MultiplicationOperator(Operator):
    """Multiplication operator."""

    def __init__(self):
        """Initialization."""
        Operator.__init__(self)
        self.name = "multiply"
        self.symbol = "*"
        self.numargs = 2
        self.operation = numpy.multiply
        return

class DivisionOperator(Operator):
    """Division operator."""

    def __init__(self):
        """Initialization."""
        Operator.__init__(self)
        self.name = "divide"
        self.symbol = "/"
        self.numargs = 2
        self.operation = numpy.divide
        return

class ExponentiationOperator(Operator):
    """Exponentiation operator."""

    def __init__(self):
        """Initialization."""
        Operator.__init__(self)
        self.name = "power"
        self.symbol = "**"
        self.numargs = 2
        self.operation = numpy.power
        return


class RemainderOperator(Operator):
    """Remainder operator."""

    def __init__(self):
        """Initialization."""
        Operator.__init__(self)
        self.name = "mod"
        self.symbol = "%"
        self.numargs = 2
        self.operation = numpy.mod
        return


# version
__id__ = "$Id$"

#
# End of file
