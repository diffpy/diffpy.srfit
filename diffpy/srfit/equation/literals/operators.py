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

These classes are combined with Arguments to create an equation. Each Operator
can be associated with Literals to create an equation that can be evaluated by
the Evaluator visitor, or otherwise inspected. Operators contain data 'value'
and 'clicker' attributes for facilitating rapid re-evaluation of an equation
tree.

The Operator class contains all the information necessary to be identified and
evaluated by a Visitor. Thus, a single onOperator method exists in the Visitor
base class. Other Operators can be derived from Operator (see AdditionOperator),
but they all identify themselves with the Visitor.onOperator method.
"""

from diffpy.srfit.equation.literals.Literal import Literal
import numpy

class Operator(Literal):
    """Class for holding a general operator.

    This holds a general operator and records its function, arguments, name and
    symbol.  The visitors should be able to process any Operator with this
    information alone.

    Attributes
    name    --  A name for this operator. e.g. "add" or "sin"
    symbol  --  The symbolic representation. e.g. "+" or "sin"
    nin     --  Number of inputs
    nout    --  Number of outputs
    args    --  List of Literal arguments, set with addLiteral
    operation   --  Function that performs the operation. e.g. numpy.add or
                    numpy.sin
    clicker --  A Clicker instance for recording change in the value
    value   --  The evaluated value of this Operator.
    """

    def __init__(self):
        """Initialization."""
        Literal.__init__(self)
        self.symbol = None
        self.nin = 2
        self.nout = 1
        self.args = []
        self.operation = None
        return

    def identify(self, visitor):
        """Identify self to a visitor."""
        visitor.onOperator(self)
        return

    def addLiteral(self, literal):
        """Add a literal to this operator.

        Note that order of operation matters. The first literal added is the
        leftmost argument. The last is the rightmost.
        """
        self.args.append(literal)
        self.clicker.addSubject(literal.clicker)
        return

# Some specified operators


class AdditionOperator(Operator):
    """Addition operator."""

    def __init__(self):
        """Initialization."""
        Operator.__init__(self)
        self.name = "add"
        self.symbol = "+"
        self.operation = numpy.add
        return

class SubtractionOperator(Operator):
    """Subtraction operator."""

    def __init__(self):
        """Initialization."""
        Operator.__init__(self)
        self.name = "subtract"
        self.symbol = "-"
        self.operation = numpy.subtract
        return

class MultiplicationOperator(Operator):
    """Multiplication operator."""

    def __init__(self):
        """Initialization."""
        Operator.__init__(self)
        self.name = "multiply"
        self.symbol = "*"
        self.operation = numpy.multiply
        return

class DivisionOperator(Operator):
    """Division operator."""

    def __init__(self):
        """Initialization."""
        Operator.__init__(self)
        self.name = "divide"
        self.symbol = "/"
        self.operation = numpy.divide
        return

class ExponentiationOperator(Operator):
    """Exponentiation operator."""

    def __init__(self):
        """Initialization."""
        Operator.__init__(self)
        self.name = "power"
        self.symbol = "**"
        self.operation = numpy.power
        return


class RemainderOperator(Operator):
    """Remainder operator."""

    def __init__(self):
        """Initialization."""
        Operator.__init__(self)
        self.name = "mod"
        self.symbol = "%"
        self.operation = numpy.mod
        return

class UfuncOperator(Operator):
    """A operator wrapper around a numpy ufunc.

    The name and symbol attributes are set equal to the ufunc.__name__
    attribute. nin and nout are also taken from the ufunc.
    
    """

    def __init__(self, op):
        """Initialization.

        Arguments
        op  --  A numpy ufunc
        """
        Operator.__init__(self)
        self.name = op.__name__
        self.symbol = op.__name__
        self.nin = op.nin
        self.nout = op.nout
        self.operation = op
        return

# version
__id__ = "$Id$"

#
# End of file
