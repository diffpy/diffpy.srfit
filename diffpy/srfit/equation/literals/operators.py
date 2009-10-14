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

Operators are combined with other Literals to create an equation. Operators are
non-leaf nodes on a Literal tree. 

"""

from .literal import Literal

import numpy

class Operator(Literal):
    """Class for holding a general operator.

    This class inherits from Literal. See the Literal documentation.

    This holds a general operator and records its function, arguments, name and
    symbol.  The visitors should be able to process any Operator with this
    information alone.

    Attributes
    args   --  List of Literal arguments, set with addLiteral.
    nin     --  Number of inputs (non-positive means this is variable)
    nout    --  Number of outputs
    operation   --  Function that performs the operation. 
                operation(*args) -> value
    symbol  --  The symbolic representation. e.g. "+" or "sin"

    """

    def __init__(self, name = None, symbol = None, operation = None, nin = 2):
        """Initialization."""
        Literal.__init__(self, name = name)
        self.symbol = symbol
        self.nin = nin
        self.nout = 1
        self.args = []
        self.operation = operation
        return

    def addLiteral(self, literal):
        """Add a literal to this operator.

        Note that order of operation matters. The first literal added is the
        leftmost argument. The last is the rightmost.

        Raises ValueError if we cannot use the argument.

        """
        if self.nin >= 0 and len(self.args) >= self.nin:
            raise ValueError("Cannot accept more arguments")
        self._validateArg(literal)

        self.args.append(literal)
        literal.addObserver(self.flush)
        self.flush(None)
        return

    def _updateValue(self):
        """Update my value if possible."""
        vals = [l.getValue() for l in self.args]
        self._value = self.operation(*vals)
        return

    def _identify(self, visitor):
        """Identify self to a visitor."""
        visitor.onOperator(self)
        return

    def __str__(self):
        if self.name:
            return "Operator(" + self.name + ")"
        return self.__repr__()

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

class NegationOperator(Operator):
    """Negation operator."""

    def __init__(self):
        """Initialization."""
        Operator.__init__(self)
        self.name = "negative"
        self.symbol = "-"
        self.nin = 1
        self.operation = numpy.negative
        return

class ConvolutionOperator(Operator):
    """Scaled version of the numpy.convolve operator.

    This calls numpy.convolve, but divides by the sum of the second argument in
    hopes of preserving the scale of the first argument.
    numpy.conolve(v1, v2, mode = "same")/sum(v2)
    It then truncates to the length of the first array.

    """

    def __init__(self):
        """Initialization."""
        Operator.__init__(self)
        self.name = "convolve"
        self.symbol = "convolve"

        def conv(v1, v2):
            """numpy.conolve(v1, v2, mode = "same")/sum(v2)"""
            c = numpy.convolve(v1, v2, mode="same")/sum(v2)
            c.resize((len(v1),))
            return c

        self.operation = conv
        return


class SumOperator(Operator):
    """numpy.sum operator."""

    def __init__(self):
        """Initialization."""
        Operator.__init__(self)
        self.name = "sum"
        self.symbol = "sum"
        self.nin = 1
        self.operation = numpy.sum
        return

class UFuncOperator(Operator):
    """A operator wrapper around a numpy ufunc.

    The name and symbol attributes are set equal to the ufunc.__name__
    attribute. nin is also taken from the ufunc.
    
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

class ListOperator(Operator):
    """Operator that will take parameters and turn them into a list."""

    def __init__(self):
        """Initialization."""
        Operator.__init__(self)
        self.name = "list"
        self.symbol = "list"
        self.nin = -1

        def makeList(*args):
            return args

        self.operation = makeList
        return

class SetOperator(Operator):
    """Operator that will take parameters and turn them into a set."""

    def __init__(self):
        """Initialization."""
        Operator.__init__(self)
        self.name = "set"
        self.symbol = "set"
        self.nin = -1

        def makeSet(*args):
            return set(args)

        self.operation = makeSet
        return

class ArrayOperator(Operator):
    """Operator that will take parameters and turn them into an array."""

    def __init__(self):
        """Initialization."""
        Operator.__init__(self)
        self.name = "array"
        self.symbol = "array"
        self.nin = -1

        def makeArray(*args):
            return numpy.array(args)

        self.operation = makeArray
        return

class PolyvalOperator(Operator):
    """Operator for numpy polyval."""

    def __init__(self):
        """Initialization."""
        Operator.__init__(self)
        self.name = "polyval"
        self.symbol = "polyval"
        self.nin = 2
        self.operation = numpy.polyval
        return

# version
__id__ = "$Id$"

#
# End of file
