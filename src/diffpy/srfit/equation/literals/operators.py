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

"""Operator classes.

Operators are combined with other Literals to create an equation. Operators are
non-leaf nodes on a Literal tree. These trees can be evaluated by the Evaluator
visitor, or otherwise inspected.

The Operator class contains all the information necessary to be identified and
evaluated by a Visitor. Thus, a single onOperator method exists in the Visitor
base class. Other Operators can be derived from Operator (see AdditionOperator),
but they all identify themselves with the Visitor.onOperator method.
"""

__all__ = ["Operator", "AdditionOperator", "SubtractionOperator",
           "MultiplicationOperator", "DivisionOperator", "ExponentiationOperator",
           "RemainderOperator", "NegationOperator", "ConvolutionOperator",
           "SumOperator", "UFuncOperator", "ListOperator", "SetOperator",
           "ArrayOperator", "PolyvalOperator"]

import numpy

from diffpy.srfit.equation.literals.abcs import OperatorABC
from diffpy.srfit.equation.literals.literal import Literal


class Operator(Literal, OperatorABC):
    """Class for holding a general operator.

    This holds a general operator and records its function, arguments, name and
    symbol.  The visitors should be able to process any Operator with this
    information alone.

    Attributes
    args    --  List of Literal arguments, set with 'addLiteral'
    name    --  A name for this operator. e.g. "add" or "sin"
    nin     --  Number of inputs (<1 means this is variable)
    nout    --  Number of outputs
    operation   --  Function that performs the operation. e.g. numpy.add.
    symbol  --  The symbolic representation. e.g. "+" or "sin".
    _value  --  The value of the Operator.
    value   --  Property for 'getValue'.

    """

    args = None
    nin = None
    nout = None
    operation = None
    symbol = None
    _value = None

    def __init__(self, name = None, symbol = None, operation = None, nin = 2,
            nout = 1):
        """Initialization."""
        Literal.__init__(self, name)
        self.symbol = symbol
        self.nin = nin
        self.nout = nout
        self.args = []
        self.operation = operation
        return

    def identify(self, visitor):
        """Identify self to a visitor."""
        return visitor.onOperator(self)

    def addLiteral(self, literal):
        """Add a literal to this operator.

        Note that order of operation matters. The first literal added is the
        leftmost argument. The last is the rightmost.

        Raises ValueError if the literal causes a self-reference.

        """
        # Make sure we don't have self-reference
        self._loopCheck(literal)
        self.args.append(literal)
        literal.addObserver(self._flush)
        self._flush(other=(self,))
        return

    def getValue(self):
        """Get or evaluate the value of the operator."""
        if self._value is None:
            vals = [l.value for l in self.args]
            self._value = self.operation(*vals)
        return self._value

    value = property(lambda self: self.getValue())

    def _loopCheck(self, literal):
        """Check if a literal causes self-reference."""
        if literal is self:
            raise ValueError("'%s' causes self-reference"%literal)

        # Check to see if I am a dependency of the literal.
        if hasattr(literal, "args"):
            for l in literal.args:
                self._loopCheck(l)
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


def _conv(v1, v2):
    # Get the full convolution
    c = numpy.convolve(v1, v2, mode="full")
    # Find the centroid of the first signal
    s1 = sum(v1)
    x1 = numpy.arange(len(v1), dtype=float)
    c1idx = numpy.sum(v1 * x1)/s1
    # Find the centroid of the convolution
    xc = numpy.arange(len(c), dtype=float)
    ccidx = numpy.sum(c * xc)/sum(c)
    # Interpolate the convolution such that the centroids line up. This
    # uses linear interpolation.
    shift = ccidx - c1idx
    x1 += shift
    c = numpy.interp(x1, xc, c)

    # Normalize
    sc = sum(c)
    if sc > 0:
        c *= s1/sc

    return c

class ConvolutionOperator(Operator):
    """Convolve two signals.

    This convolves two signals such that centroid of the first array is not
    altered by the convolution. Furthermore, the integrated amplitude of the
    convolution is scaled to be that of the first signal. This is mean to act
    as a convolution of a signal by a probability distribution.

    Note that this is only possible when the signals are computed over the same
    range.

    """

    def __init__(self):
        """Initialization."""
        Operator.__init__(self)
        self.name = "convolve"
        self.symbol = "convolve"
        self.operation = _conv
        return

class SumOperator(Operator):
    """numpy.sum operator."""

    def __init__(self):
        """Initialization."""
        Operator.__init__(self)
        self.name = "sum"
        self.symbol = "sum"
        self.nin = 1
        self.nout = 1
        self.operation = numpy.sum
        return

class UFuncOperator(Operator):
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

def _makeList(*args):
    return args

class ListOperator(Operator):
    """Operator that will take parameters and turn them into a list."""

    def __init__(self):
        """Initialization."""
        Operator.__init__(self)
        self.name = "list"
        self.symbol = "list"
        self.nin = -1
        self.operation = _makeList
        return

def _makeSet(*args):
    return set(args)

class SetOperator(Operator):
    """Operator that will take parameters and turn them into a set."""

    def __init__(self):
        """Initialization."""
        Operator.__init__(self)
        self.name = "set"
        self.symbol = "set"
        self.nin = -1
        self.operation = _makeSet
        return

def _makeArray(*args):
    return numpy.array(args)

class ArrayOperator(Operator):
    """Operator that will take parameters and turn them into an array."""

    def __init__(self):
        """Initialization."""
        Operator.__init__(self)
        self.name = "array"
        self.symbol = "array"
        self.nin = -1
        self.operation = _makeArray
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

# End of file
