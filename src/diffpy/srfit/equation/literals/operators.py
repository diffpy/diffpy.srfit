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

    FIXME

    args    --  List of Literal arguments, set with 'addLiteral'
    name    --  A name for this operator. e.g. "add" or "sin"
    nin     --  Number of inputs (-1 means this is variable)
    nout    --  Number of outputs
    operation   --  Function that performs the operation. e.g. numpy.add.
    symbol  --  The symbolic representation. e.g. "+" or "sin".
    _value  --  The value of the Operator.
    value   --  Property for 'getValue'.
    """

    args = None
    _value = None


    def __init__(self, name=None):
        """Initialization."""
        Literal.__init__(self, name)
        self.args = []
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


class BinaryOperator(Operator):

    nin = 2
    nout = 1
    pass


class CustomOperator(Operator):

    # non-abstract user-defined operator
    nin = None
    nout = None
    operation = None
    symbol = None


def makeOperator(name, symbol, operation, nin, nout):
    op = CustomOperator(name=name)
    op.symbol = symbol
    op.operation = operation
    op.nin = nin
    op.nout = nout
    return op

# Some specified operators


class AdditionOperator(BinaryOperator):
    """Addition operator."""

    name = "add"
    symbol = "+"
    operation = staticmethod(numpy.add)
    pass


class SubtractionOperator(BinaryOperator):
    """Subtraction operator."""

    name = "subtract"
    symbol = "-"
    operation = staticmethod(numpy.subtract)
    pass


class MultiplicationOperator(BinaryOperator):
    """Multiplication operator."""

    name = "multiply"
    symbol = "*"
    operation = staticmethod(numpy.multiply)
    pass


class DivisionOperator(BinaryOperator):
    """Division operator."""

    name = "divide"
    symbol = "/"
    operation = staticmethod(numpy.divide)
    pass


class ExponentiationOperator(BinaryOperator):
    """Exponentiation operator."""

    name = "power"
    symbol = "**"
    operation = staticmethod(numpy.power)
    pass


class RemainderOperator(BinaryOperator):
    """Remainder operator."""

    name = "mod"
    symbol = "%"
    operation = staticmethod(numpy.mod)
    pass


class NegationOperator(Operator):
    """Negation operator."""

    name = "negative"
    symbol = "-"
    nin = 1
    nout = 1
    operation = staticmethod(numpy.negative)
    pass


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


class ConvolutionOperator(BinaryOperator):
    """Convolve two signals.

    This convolves two signals such that centroid of the first array is not
    altered by the convolution. Furthermore, the integrated amplitude of the
    convolution is scaled to be that of the first signal. This is mean to act
    as a convolution of a signal by a probability distribution.

    Note that this is only possible when the signals are computed over the same
    range.
    """

    name = "convolve"
    symbol = "convolve"
    operation = staticmethod(_conv)
    pass


class SumOperator(Operator):
    """numpy.sum operator."""

    name = "sum"
    symbol = "sum"
    nin = 1
    nout = 1
    operation = staticmethod(numpy.sum)


class UFuncOperator(Operator):
    """A operator wrapper around a numpy ufunc.

    The name and symbol attributes are set equal to the ufunc.__name__
    attribute. nin and nout are also taken from the ufunc.
    """

    symbol = None
    nin = None
    nout = None
    operation = None

    def __init__(self, op):
        """Initialization.

        Arguments

        op  --  A numpy ufunc
        """
        Operator.__init__(self, name=op.__name__)
        self.symbol = op.__name__
        self.nin = op.nin
        self.nout = op.nout
        self.operation = op
        return


class ListOperator(Operator):
    """Operator that will take parameters and turn them into a list."""

    name = "list"
    symbol = "list"
    nin = -1
    nout = 1

    @staticmethod
    def operation(*args):
        """Convert arguments into a list."""
        return args


class SetOperator(Operator):
    """Operator that will take parameters and turn them into a set."""

    name = "set"
    symbol = "set"
    nin = -1
    nout = 1

    @staticmethod
    def operation(*args):
        """Convert arguments into a set."""
        return set(args)


class ArrayOperator(Operator):
    """Operator that will take parameters and turn them into an array."""

    name = "array"
    symbol = "array"
    nin = -1
    nout = 1

    @staticmethod
    def operation(*args):
        return numpy.array(args)


class PolyvalOperator(BinaryOperator):
    """Operator for numpy polyval."""

    name = "polyval"
    symbol = "polyval"
    operation = staticmethod(numpy.polyval)
    pass

# End of file
