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

Operators are combined with other Literals to create an equation.
Operators are non-leaf nodes on a Literal tree. These trees can be
evaluated by the Evaluator visitor, or otherwise inspected.

The Operator class contains all the information necessary to be
identified and evaluated by a Visitor. Thus, a single onOperator method
exists in the Visitor base class. Other Operators can be derived from
Operator (see AdditionOperator), but they all identify themselves with
the Visitor.onOperator method.
"""

__all__ = [
    "Operator",
    "AdditionOperator",
    "SubtractionOperator",
    "MultiplicationOperator",
    "DivisionOperator",
    "ExponentiationOperator",
    "RemainderOperator",
    "NegationOperator",
    "ConvolutionOperator",
    "SumOperator",
    "UFuncOperator",
    "ArrayOperator",
    "PolyvalOperator",
]

import numpy

from diffpy.srfit.equation.literals.abcs import OperatorABC
from diffpy.srfit.equation.literals.literal import Literal


class Operator(Literal, OperatorABC):
    """Abstract class for specifying a general operator.

    This class provides several methods that are common to a derived
    classes for concrete operations.

    Class Attributes
    ----------------
    nin : int, abstract
        Number of input arguments for the operator.  Any number of
        arguments is allowed when -1.  This attribute must be defined
        in a derived class.
    nout : int, abstract
        Number of outputs returned by the `operation`.  This attribute
        must be defined in a derived class.
    operation : callable, abstract
        Function that performs the operation, e.g., `numpy.add`.
        This must be defined in a derived class.
    symbol : str, abstract
        The symbolic representation for the operator such as
        "+" or "sin".  This attribute must be defined in a derived
        class.

    Attributes
    ----------
    args : list
        The list of `Literal` arguments.  Read-only, use the
        `addLiteral` method to change its content.
    """

    # Private Attributes
    # ------------------
    # _value : float, numpy.ndarray or None
    #     The last value of the operator or None.

    # We must declare the abstract `args` here.
    args = None
    # default for the value
    _value = None

    def __init__(self, name=None):
        """Initialize the operator object with the specified name.

        Parameters
        ----------
        name : str, optional
            Name for the operator object.  When not specified,
            use the class attribute `name`.
        """
        Literal.__init__(self, name)
        self.args = []
        return

    def identify(self, visitor):
        """Identify self to a visitor."""
        return visitor.onOperator(self)

    def addLiteral(self, literal):
        """Add a literal to this operator.

        Note that order of operation matters. The first literal added is
        the leftmost argument. The last is the rightmost.

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
            vals = [arg.value for arg in self.args]
            self._value = self.operation(*vals)
        return self._value

    value = property(lambda self: self.getValue())

    def _loopCheck(self, literal):
        """Check if a literal causes self-reference."""
        if literal is self:
            raise ValueError("'%s' causes self-reference" % literal)

        # Check to see if I am a dependency of the literal.
        if hasattr(literal, "args"):
            for lit_arg in literal.args:
                self._loopCheck(lit_arg)
        return


class UnaryOperator(Operator):
    """Abstract class for an unary operator with one input and one result.

    This base class defines the `nin` and `nout` attributes.  The derived
    concrete operator must provide the remaining abstract attributes
    of the `Operator` class.
    """

    nin = 1
    nout = 1
    pass


class BinaryOperator(Operator):
    """Abstract class for a binary operator with two inputs and one result.

    This base class defines the `nin` and `nout` attributes.  The derived
    concrete operator must define the remaining abstract attributes
    of the `Operator` class.
    """

    nin = 2
    nout = 1
    pass


class CustomOperator(Operator):
    """Concrete class for a user-defined Operator.

    Use the `makeOperator` factory function to create an instance.
    """

    # declare all abstract attributes from the Operator base.
    nin = None
    nout = None
    operation = None
    symbol = None
    pass


def makeOperator(name, symbol, operation, nin, nout):
    """Return a new custom operator object.

    Parameters
    ----------
    name : str
        Name of the custom operator object.
    symbol : str
        The symbolic representation for the operator such as
        "+" or "sin".
    operation : callable
        Function that performs the operation, e.g., `numpy.add`.
    nin : int
        Number of input arguments for the operator.  Any number of
        arguments is allowed when -1.
    nout : in
        Number of outputs returned by the `operation`.

    Returns
    -------
    CustomOperator
        The new custom operator object.
    """
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


class NegationOperator(UnaryOperator):
    """Negation operator."""

    name = "negative"
    symbol = "-"
    operation = staticmethod(numpy.negative)
    pass


def _conv(v1, v2):
    # Get the full convolution
    c = numpy.convolve(v1, v2, mode="full")
    # Find the centroid of the first signal
    s1 = sum(v1)
    x1 = numpy.arange(len(v1), dtype=float)
    c1idx = numpy.sum(v1 * x1) / s1
    # Find the centroid of the convolution
    xc = numpy.arange(len(c), dtype=float)
    ccidx = numpy.sum(c * xc) / sum(c)
    # Interpolate the convolution such that the centroids line up. This
    # uses linear interpolation.
    shift = ccidx - c1idx
    x1 += shift
    c = numpy.interp(x1, xc, c)

    # Normalize
    sc = sum(c)
    if sc > 0:
        c *= s1 / sc

    return c


class ConvolutionOperator(BinaryOperator):
    """Convolve two signals.

    This convolves two signals such that centroid of the first array is
    not altered by the convolution. Furthermore, the integrated
    amplitude of the convolution is scaled to be that of the first
    signal. This is mean to act as a convolution of a signal by a
    probability distribution.

    Note that this is only possible when the signals are computed over
    the same range.
    """

    name = "convolve"
    symbol = "convolve"
    operation = staticmethod(_conv)
    pass


class SumOperator(UnaryOperator):
    """numpy.sum operator."""

    name = "sum"
    symbol = "sum"
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
