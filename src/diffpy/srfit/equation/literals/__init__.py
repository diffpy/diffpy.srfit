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
"""Building blocks for defining a lazy evaluation network.

Literals are the building blocks of the evaluation network. An Argument
holds the name and value of an equation variable. Operators are used to
compose other Literals to produce a new value.

Literal networks can be evaluated or have other actions performed on
them by Visitors (in diffpy.srfit.equation.visitors). The Literal-
Visitor relationship is that described by the Visitor pattern (

http://en.wikipedia.org/wiki/Visitor_pattern).
"""

__all__ = [
    "Argument",
    "Operator",
    "BinaryOperator",
    "CustomOperator",
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
    "makeOperator",
]


# Import the operators

from diffpy.srfit.equation.literals.argument import Argument
from diffpy.srfit.equation.literals.operators import (
    AdditionOperator,
    ArrayOperator,
    BinaryOperator,
    ConvolutionOperator,
    CustomOperator,
    DivisionOperator,
    ExponentiationOperator,
    MultiplicationOperator,
    NegationOperator,
    Operator,
    PolyvalOperator,
    RemainderOperator,
    SubtractionOperator,
    SumOperator,
    UFuncOperator,
    makeOperator,
)

# End of file
