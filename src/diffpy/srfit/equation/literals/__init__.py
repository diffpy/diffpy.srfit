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

Literals are the building blocks of the evaluation network. An Argument holds
the name and value of an equation variable. Operators are used to compose
other Literals to produce a new value.

Literal networks can be evaluated or have other actions performed on them by
Visitors (in diffpy.srfit.equation.visitors). The Literal-Visitor relationship
is that described by the Visitor pattern
(http://en.wikipedia.org/wiki/Visitor_pattern).
"""

__all__ = ["Argument", "Operator", "BinaryOperator", "CustomOperator",
           "AdditionOperator", "SubtractionOperator",
           "MultiplicationOperator", "DivisionOperator", "ExponentiationOperator",
           "RemainderOperator", "NegationOperator", "ConvolutionOperator",
           "SumOperator", "UFuncOperator", "ArrayOperator", "PolyvalOperator",
           "makeOperator"]


# Import the operators

from diffpy.srfit.equation.literals.argument import Argument
from diffpy.srfit.equation.literals.operators import Operator
from diffpy.srfit.equation.literals.operators import BinaryOperator
from diffpy.srfit.equation.literals.operators import CustomOperator
from diffpy.srfit.equation.literals.operators import AdditionOperator
from diffpy.srfit.equation.literals.operators import SubtractionOperator
from diffpy.srfit.equation.literals.operators import MultiplicationOperator
from diffpy.srfit.equation.literals.operators import DivisionOperator
from diffpy.srfit.equation.literals.operators import ExponentiationOperator
from diffpy.srfit.equation.literals.operators import RemainderOperator
from diffpy.srfit.equation.literals.operators import NegationOperator
from diffpy.srfit.equation.literals.operators import ConvolutionOperator
from diffpy.srfit.equation.literals.operators import UFuncOperator
from diffpy.srfit.equation.literals.operators import SumOperator
from diffpy.srfit.equation.literals.operators import ArrayOperator
from diffpy.srfit.equation.literals.operators import PolyvalOperator
from diffpy.srfit.equation.literals.operators import makeOperator

# End of file
