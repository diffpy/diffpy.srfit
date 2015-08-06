#!/usr/bin/env python

########################################################################
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
########################################################################
"""Building blocks for defining a lazy evaluation network.

Literals are the building blocks of the evaluation network. An Argument holds
the name and value of an equation variable. Operators are used to compose
other Literals to produce a new value.

Literal networks can be evaluated or have other actions performed on them by
Visitors (in diffpy.srfit.equation.visitors). The Literal-Visitor relationship
is that described by the Visitor pattern
(http://en.wikipedia.org/wiki/Visitor_pattern).

"""
from __future__ import print_function
import six

__all__ = ["Argument", "Operator", "AdditionOperator", "SubtractionOperator",
        "MultiplicationOperator", "DivisionOperator", "ExponentiationOperator",
        "RemainderOperator", "NegationOperator", "ConvolutionOperator",
        "SumOperator", "UFuncOperator", "ListOperator", "SetOperator",
        "ArrayOperator", "PolyvalOperator"]


# Import the operators

from diffpy.srfit.equation.literals.argument import Argument
from diffpy.srfit.equation.literals.operators import Operator
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
from diffpy.srfit.equation.literals.operators import ListOperator
from diffpy.srfit.equation.literals.operators import SetOperator
from diffpy.srfit.equation.literals.operators import ArrayOperator
from diffpy.srfit.equation.literals.operators import PolyvalOperator

# Try some optimizations on these classes
try:
    import psyco
    psyco.profile()
except ImportError:
    pass
