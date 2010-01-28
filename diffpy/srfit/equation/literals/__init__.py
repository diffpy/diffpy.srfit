
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
"""Building blocks for defining a lazy evaluation network.

Literals are the building blocks of the evaluation network. An Argument holds
the name and value of an equation variable. Operators are used to compose
other Literals to produce a new value.

Literal networks can be evaluated or have other actions performed on them by
Visitors (in diffpy.srfit.equation.visitors). The Literal-Visitor relationship
is that described by the Visitor pattern
(http://en.wikipedia.org/wiki/Visitor_pattern).

"""

__all__ = ["Argument", "Operator", "AdditionOperator", "SubtractionOperator",
        "MultiplicationOperator", "DivisionOperator", "ExponentiationOperator",
        "RemainderOperator", "NegationOperator", "ConvolutionOperator",
        "SumOperator", "UFuncOperator", "ListOperator", "SetOperator",
        "ArrayOperator", "PolyvalOperator"]


# package version
from diffpy.srfit.version import __version__

# Import the operators

from .argument import Argument
from .operators import Operator
from .operators import AdditionOperator
from .operators import SubtractionOperator
from .operators import MultiplicationOperator
from .operators import DivisionOperator
from .operators import ExponentiationOperator
from .operators import RemainderOperator
from .operators import NegationOperator
from .operators import ConvolutionOperator
from .operators import UFuncOperator
from .operators import SumOperator
from .operators import ListOperator
from .operators import SetOperator
from .operators import ArrayOperator
from .operators import PolyvalOperator

# Try some optimizations on these classes
try:
    import psyco
    psyco.profile()
except ImportError:
    pass
