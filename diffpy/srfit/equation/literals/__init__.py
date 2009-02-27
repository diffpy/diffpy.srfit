
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


"""Literals...

Explanation of modules and how they work.
"""

# package version
from diffpy.srfit.version import __version__

# Import the operators

from .Argument import Argument
from .operators import Operator
from .operators import AdditionOperator
from .operators import SubtractionOperator
from .operators import MultiplicationOperator
from .operators import DivisionOperator
from .operators import ExponentiationOperator
from .operators import RemainderOperator
from .operators import NegationOperator
from .operators import UfuncOperator
