
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

# Make a clicker for this package
from diffpy.srfit.equation.clicker import clickerFactory
Clicker = clickerFactory()

# Import the operators

from diffpy.srfit.equation.literals.Argument import Argument
from diffpy.srfit.equation.literals.operators import Operator
from diffpy.srfit.equation.literals.operators import AdditionOperator
from diffpy.srfit.equation.literals.operators import SubtractionOperator
from diffpy.srfit.equation.literals.operators import MultiplicationOperator
from diffpy.srfit.equation.literals.operators import DivisionOperator
from diffpy.srfit.equation.literals.operators import ExponentiationOperator
from diffpy.srfit.equation.literals.operators import RemainderOperator
