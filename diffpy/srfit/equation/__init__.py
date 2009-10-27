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

"""The core equation evaluator for diffpy.srfit.

This package contains modules and subpackages that are used to create Equation
objects. An Equation is a functor that remembers its state, and that can
quickly re-evaluate its value based on changes in its variables. Equations can
be used to encapsulate simple expressions or complex signal generators.
provided one gets to know how to create and use Literal classes (from the
literals subpackage), that are the basic building blocks of an Equation.

"""

# package version
from diffpy.srfit.version import __version__

from .equationmod import Equation


# End of file
