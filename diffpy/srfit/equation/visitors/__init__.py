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
"""Visitors that perform on Literal trees.

Visitors are designed to traverse and extract infromation from Literal trees
(diffpy.srfit.equation.literals). Visitors are useful for validating, printing
and extracting Argument literals from the Literal tree.

The Literal-Visitor relationship is that described by the Visitor pattern
(http://en.wikipedia.org/wiki/Visitor_pattern), execept that Literals contain
attributes that are used specifically by the Evaluator visitor.

"""

# package version
from diffpy.srfit.version import __version__

from .argfinder import ArgFinder, getArguments
from .printer import Printer
from .validator import Validator
