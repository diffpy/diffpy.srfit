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


"""Core equation evaluator used in SrFit.

Explanation of modules and how they work.
"""

# package version
from diffpy.srfit.version import __version__

from .clicker import clickerFactory

Clicker = clickerFactory()

from .Equation import Equation

# End of file
