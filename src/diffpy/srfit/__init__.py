#!/usr/bin/env python
##############################################################################
#
# (c) 2008-2025 The Trustees of Columbia University in the City of New York.
# All rights reserved.
#
# File coded by: Christopher Farrow, Pavol Juhas, and members of the
# Billinge Group.
#
# See GitHub contributions for a more detailed list of contributors.
# https://github.com/diffpy/diffpy.srfit/graphs/contributors
#
# See LICENSE.rst for license information.
#
##############################################################################
"""Complex modeling framework for structure refinement and solution.

SrFit is a tool for coherently combining known information about a
material to derive other properties, in particular material structure.
SrFit allows the customization and creation of structure
representations, profile calculators, constraints, restraints and file
input parsers. The customized pieces can be glued together within SrFit
to optimize a structure, or other physically relevant information from
one or more experimental profiles. Other known information about the
system of interest can be included with arbitrarily complex constraints
and restraints. In this way, the end user creates a customized fitting
application that suits the problem to the available information.

The subpackages herein define various pieces of the SrFit framework.
Developers are encouraged to work through the examples described in the
documentation to learn how to use and customize the various parts of
SrFit.
"""

# package version
from diffpy.srfit.version import __version__

# silence the pyflakes syntax checker
assert __version__ or True

# End of file
