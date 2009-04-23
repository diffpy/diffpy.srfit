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

"""The base fit classes for diffpy.srfit.

This package contains modules and subpackages that are used to define a fit
problem in SrFit. Unaltered, these classes will help set up a fit problem that
can be optimized within a fitting framework. They provide the basic framework
for defining a forward calculator, and from that defining a fit problem with
data, constraints and restraints. The classes involved in this collaboration
can be tied to a fitting framework at various levels through inheritance. At
the very least, one must derive from FIXME to create a forward calculator,
which must also be made to conform to the requirements of the fitting
framework.

Packages:

Modules:

Classes:

Methods:

"""

# package version
from diffpy.srfit.version import __version__

from .calculator import Calculator
from .contribution import Contribution
from .fitmodel import FitModel
from .parameter import Parameter
from .profile import Profile

# End of file
