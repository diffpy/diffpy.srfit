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

"""The base fitting classes for diffpy.srfit.

This package contains modules and subpackages that are used to define a fit
problem in SrFit. Unaltered, these classes will help set up a fit problem that
can be optimized within a fitting framework or with a standalone optimizer.
They provide the basic framework for defining a forward calculator, and from
that defining a fit problem with data, constraints and restraints. The classes
involved in this collaboration can be tied to a fitting framework at various
levels through inheritance. One can create a fitting problem using the
FitContribution, Profile and FitRecipe classes.

Various code and design taken from Paul Kienzle's PARK package.
http://www.reflectometry.org/danse/park.html

"""

__all__ = ['Calculator', 'FitContribution', 'FitHook', 'FitRecipe',
'FitResults', 'initalizeRecipe', 'PlotFitHook', 'Profile', 'ProfileGenerator']

# package version
from diffpy.srfit.version import __version__

from .calculator import Calculator
from .fitcontribution import FitContribution
from .fithook import FitHook, PlotFitHook
from .fitrecipe import FitRecipe
from .simplerecipe import SimpleRecipe
from .fitresults import FitResults, initializeRecipe
from .profile import Profile
from .profilegenerator import ProfileGenerator

__id__ = "$Id$"
# End of file
