#!/usr/bin/env python
##############################################################################
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
##############################################################################

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
           'FitResults', 'initializeRecipe', 'PlotFitHook', 'Profile',
           'ProfileGenerator', 'SimpleRecipe']

from diffpy.srfit.fitbase.calculator import Calculator
from diffpy.srfit.fitbase.fitcontribution import FitContribution
from diffpy.srfit.fitbase.fithook import FitHook, PlotFitHook
from diffpy.srfit.fitbase.fitrecipe import FitRecipe
from diffpy.srfit.fitbase.simplerecipe import SimpleRecipe
from diffpy.srfit.fitbase.fitresults import FitResults, initializeRecipe
from diffpy.srfit.fitbase.profile import Profile
from diffpy.srfit.fitbase.profilegenerator import ProfileGenerator

# End of file
