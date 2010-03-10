########################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2010 Trustees of the Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Chris Farrow
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
########################################################################

"""Usability interface for SrFit.

The classes and functions in this package are designed to unobtrusively mix
with base SrFit objects and provide them with interface enhancements for
scripting. To use the interface, import this module and call the 'use' method
before importing any other part of SrFit.

"""

_parameter_interface = object
_recipeorganizer_interface = object
_fitcontribution_interface = object
_fitrecipe_interface = object

from .interface import ParameterFactory

# package version
from diffpy.srfit.version import __version__

def use():
    """Use the interface enhancements."""
    from .interface import ParameterInterface
    global _parameter_interface
    _parameter_interface = ParameterInterface

    from .interface import RecipeOrganizerInterface
    global _recipeorganizer_interface
    _recipeorganizer_interface = RecipeOrganizerInterface

    from .interface import FitContributionInterface
    global _fitcontribution_interface
    _fitcontribution_interface = FitContributionInterface

    from .interface import FitRecipeInterface
    global _fitrecipe_interface
    _fitrecipe_interface = FitRecipeInterface

    return

__id__ = "$Id$"
# End of file
