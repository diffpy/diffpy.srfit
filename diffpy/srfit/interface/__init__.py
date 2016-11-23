#!/usr/bin/env python
##############################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2010 The Trustees of Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Chris Farrow
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
##############################################################################

"""Usability interface for SrFit.

The classes and functions in this package are designed to unobtrusively mix
with base SrFit objects and provide them with interface enhancements for
scripting.
"""


from diffpy.srfit.interface.interface import ParameterInterface
_parameter_interface = ParameterInterface

from diffpy.srfit.interface.interface import RecipeOrganizerInterface
_recipeorganizer_interface = RecipeOrganizerInterface

from diffpy.srfit.interface.interface import FitRecipeInterface
_fitrecipe_interface = FitRecipeInterface


# End of file
