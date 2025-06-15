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
"""The core equation evaluator for diffpy.srfit.

This package contains modules and subpackages that are used to create
Equation objects. The Equation class is an encapsulation of a lazy
evaluation network that is used throughout SrFit. The EquationsFactory
class is used to create Equation objects from strings and can
incorporate user-defined functions as well as default operations.

The subpackages define various pieces of the evaluation network.
"""

__all__ = ["Equation"]


from diffpy.srfit.equation.equationmod import Equation

# End of file
