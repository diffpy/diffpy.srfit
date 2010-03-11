#!/usr/bin/env python
########################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2009 Trustees of the Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Chris Farrow
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
########################################################################
"""Example of simplified fitting .

This is like gaussianrecipe.py, but it uses the SimpleRecipe, which integrates
the FitContribution and Profile objects for simple recipe creation.

"""

import numpy

from diffpy.srfit.fitbase import SimpleRecipe, FitResults

####### Example Code

def main():
    """Set up a simple recipe in a few lines."""

    # The SimpleRecipe class is a type of FitRecipe
    recipe = SimpleRecipe()

    # SimpleRecipe adopts Profile methods, like 'loadtxt'
    recipe.loadtxt("data/gaussian.dat", float, '#', None, None, 0, None, False)

    # SimpleRecipe adopts FitContribution methods, like 'setEquation'
    recipe.setEquation("A * exp(-0.5*(x-x0)**2/sigma**2)")

    # SimpleRecipe has a 'vary' method that varies a parameter.
    recipe.vary('A', 0.5)
    recipe.vary('x0', 5)
    recipe.vary('sigma', 1)

    from gaussianrecipe import scipyOptimize
    scipyOptimize(recipe)

    res = FitResults(recipe)
    res.printResults()

    return

if __name__ == "__main__":

    main()

# End of file
