#!/usr/bin/env python
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
"""Example of using ProfileGenerators in FitContributions.

This is an example of building a ProfileGenerator and using it in a
FitContribution in order to fit theoretical intensity data.

The IntensityGenerator class is an example of a ProfileGenerator that can be
used by a FitContribution to help generate a signal.

The makeRecipe function shows how to build a FitRecipe that uses the
IntensityGenerator.

"""

import os

import numpy

from diffpy.Structure import Structure
from diffpy.srfit.pdf import PDFGenerator
from diffpy.srfit.fitbase import Profile
from diffpy.srfit.fitbase import FitContribution, FitRecipe
from diffpy.srfit.fitbase import FitResults
from diffpy.srfit.structure.diffpystructure import StructureParSet

from gaussianrecipe import scipyOptimize, parkOptimize

####### Example Code

def makeRecipe(strufile, datname):
    """Create a recipe that uses the IntensityGenerator.

    This will create a FitContribution that uses the IntensityGenerator,
    associate this with a Profile, and use this to define a FitRecipe.

    """

    ## The Profile
    profile = Profile()

    # Load data and add it to the profile
    x, y, junk, u = numpy.loadtxt(datname, unpack=True)
    profile.setObservedProfile(x, y, u)
    profile.setCalculationRange(0, 20, 0.05)

    ## The ProfileGenerator
    generator = PDFGenerator("G")
    stru = Structure()
    stru.read(strufile)
    generator.setPhase(stru)
    
    ## The FitContribution
    contribution = FitContribution("nickel")
    contribution.addProfileGenerator(generator)
    contribution.setProfile(profile, xname = "r")

    # Make the FitRecipe and add the FitContribution.
    recipe = FitRecipe()
    recipe.addContribution(contribution)

    phase = generator.phase

    recipe.addVar(generator.scale, 1)
    recipe.addVar(generator.qdamp, 0.01)
    recipe.addVar(generator.delta2, 5)
    lattice = phase.getLattice()
    a = lattice.a
    recipe.addVar(a)
    recipe.constrain(lattice.b, a)
    recipe.constrain(lattice.c, a)
    # We want to refine the thermal paramters as well. We will add a new
    # Variable that we call "Uiso" and constrain the atomic Uiso values to
    # this.
    Uiso = recipe.newVar("Uiso")
    for scatterer in phase.getScatterers():
        recipe.constrain(scatterer.Uiso, Uiso)

    # Give the recipe away so it can be used!
    return recipe

def plotResults(recipe):
    """Plot the results contained within a refined FitRecipe."""

    # All this should be pretty familiar by now.
    names = recipe.getNames()
    vals = recipe.getValues()

    r = recipe.nickel.profile.x

    g = recipe.nickel.profile.y
    gcalc = recipe.nickel.profile.ycalc
    diff = g - gcalc - 0.5 * max(g)

    import pylab
    pylab.plot(r,g,'bo',label="G(r) Data")
    pylab.plot(r, gcalc,'r-',label="G(r) Fit")
    pylab.plot(r,diff,'g-',label="G(r) diff")
    pylab.xlabel("$r (\AA)$")
    pylab.ylabel("$G (\AA^{-2})$")
    pylab.legend(loc=1)

    pylab.show()
    return

if __name__ == "__main__":

    # Make the data and the recipe
    strufile = "data/ni.stru"
    data = "data/ni.dat"

    # Make the recipe
    recipe = makeRecipe(strufile, data)

    # Optimize
    scipyOptimize(recipe)
    #parkOptimize(recipe)

    # Generate and print the FitResults
    res = FitResults(recipe)
    res.printResults()

    # Plot!
    plotResults(recipe)

# End of file
