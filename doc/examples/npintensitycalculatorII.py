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
"""Example of extracting information from multiple data simultaneously.

This example builds on npintensitycalculator.py, and uses IntensityCalculator
from that example, and simultaneoulsy refines a recipe to two data sets
generated from the same structure.

The makeRecipe function shows how to build a FitRecipe that uses the
IntensityCalculator and refines it to two data sets.

"""

import os

import numpy

from diffpy.srfit.fitbase import FitContribution, FitRecipe, Profile, FitResults
from npintensitycalculator import IntensityCalculator, scipyOptimize
from npintensitycalculator import iofq, parkOptimize, makeData

from gaussianrecipe import scipyOptimize, parkOptimize

####### Example Code

def makeRecipe(strufile, datname1, datname2):
    """Create a recipe that uses the IntensityCalculator.

    We will create two FitContributions that use the IntensityCalculator from
    npintensitycalculator.py and associate each of these with a Profile, and
    use this to define a FitRecipe.

    Both simulated data sets come from the same structure. We're going to make
    two FitContributions that are identical, except for the profile that is
    held in each. We're going to assure that the structures are identical by
    using the same StructureParSet (which is generated by the
    IntensityCalculator when we load the structure) in both calculators.

    """

    ## The Profiles
    # Create two Profiles for the two FitContributions.
    profile1 = Profile()
    profile2 = Profile()

    # Load data into the Profiles
    x, y, u = numpy.loadtxt(datname1, unpack=True)
    profile1.setObservedProfile(x, y, u)
    x, y, u = numpy.loadtxt(datname2, unpack=True)
    profile2.setObservedProfile(x, y, u)

    ## The Calculators
    # Create two IntensityCalculators named "I". There will not be a name
    # conflict, since the name is only meaningful within the FitContribution
    # that holds the Calculator.  Load the structure into one and make sure
    # that the second calculator is using the same StructureParSet.  This will
    # assure that both Calculators are using the exact same Parameters, and
    # therefore the same structural information.
    calculator1 = IntensityCalculator("I")
    calculator1.setStructure(strufile)
    calculator2 = IntensityCalculator("I")
    calculator2.addParameterSet(calculator1.structure)
    
    ## The FitContributions
    # Create the FitContributions. Note that the "q" in each of the FitContributions
    # refers to that FitContributions' Profile.
    contribution1 = FitContribution("bucky1")
    contribution1.setCalculator(calculator1)
    contribution1.setProfile(profile1, xname = "q")
    contribution2 = FitContribution("bucky2")
    contribution2.setCalculator(calculator2)
    contribution2.setProfile(profile2, xname = "q")

    # Now we're ready to define the fitting equation for each FitContribution.
    # The functions registered below will be independent, even though they take
    # the same form and use the same Parameter names.  By default, Parameters
    # in different contributions are different Parameters even if they have the
    # same names.  FitContributions are isolated namespaces than only share
    # information if you tell them to by using addParameter or addParameterSet.
    #
    # Unlike last time, we'll use the "polyval" function from numpy.
    bkgdstr = "b0 + b1*q + b2*q**2 + b3*q**3 + b4*q**4 + b5*q*5 + b6*q**6 +\
               b7*q**7 +b8*q**8 + b9*q**9"

    contribution1.registerStringFunction(bkgdstr, "bkgd")
    contribution2.registerStringFunction(bkgdstr, "bkgd")

    # We will create the broadening function by registering a python function.
    pi = numpy.pi
    exp = numpy.exp
    def gaussian(q, q0, width):
        return 1/(2*pi*width**2)**0.5 * exp(-0.5 * ((q-q0)/width)**2)

    contribution1.registerFunction(gaussian)
    contribution2.registerFunction(gaussian)
    # Center the gaussian
    contribution1.q0.setValue(x[len(x)/2])
    contribution2.q0.setValue(x[len(x)/2])

    # Now we can incorporate the scale and bkgd into our calculation. We also
    # convolve the signal with the gaussian to broaden it.
    contribution1.setEquation("scale * convolve(I, gaussian) + bkgd")
    contribution2.setEquation("scale * convolve(I, gaussian) + bkgd")

    # Make a FitRecipe and associate the FitContributions.
    recipe = FitRecipe()
    recipe.addContribution(contribution1)
    recipe.addContribution(contribution2)

    # Specify which Parameters we want to refine. We want to refine the
    # background that we just defined in the FitContributions. We have to do
    # this separately for each fitcontribution.
    recipe.addVar(contribution1.b0, 0, name = "b1_0")
    recipe.addVar(contribution1.b1, 0, name = "b1_1")
    recipe.addVar(contribution1.b2, 0, name = "b1_2")
    recipe.addVar(contribution1.b3, 0, name = "b1_3")
    recipe.addVar(contribution1.b4, 0, name = "b1_4")
    recipe.addVar(contribution1.b5, 0, name = "b1_5")
    recipe.addVar(contribution1.b6, 0, name = "b1_6")
    recipe.addVar(contribution1.b7, 0, name = "b1_7")
    recipe.addVar(contribution1.b8, 0, name = "b1_8")
    recipe.addVar(contribution1.b9, 0, name = "b1_9")
    recipe.addVar(contribution2.b0, 0, name = "b2_0")
    recipe.addVar(contribution2.b1, 0, name = "b2_1")
    recipe.addVar(contribution2.b2, 0, name = "b2_2")
    recipe.addVar(contribution2.b3, 0, name = "b2_3")
    recipe.addVar(contribution2.b4, 0, name = "b2_4")
    recipe.addVar(contribution2.b5, 0, name = "b2_5")
    recipe.addVar(contribution2.b6, 0, name = "b2_6")
    recipe.addVar(contribution2.b7, 0, name = "b2_7")
    recipe.addVar(contribution2.b8, 0, name = "b2_8")
    recipe.addVar(contribution2.b9, 0, name = "b2_9")

    # We also want to adjust the scale and the convolution width
    recipe.addVar(contribution1.scale, 1, name = "scale1")
    recipe.addVar(contribution1.width, 0.1, name = "width1")
    recipe.addVar(contribution2.scale, 1, name = "scale2")
    recipe.addVar(contribution2.width, 0.1, name = "width2")

    # We can also refine structural parameters. We only have to do this once,
    # since each calculator holds the same StructureParSet.
    structure = calculator1.structure
    a = recipe.addVar(structure.lattice.a)
    # We want to allow for isotropic expansion, so we'll make constraints for
    # that.
    recipe.constrain(structure.lattice.b, a)
    recipe.constrain(structure.lattice.c, a)
    # We want to refine the thermal paramters as well. We will add a new
    # variable that we call "Uiso" and constrain the atomic Uiso values to
    # this.
    Uiso = recipe.newVar("Uiso", 0.004)
    for atom in structure.atoms:
        recipe.constrain(atom.Uiso, Uiso)

    # Give the recipe away so it can be used!
    return recipe

def plotResults(recipe):
    """Plot the results contained within a refined FitRecipe."""

    # For the basic info about the fit, we can use the FitRecipe directly
    names = recipe.getNames()
    vals = recipe.getValues()

    q = recipe.bucky1.profile.x

    # Plot this for fun.
    I1 = recipe.bucky1.profile.y
    Icalc1 = recipe.bucky1.profile.ycalc
    bkgd1 = recipe.bucky1.evaluateEquation("bkgd")
    diff1 = I1 - Icalc1
    I2 = recipe.bucky2.profile.y
    Icalc2 = recipe.bucky2.profile.ycalc
    bkgd2 = recipe.bucky2.evaluateEquation("bkgd")
    diff2 = I2 - Icalc2
    offset = 1.2*max(I2)
    I1 += offset
    Icalc1 += offset
    bkgd1 += offset
    diff1 += offset

    import pylab
    pylab.plot(q,I1,'o',label="I1(Q) Data")
    pylab.plot(q,Icalc1,label="I1(Q) Fit")
    pylab.plot(q,diff1,label="I1(Q) diff")
    pylab.plot(q,bkgd1,label="Bkgd.1 Fit")
    pylab.plot(q,I2,'o',label="I2(Q) Data")
    pylab.plot(q,Icalc2,label="I2(Q) Fit")
    pylab.plot(q,diff2,label="I2(Q) diff")
    pylab.plot(q,bkgd2,label="Bkgd.2 Fit")
    pylab.xlabel("$Q (\AA^{-1})$")
    pylab.ylabel("Intensity (arb. units)")
    pylab.legend(loc=1)

    pylab.show()
    return

if __name__ == "__main__":

    # Make two different data sets, each from the same structure, but with
    # different scale, noise, broadening and background.
    strufile = "data/C60.stru"
    q = numpy.arange(1, 20, 0.05)
    makeData(strufile, q, "C60_1.iq", 8.1, 100.68, 0.005, 0.12, 2, 0.01)
    makeData(strufile, q, "C60_2.iq", 8.8, 100.68, 0.005, 0.003, 0, 1)

    # Make the recipe
    recipe = makeRecipe(strufile, "C60_1.iq", "C60_2.iq")

    # Optimize
    scipyOptimize(recipe)
    #parkOptimize(recipe)

    # Generate and print the FitResults
    res = FitResults(recipe)
    res.printResults()

    # Plot!
    plotResults(recipe)


# End of file
