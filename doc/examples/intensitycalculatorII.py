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
"""Example of using Calculators in FitModels.

This is an example of building a FitModel in order to fit theoretical intensity
data. To understand the basics of FitModels, see the debyemodel.py.

The IntensityCalculator class is an example of a Calculator that can be used by
a Contribution to help generate a signal.

The makeModel function shows how to build a FitModel that uses the
IntensityCalculator.  The scipyOptimize and parkOptimize functions show two
different ways of refining the model, using scipy and park, respectively.
"""

import os

import numpy

from diffpy.srfit.fitbase import Contribution, FitModel, Profile
from intensitycalculator import IntensityCalculator, scipyOptimize
from intensitycalculator import iofq, parkOptimize, makeData

from debyemodel import scipyOptimize, parkOptimize

####### Example Code

def makeModel(strufile, datname1, datname2):
    """Create a model that uses the IntensityCalculator.

    This will create a two Contributions that use the IntensityCalculator,
    associate this with a Profile, and use this to define a FitModel.

    The simulated data comes from the same structure. We're going to make two
    contributions that are identical, except for the profile that is held in
    each. We're going to assure that the structures are identical by using the
    same StructureParSet in both calculators.

    """

    ## The Profiles
    # Create a Profile2.
    profile1 = Profile()
    profile2 = Profile()

    # Load data and add it to the profile
    x, y, u = numpy.loadtxt(datname1, unpack=True)
    profile1.setObservedProfile(x, y, u)
    x, y, u = numpy.loadtxt(datname2, unpack=True)
    profile2.setObservedProfile(x, y, u)

    ## The Calculators
    # Create two IntensityCalculators named "I".
    # Load the structure into one and set the second calculator's structure
    # equal to this one.
    calculator1 = IntensityCalculator("I")
    calculator1.setStructure(strufile)
    calculator2 = IntensityCalculator("I")
    calculator2.addParameterSet(calculator1.structure)
    
    ## The Contributions
    # Create the Contributions.
    contribution1 = Contribution("bucky1")
    contribution1.setCalculator(calculator1)
    contribution1.setProfile(profile1, xname = "q")
    contribution2 = Contribution("bucky2")
    contribution2.setCalculator(calculator2)
    contribution2.setProfile(profile2, xname = "q")

    # Now we're ready to define the contribution equation. We need to modify
    # each Calcultor. The functions registered below will be independent, even
    # though they take the same form and use the same parameter names.  By
    # default, parameters in different contributions are different parameters
    # even if they have the same names.

    # We will define the background as a string.
    bkgdstr = "b0 + b1*q + b2*q**2 + b3*q**3 + b4*q**4 + b5*q*5 + b6*q**6 +\
               b7*q**7 +b8*q**8 + b9*q**9"

    contribution1.registerStringFunction(bkgdstr, "bkgd")
    contribution2.registerStringFunction(bkgdstr, "bkgd")

    # We will create the broadening function by registering a python function.
    pi = numpy.pi
    exp = numpy.exp
    def gaussian(q, q0, width):
        return 1/(width*(2*pi)**0.5) * exp(-0.5 * ((q-q0)/width)**2)

    contribution1.registerFunction(gaussian)
    contribution2.registerFunction(gaussian)
    # Center the gaussian
    contribution1.q0.setValue(x[len(x)/2])
    contribution2.q0.setValue(x[len(x)/2])

    # Now we can incorporate the scale and bkgd into our calculation. We also
    # convolve the signal with the gaussian to broaden it.
    contribution1.setEquation("scale * convolve(I, gaussian) + bkgd")
    contribution2.setEquation("scale * convolve(I, gaussian) + bkgd")

    # Make a FitModel where we can create variables, constraints and
    # restraints. If we had multiple profiles to fit simultaneously, the
    # contribution from each could be added to the model.
    model = FitModel()
    model.addContribution(contribution1)
    model.addContribution(contribution2)

    # Specify which parameters we want to refine. We can give them initial
    # values in the process. We want to refine the background variables that we
    # just defined in the contribution. We have to do this separately for each
    # contribution.
    model.addVar(contribution1.b0, 0, name = "b1_0")
    model.addVar(contribution1.b1, 0, name = "b1_1")
    model.addVar(contribution1.b2, 0, name = "b1_2")
    model.addVar(contribution1.b3, 0, name = "b1_3")
    model.addVar(contribution1.b4, 0, name = "b1_4")
    model.addVar(contribution1.b5, 0, name = "b1_5")
    model.addVar(contribution1.b6, 0, name = "b1_6")
    model.addVar(contribution1.b7, 0, name = "b1_7")
    model.addVar(contribution1.b8, 0, name = "b1_8")
    model.addVar(contribution1.b9, 0, name = "b1_9")
    model.addVar(contribution2.b0, 0, name = "b2_0")
    model.addVar(contribution2.b1, 0, name = "b2_1")
    model.addVar(contribution2.b2, 0, name = "b2_2")
    model.addVar(contribution2.b3, 0, name = "b2_3")
    model.addVar(contribution2.b4, 0, name = "b2_4")
    model.addVar(contribution2.b5, 0, name = "b2_5")
    model.addVar(contribution2.b6, 0, name = "b2_6")
    model.addVar(contribution2.b7, 0, name = "b2_7")
    model.addVar(contribution2.b8, 0, name = "b2_8")
    model.addVar(contribution2.b9, 0, name = "b2_9")

    # We also want to adjust the scale and the convolution width
    model.addVar(contribution1.scale, 1, name = "scale1")
    model.addVar(contribution1.width, 0.1, name = "width1")
    model.addVar(contribution2.scale, 1, name = "scale2")
    model.addVar(contribution2.width, 0.1, name = "width2")

    # We can also refine structural parameters. We only have to do these
    # equations for one of the calculators, since they hold the same structure
    # instance.
    structure = calculator1.structure
    a = structure.lattice.a
    model.addVar(a)
    # We want to allow for isotropic expansion, so we'll make constraints for
    # that.
    model.constrain(structure.lattice.b, a)
    model.constrain(structure.lattice.c, a)
    # We want to refine the thermal paramters as well. We will add a new
    # variable that we call "Uiso" and constrain the atomic Uiso values to
    # this.
    Uiso = model.newVar("Uiso", 0.004)
    for atom in structure.atoms:
        model.constrain(atom.Uiso, Uiso)

    # Give the model away so it can be used!
    return model

def displayResults(model):
    """Display the results contained within a refined FitModel."""

    # For the basic info about the fit, we can use the FitModel directly
    chiv = model.residual()
    names = model.getNames()
    vals = model.getValues()

    q = model.bucky1.profile.x

    print "Fit using the scipy LM optimizer"
    chi2 = numpy.dot(chiv, chiv)
    rchi2 = chi2 / (len(q) - len(vals))
    print "Chi^2 =", chi2
    print "Red. Chi^2 =", rchi2
    
    for name, val in zip(names, vals):
        print "%s = %f" % (name, val)

    # Plot this for fun.
    I1 = model.bucky1.profile.y
    Icalc1 = model.bucky1.profile.ycalc
    bkgd1 = model.bucky1.evaluateEquation("bkgd()")
    diff1 = I1 - Icalc1
    I2 = model.bucky2.profile.y
    Icalc2 = model.bucky2.profile.ycalc
    bkgd2 = model.bucky2.evaluateEquation("bkgd()")
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

    # Make the data and the model
    strufile = "data/C60.stru"
    q = numpy.arange(1, 20, 0.05)
    # Make two different data sets, each from the same structure, but with
    # different noise, broadening and background.
    makeData(strufile, q, "C60_1.iq", 8, 10.068, 0.005, 0.12, 2)
    makeData(strufile, q, "C60_2.iq", 8, 10.068, 0.005, 0.003, 0)

    model = makeModel(strufile, "C60_1.iq", "C60_2.iq")
    scipyOptimize(model)
    displayResults(model)

    model = makeModel(strufile, "C60_1.iq", "C60_1.iq")
    parkOptimize(model)
    displayResults(model)

# End of file
