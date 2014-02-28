#!/usr/bin/env python
########################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2011 The Trustees of Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Chris Farrow
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
########################################################################
"""Example of using simulated annealing for refining a structure.

This example uses the anneal optimizer from scipy.optimize and tight bounds on
the atom positions to determine the structure of C60.

"""

import random

import numpy

from diffpy.srfit.fitbase import Profile
from diffpy.srfit.fitbase import FitContribution, FitRecipe
from diffpy.srfit.fitbase import FitResults
from diffpy.srfit.pdf.debyepdfgenerator import DebyePDFGenerator

####### Example Code

def makeRecipe(molecule, datname):
    """Create a recipe that uses the DebyePDFGenerator."""

    ## The Profile
    profile = Profile()

    # Load data and add it to the profile
    profile.loadtxt(datname)
    profile.setCalculationRange(xmin=1.2, xmax=8)

    ## The ProfileGenerator
    # Create a DebyePDFGenerator named "G".
    generator = DebyePDFGenerator("G")
    generator.setStructure(molecule)
    # These are metadata needed by the generator
    generator.setQmin(0.68)
    generator.setQmax(22)

    ## The FitContribution
    contribution = FitContribution("bucky")
    contribution.addProfileGenerator(generator)
    contribution.setProfile(profile, xname = "r")

    # Make a FitRecipe.
    recipe = FitRecipe()
    recipe.addContribution(contribution)

    # Specify which parameters we want to refine.
    c60 = generator.phase

    # We're not going to refine the ADPs. However, every atom must have a
    # small, but finite ADP for the PDF calculator to work properly.
    atoms = c60.getScatterers()
    for atom in atoms:
        atom.Uiso.value = 0.001

    # Scale factor. We cannot optimize this efficiently, so we take the value
    # from a previous refinement. In general, care must be taken to properly
    # determine the scale of the profile, or to make sure that the residual is
    # not affected by the scale.
    generator.scale.value = 1.24457360e+4

    # Allow every atom to move.  We define the bounds to be a window of radius
    # 0.1 centered on the current value of the position.
    win = 0.1
    for idx, atom in enumerate(atoms):
        xname, yname, zname = getXYZNames(idx)
        recipe.addVar(atom.x, name = xname).boundWindow(win)
        recipe.addVar(atom.y, name = yname).boundWindow(win)
        recipe.addVar(atom.z, name = zname).boundWindow(win)

    return recipe

def plotResults(recipe):
    """Plot the results contained within a refined FitRecipe."""

    r = recipe.bucky.profile.x

    # Plot this.
    g = recipe.bucky.profile.y
    gcalc = recipe.bucky.profile.ycalc
    diffzero =  -0.8 * max(g) * numpy.ones_like(g)
    diff = g - gcalc + diffzero

    import pylab
    pylab.plot(r,g,'ob',label="G(r) Data")
    pylab.plot(r,gcalc,'-r',label="G(r) Fit")
    pylab.plot(r,diff,'-g',label="G(r) diff")
    pylab.plot(r,diffzero,'-k')
    pylab.xlabel("$r (\AA)$")
    pylab.ylabel("$G (\AA^{-2})$")
    pylab.legend(loc=1)

    pylab.show()
    return

def groupAnneal(recipe, groups):
    """Simple annealing optmizer that works with the above recipe.

    Optimize variables in related groups using scipy.optimize.anneal.  This
    function does the work of selecting groups of variables from the groups
    list at random and schedules the temperature down-ramp. When the groups
    contain atom coordinates, this essentially the Reverse Monte Carlo (RMC
    [1]) refinement approach. With constraint and restraint flexibility
    afforded by SrFit, we can include as much information about the system as
    we have.

    recipe  --  The recipe to optimize.
    groups  --  List of groups of variables or varible names. During
                a refinement step, a group is selected for perturbation. For
                RMC-like refinement, these groups represent atom positions (x,
                y, z), but in general, they can be anything.

    [1] RL McGreevy and L Pusztai, Reverse Monte Carlo Simulation: A New
        Technique for the Determination of Disordered Structures, Molecular
        Simulation 1, 359-367, 1988

    """
    from scipy.optimize import anneal

    # Maximum number of iterations
    maxiter = 10000
    # Counter for above
    niter = 0
    # How long to dwell at each temperature
    dwell = len(groups)
    # How long to wait for an improvement before exiting
    maxwait = 5 * dwell
    # Counter for the above
    waitcount = 0
    # Temperature
    T = T0 = 300
    # Current error
    err = 0
    # Minimum error
    minerr = None

    # Fix all parameters
    recipe.fix("all")

    while 1:

        if niter >= maxiter:
            print "*** Maximum interations exceeded ***"
            return
        if waitcount >= maxwait:
            print "*** Improvement waiting period exceeded ***"
            return

        niter += 1

        # Pick an atom and get the names of its positions. Free these
        # variables.
        pars = random.choice(groups)
        print "iter =", niter
        print ("%s " * len(pars)) % pars
        recipe.free(*pars)

        # Get the bounds for the variables. These are used by anneal.
        lb, ub = recipe.getBounds2()
        out = anneal(recipe, recipe.values, T0 = T,
                lower = lb, upper = ub,
                maxeval = 1, full_output = True, dwell = 1)

        err = out[1]
        if minerr is None: minerr = err

        # Record the minimum error. If we go waitcount steps without a
        # reduction in error, then we're out of here.
        if err < minerr:
            minerr = err
            waitcount = 0
        else:
            waitcount += 1

        # Update temperature
        if not niter%dwell:
            T = T0 / numpy.log(10 + niter * 1.0 / dwell)
        print "T =", T
        if out[5]:
            print "move accepted"
            reject = 0
        else:
            print "move rejected"
            reject += 1
        print "residual =", err

        # Clean up for next iteration
        recipe._applyValues(out[0])
        recipe.fix(*pars)

    return out


def main():
    """Refine the C60 molecule.

    From previous examples we know that the radius of the model structure is
    slightly too small. The annealing algorithm has to determine the proper
    radius, and must also account for thermal vibrations in the PDF.

    """

    molecule = makeC60()
    recipe = makeRecipe(molecule, "data/C60.gr")
    recipe.fithooks[0].verbose = 0

    # Optimize
    # Group related atomic positions for the optimizer
    parlist = [getXYZNames(i) for i in xrange(60)]
    groupAnneal(recipe, parlist)

    # Print results
    recipe.fix("all")
    res = FitResults(recipe, showfixed = False)
    res.printResults()

    # Save the structure
    molecule.write("C60_refined.stru", "pdffit")

    # Plot results
    plotResults(recipe)

    return

c60xyz = \
"""
3.451266498   0.685000000   0.000000000
3.451266498  -0.685000000   0.000000000
-3.451266498   0.685000000   0.000000000
-3.451266498  -0.685000000   0.000000000
0.685000000   0.000000000   3.451266498
-0.685000000   0.000000000   3.451266498
0.685000000   0.000000000  -3.451266498
-0.685000000   0.000000000  -3.451266498
0.000000000   3.451266498   0.685000000
0.000000000   3.451266498  -0.685000000
0.000000000  -3.451266498   0.685000000
0.000000000  -3.451266498  -0.685000000
3.003809890   1.409000000   1.171456608
3.003809890   1.409000000  -1.171456608
3.003809890  -1.409000000   1.171456608
3.003809890  -1.409000000  -1.171456608
-3.003809890   1.409000000   1.171456608
-3.003809890   1.409000000  -1.171456608
-3.003809890  -1.409000000   1.171456608
-3.003809890  -1.409000000  -1.171456608
1.409000000   1.171456608   3.003809890
1.409000000  -1.171456608   3.003809890
-1.409000000   1.171456608   3.003809890
-1.409000000  -1.171456608   3.003809890
1.409000000   1.171456608  -3.003809890
1.409000000  -1.171456608  -3.003809890
-1.409000000   1.171456608  -3.003809890
-1.409000000  -1.171456608  -3.003809890
1.171456608   3.003809890   1.409000000
-1.171456608   3.003809890   1.409000000
1.171456608   3.003809890  -1.409000000
-1.171456608   3.003809890  -1.409000000
1.171456608  -3.003809890   1.409000000
-1.171456608  -3.003809890   1.409000000
1.171456608  -3.003809890  -1.409000000
-1.171456608  -3.003809890  -1.409000000
2.580456608   0.724000000   2.279809890
2.580456608   0.724000000  -2.279809890
2.580456608  -0.724000000   2.279809890
2.580456608  -0.724000000  -2.279809890
-2.580456608   0.724000000   2.279809890
-2.580456608   0.724000000  -2.279809890
-2.580456608  -0.724000000   2.279809890
-2.580456608  -0.724000000  -2.279809890
0.724000000   2.279809890   2.580456608
0.724000000  -2.279809890   2.580456608
-0.724000000   2.279809890   2.580456608
-0.724000000  -2.279809890   2.580456608
0.724000000   2.279809890  -2.580456608
0.724000000  -2.279809890  -2.580456608
-0.724000000   2.279809890  -2.580456608
-0.724000000  -2.279809890  -2.580456608
2.279809890   2.580456608   0.724000000
-2.279809890   2.580456608   0.724000000
2.279809890   2.580456608  -0.724000000
-2.279809890   2.580456608  -0.724000000
2.279809890  -2.580456608   0.724000000
-2.279809890  -2.580456608   0.724000000
2.279809890  -2.580456608  -0.724000000
-2.279809890  -2.580456608  -0.724000000
"""

def getXYZNames(idx):
    """Get names for x, y, z variables based on index."""
    x = "x_%i" % idx
    y = "y_%i" % idx
    z = "z_%i" % idx
    return x, y, z

def makeC60():
    """Make the C60 molecule using diffpy.Structure."""
    from diffpy.Structure import Structure
    stru = Structure()
    for line in c60xyz.splitlines():
        if not line: continue
        xyz = map(float, line.split())
        stru.addNewAtom("C", xyz)
    return stru


if __name__ == "__main__":

    main()

# End of file
