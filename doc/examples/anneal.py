#!/usr/bin/env python
########################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2011 Trustees of the Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Chris Farrow
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
########################################################################
"""Example of using simulated annealing for refining a structure.

This example uses the anneal optimizer from scipy.optimize and tight bounds on
the atom positions to determine the structure of C60.

"""

import random

import numpy
import scipy

from diffpy.srfit.fitbase import ProfileGenerator, Profile
from diffpy.srfit.fitbase import FitContribution, FitRecipe
from diffpy.srfit.fitbase import FitResults
from diffpy.srfit.pdf.debyepdfgenerator import DebyePDFGenerator
from diffpy.srfit.pdf import PDFGenerator

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
        atom.Biso.value = 0.001

    # Scale factor. We cannot optimize this efficiently, so we take the value
    # from another optimization. In general, care must be taken to properly
    # determine the scale of the profile, or to make sure that the residual is
    # not affected by the scale.
    generator.scale.value = 1.24457360e+04

    # Allow every atom to move
    for idx, atom in enumerate(atoms):
        # We put tight bounds on the atoms since the optimizer will not
        # converge otherwise.
        var = recipe.addVar(atom.x, name = "x_%i"%idx)
        val = var.value
        var.bounds = (val - 0.1, val + 0.1)
        var = recipe.addVar(atom.y, name = "y_%i"%idx)
        val = var.value
        var.bounds = (val - 0.1, val + 0.1)
        var = recipe.addVar(atom.z, name = "z_%i"%idx)
        val = var.value
        var.bounds = (val - 0.1, val + 0.1)

    return recipe

def plotResults(recipe):
    """Plot the results contained within a refined FitRecipe."""

    names = recipe.getNames()
    vals = recipe.getValues()

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

def atomAnneal(recipe, natoms):
    """Simple annealing optmizer that works with the above recipe.

    Optimize the atomic positions of the recipe using scipy.optimize.anneal.
    This function does the work of selecting atoms at random and schedules the
    temperature down-ramp. Since this is not a general purpose annealing
    optimizer, the scheduling is rather crude.

    It is assumed that the atomic positions in the structure are variables in
    the recipe and are named "x_i", "y_i" and "z_i", where 'i' is the atom
    index.

    recipe  --  The recipe to optimize.
    natoms  --  The number of atoms in the structure to optimize. This is used
                to index atom positions.
    
    """
    from scipy.optimize import anneal

    # Maximum number of iterations
    maxiter = 10000
    # Counter for above
    niter = 0
    # How long to dwell at each temperature
    dwell = natoms
    # How long to wait for an improvement before exiting
    maxwait = dwell
    # Counter for the above
    waitcount = 0
    # Temperature
    T = T0 = 300
    # Current error
    err = 0
    # Minimum error
    minerr = None

    def getVarNames(idx):
        x = "x_%i" % idx
        y = "y_%i" % idx
        z = "z_%i" % idx
        return x, y, z

    def updateTemp():
        T = T0 / numpy.log(1 + niter * 1.0 / dwell)

    # Fix all parameters
    recipe.fix("all")

    while niter < maxiter and waitcount < maxwait:

        niter += 1

        # Pick an atom and get the names of its positions. Free these
        # variables.
        idx = random.randint(0, natoms-1)
        x, y, z = getVarNames(idx)
        print "iter =", niter
        print x, y, z
        recipe.free(x, y, z)

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
        if out[5]: 
            print "move accepted"
            reject = 0
        else:
            print "move rejected"
            reject += 1
        if not niter%dwell:
            T = T0 / numpy.log(10 + niter * 1.0 / dwell)
        print "T =", T
        print "residual =", err

        # Clean up for next iteration
        recipe._applyValues(out[0])
        recipe.fix(x, y, z)

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
    atomAnneal(recipe, len(molecule))

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
