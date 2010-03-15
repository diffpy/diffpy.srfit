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
"""Example of fitting a crystal-like nanoparticle PDF.

This is an example of modeling the PDF from a nanocrystal as an attenuated bulk
PDF. This involves a crystal PDF calculation and a nanoparticle form factor.
The equation we model is
Gnano(r) = f(r) * Gbulk(r),
where f(r) is the nanoparticle form factor (or characteristic function) for the
nanoparticle shape. Functions for calculating the nanoparticle form factor in
the diffpy.srfit.pdf.nanoformfactors module.

"""

import numpy

from pyobjcryst.crystal import CreateCrystalFromCIF

from diffpy.srfit.pdf import PDFGenerator, PDFParser
from diffpy.srfit.sas import PrCalculator, SASParser, SASGenerator
from diffpy.srfit.fitbase import Profile
from diffpy.srfit.fitbase import FitContribution, FitRecipe
from diffpy.srfit.fitbase import FitResults

from gaussianrecipe import scipyOptimize

def makeRecipe(ciffile, grdata):
    """Make a recipe to model a crystal-like nanoparticle PDF."""

    # Set up a PDF fit as has been done in other examples.
    pdfprofile = Profile()

    pdfparser = PDFParser()
    pdfparser.parseFile(grdata)
    pdfprofile.loadParsedData(pdfparser)
    pdfprofile.setCalculationRange(xmin = 0.1, xmax = 20)

    pdfcontribution = FitContribution("pdf")
    pdfcontribution.setProfile(pdfprofile, xname = "r")

    pdfgenerator = PDFGenerator("G")
    pdfgenerator.setQmax(30.0)
    stru = CreateCrystalFromCIF(file(ciffile))
    pdfgenerator.setPhase(stru)
    pdfcontribution.addProfileGenerator(pdfgenerator)

    # Register the nanoparticle shape factor.
    from diffpy.srfit.pdf.nanoformfactors import sphericalFF
    pdfcontribution.registerFunction(sphericalFF, name = "f")

    # Now we set up the fitting equation.
    pdfcontribution.setEquation("f * G")

    # Now make the recipe. Make sure we fit the form factor shape parameters,
    # in this case 'psize', which is the diameter of the particle.
    recipe = FitRecipe()
    recipe.addContribution(pdfcontribution)

    phase = pdfgenerator.phase
    lattice = phase.getLattice()
    recipe.addVar(lattice.a)
    Biso = recipe.newVar("Biso", 0.5)
    for scatterer in phase.getScatterers():
        recipe.constrain(scatterer.Biso, Biso)

    recipe.addVar(pdfcontribution.psize, 20)
    recipe.addVar(pdfgenerator.scale, 1)
    recipe.addVar(pdfgenerator.delta2, 0)

    return recipe

def plotResults(recipe):
    """Plot the results contained within a refined FitRecipe."""

    # All this should be pretty familiar by now.
    names = recipe.getNames()
    vals = recipe.getValues()

    r = recipe.pdf.profile.x

    g = recipe.pdf.profile.y
    gcalc = recipe.pdf.profile.ycalc
    offset = -0.8 * max(g)
    diff = g - gcalc + offset

    gcryst = recipe.pdf.evaluateEquation("G")
    gcryst /= recipe.scale.value

    fr = recipe.pdf.evaluateEquation("f")
    fr *= max(g) / fr[0]

    import pylab
    pylab.plot(r,g,'bo',label="G(r) Data")
    pylab.plot(r, gcryst,'y--',label="G(r) Crystal")
    pylab.plot(r, fr,'k--',label="f(r) calculated (scaled)")
    pylab.plot(r, gcalc,'r-',label="G(r) Fit")
    pylab.plot(r,numpy.ones_like(r)*offset,'k-')
    pylab.plot(r,diff,'g-',label="G(r) diff")
    pylab.xlabel("$r (\AA)$")
    pylab.ylabel("$G (\AA^{-2})$")
    pylab.legend(loc=1)

    pylab.show()
    return

if __name__ == "__main__":

    ciffile = "data/pb.cif"
    grdata = "data/pb_100_qmin1.gr"

    recipe = makeRecipe(ciffile, grdata)
    scipyOptimize(recipe)

    res = FitResults(recipe)
    res.printResults()

    plotResults(recipe)

# End of file
