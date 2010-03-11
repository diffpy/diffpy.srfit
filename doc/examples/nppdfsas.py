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
"""Example of combining PDF and SAS nanoparticles data. 

This is an example of using both PDF and SAS data in the same fit. This example
refines a crystal structure (Pb) to nanoparticle data by approximating the PDF
as 
Gnano(r) = f(r) * Gcrystal(r).
(See Acta Cryst. A, 65 p. 232 (2009) and references therein.) This equation
assumes that the nanoparticle is crystal-like. (The hypothetical structure from
which the PDF and SAS data were calculated, data/pb_100.xyz, is indeed a cut
from a crystal.)

We use the PDFGenerator to calculate Gcrystal(r), and retrieve f(r) from the
SAS data. The PrCalculator class can calculate P(r), the nanoparticle density
factor, from the nanoparticle I(Q), using the Invertor object from
sans.pr.invertor. P(r) is related to f(r) via
P(r) = 4 * pi * rho * r**2 f(r).
P(r) is the radial distribution function (RDF) for particle with uniform
density.  The Invertor class performs an indirect transform of I(Q) to obtain
P(r). See the class documentation for more details.

Below we use both a PDFGenerator and PrCalculator to calculate Gnano(r) and
refine the crystal structure of lead to fit the nanoparticle PDF data.

"""

import numpy

from pyobjcryst.crystal import CreateCrystalFromCIF

from diffpy.srfit.pdf import PDFGenerator, PDFParser
from diffpy.srfit.sas import PrCalculator, SASParser, SASGenerator
from diffpy.srfit.fitbase import Profile
from diffpy.srfit.fitbase import FitContribution, FitRecipe
from diffpy.srfit.fitbase import FitResults

from gaussianrecipe import scipyOptimize

def makeRecipe(ciffile, grdata, iqdata):
    """Make a recipe to combine PDF and SAS data in a nanoparticle fit."""

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
    pdfgenerator.setQmax(30.0)
    pdfcontribution.addProfileGenerator(pdfgenerator)

    # Load in the SAS data for the nanoparticle. For convenience, we hold this
    # in a Profile.
    sasprofile = Profile()
    sasparser = SASParser()
    sasparser.parseFile(iqdata)
    sasprofile.loadParsedData(sasparser)

    # Create the PrCalculator. The PrCalculator requires parameters q, iq and
    # diq to be specified. These represent the SAS Q, I(Q) and uncertainty in
    # I(Q). These are held in the sasprofile.
    prcalculator = PrCalculator("P")
    prcalculator.q.setValue(sasprofile.x)
    prcalculator.iq.setValue(sasprofile.y)
    prcalculator.diq.setValue(sasprofile.dy)

    # Now we register the calculator with pdfcontribution. This allows us to
    # use it in the fitting equation. The nanoparticle fitting equation is 
    # f(r) * G(r), where f(r) = P(r) / (4 * pi * r**2).
    pdfcontribution.registerCalculator(prcalculator)
    pdfcontribution.setEquation("P/(4 * pi * r**2) * G")

    # Now we can move on as before. The use of the PrCalculator does not add
    # any fittable parameters to the 
    recipe = FitRecipe()
    recipe.addContribution(pdfcontribution)

    phase = pdfgenerator.phase
    lattice = phase.getLattice()
    recipe.addVar(lattice.a)
    Biso = recipe.newVar("Biso", 0.5)
    for scatterer in phase.getScatterers():
        recipe.constrain(scatterer.Biso, Biso)

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

    pr = recipe.pdf.evaluateEquation("P")
    fr = pr / r**2
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
    iqdata = "data/pb_100_qmax1.1.iq"

    recipe = makeRecipe(ciffile, grdata, iqdata)
    scipyOptimize(recipe)

    res = FitResults(recipe)
    res.printResults()

    plotResults(recipe)

# End of file
