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
"""Example of a PDF refinement using pyobjcryst and PDFGenerator.

This example is similar to crystapdf.py, except that here we refine a
pyobjcryst crystal object. In this example we use internal constraints provided
by the ObjCrystParSet structure adapter.

"""

import os
import numpy

from pyobjcryst.crystal import CreateCrystalFromCIF

from diffpy.srfit.pdf import PDFGenerator, PDFParser
from diffpy.srfit.fitbase import Profile
from diffpy.srfit.fitbase import FitContribution, FitRecipe
from diffpy.srfit.fitbase import FitResults

from gaussianrecipe import scipyOptimize
from crystalpdf import plotResults

####### Example Code

def makeRecipe(ciffile, datname):
    """Create a fitting recipe for crystalline PDF data."""

    ## The Profile
    # This will be used to store the observed and calculated PDF profile.
    profile = Profile()

    # Load data and add it to the Profile. As before we use a PDFParser. The
    # metadata is still passed to the PDFGenerator later on. The interaction
    # between the PDFGenerator and the metadata does not depend on type of
    # structure being refined.
    parser = PDFParser()
    parser.parseFile(datname)
    profile.loadParsedData(parser)
    profile.setCalculationRange(xmax = 20)

    ## The ProfileGenerator
    # This time we use the CreateCrystalFromCIF method of pyobjcryst.crystal to
    # create a Crystal object. That object is passed to the PDFGenerator as in
    # the previous example.
    generator = PDFGenerator("G")
    stru = CreateCrystalFromCIF(file(ciffile))
    generator.setPhase(stru)
    generator.setQmax(40.0)
    
    ## The FitContribution
    contribution = FitContribution("nickel")
    contribution.addProfileGenerator(generator)
    contribution.setProfile(profile, xname = "r")

    # Make the FitRecipe and add the FitContribution.
    recipe = FitRecipe()
    recipe.addContribution(contribution)

    ## Configure the fit variables

    # As before, we get a handle to the structure parameter set. In this case,
    # it is a ObjCrystParSet instance that was created when we called
    # 'setPhase' above. The ObjCrystParSet has different Parameters and options
    # than the StructureParSet used in the last example. See its documentation
    # in diffpy.srfit.structure.objcryststructure.
    phase = generator.phase

    # Here is where we created space group constraints in the previous example.
    # The difference in this example is that the ObjCrystParSet is aware of
    # space groups, and the StructureParSet is not. Constraints are created
    # internally when the ObjCrystParSet is created. (The CreateCrystalFromCIF
    # read the space group information from the cif file.) What this means is
    # that we cannot blindly constrain free Parameters from the ObjCrystParSet,
    # but must know what the space group is doing. If the ObjCrystParSet had
    # P1 symmetry, we could use the 'constrainAsSpaceGroup' method as before.
    #
    # As before, we have one free lattice parameter ('a').
    lattice = phase.getLattice()
    recipe.addVar(lattice.a)
    # And there is a free isotropic ADP. ObjCryst uses B-factors as opposed to
    # U-factors.
    Biso = recipe.newVar("Biso", 0.5)
    for scatterer in phase.getScatterers():
        recipe.constrain(scatterer.Biso, Biso)

    # We now select non-structural parameters to refine.
    # This controls the scaling of the PDF.
    recipe.addVar(generator.scale, 1)
    # This is a peak-damping resolution term.
    recipe.addVar(generator.qdamp, 0.01)
    # This is a vibrational correlation term that sharpens peaks at low-r.
    recipe.addVar(generator.delta2, 5)

    # Give the recipe away so it can be used!
    return recipe

if __name__ == "__main__":

    # Make the data and the recipe
    ciffile = "data/si.cif"
    data = "data/si-q27r60-xray.gr"

    # Make the recipe
    recipe = makeRecipe(ciffile, data)

    # Optimize
    scipyOptimize(recipe)

    # Generate and print the FitResults
    res = FitResults(recipe)
    res.printResults()

    # Plot!
    plotResults(recipe)

# End of file
