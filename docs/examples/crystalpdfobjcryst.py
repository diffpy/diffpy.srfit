#!/usr/bin/env python
########################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2009 The Trustees of Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Chris Farrow
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
########################################################################
"""Example of a PDF refinement using pyobjcryst and PDFGenerator.

This example is similar to crystalpdf.py, except that here we refine a
pyobjcryst crystal object. In this example we use internal constraints
provided by the ObjCrystCrystalParSet structure adapter.
"""

from crystalpdf import plotResults
from gaussianrecipe import scipyOptimize
from pyobjcryst import loadCrystal

from diffpy.srfit.fitbase import (
    FitContribution,
    FitRecipe,
    FitResults,
    Profile,
    ProfileParser,
)
from diffpy.srfit.pdf import PDFGenerator

######
#  Example Code


def makeRecipe(ciffile, datname):
    """Create a fitting recipe for crystalline PDF data."""

    # The Profile
    # This will be used to store the observed and calculated PDF profile.
    profile = Profile()

    # Load data and add it to the Profile. As before we use a ProfileParser.
    # The metadata is still passed to the PDFGenerator later on.
    # The interaction between the PDFGenerator and the metadata does not
    # depend on type of structure being refined.
    parser = ProfileParser()
    parser.parseFile(datname)
    profile.load_parsed_data(parser)
    profile.set_calculation_range(xmax=20)

    # The ProfileGenerator
    # This time we use the CreateCrystalFromCIF method of pyobjcryst.crystal to
    # create a Crystal object. That object is passed to the PDFGenerator as in
    # the previous example.
    generator = PDFGenerator("G")
    stru = loadCrystal(ciffile)
    generator.setStructure(stru)
    generator.setQmax(40.0)

    # The FitContribution
    contribution = FitContribution("nickel")
    contribution.add_profile_generator(generator)
    contribution.set_profile(profile, xname="r")

    # Make the FitRecipe and add the FitContribution.
    recipe = FitRecipe()
    recipe.add_contribution(contribution)

    # Configure the fit variables

    # As before, we get a handle to the structure parameter set. In this case,
    # it is a ObjCrystCrystalParSet instance that was created when we called
    # 'setStructure' above. The ObjCrystCrystalParSet has different Parameters
    # and options than the DiffpyStructureParSet used in the last example. See
    # its documentation in diffpy.srfit.structure.objcrystparset.
    phase = generator.phase

    # Here is where we created space group constraints in the previous example.
    # The difference in this example is that the ObjCrystCrystalParSet is aware
    # of space groups, and the DiffpyStructureParSet is not. Constraints are
    # created internally when "sgpars" attribute is called for. These
    # constraints get enforced within the ObjCrystCrystalParSet. Free
    # Parameters are stored within the 'sgpars' member of the
    # ObjCrystCrystalParSet, which is the same as the object returned from
    # 'constrainAsSpaceGroup'.
    #
    # As before, we have one free lattice parameter ('a'). We can simplify
    # things by iterating through all the sgpars.
    for par in phase.sgpars:
        recipe.add_variable(par)
    # set the initial thermal factor to a non-zero value
    assert hasattr(recipe, "B11_0")
    recipe.B11_0 = 0.1

    # We now select non-structural parameters to refine.
    # This controls the scaling of the PDF.
    recipe.add_variable(generator.scale, 1)
    # This is a peak-damping resolution term.
    recipe.add_variable(generator.qdamp, 0.01)
    # This is a vibrational correlation term that sharpens peaks at low-r.
    recipe.add_variable(generator.delta2, 5)

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
    res.print_results()

    # Plot!
    plotResults(recipe)

# End of file
