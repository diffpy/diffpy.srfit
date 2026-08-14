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

"""Example of a PDF refinement using diffpy.structure and PDFGenerator.

This is example of fitting the fcc nickel structure to measured PDF
data. The purpose of this example is to demonstrate and describe the
classes in configuration options involved with setting up a fit in this
way. The main benefit of using SrFit for PDF refinement is the
flexibility of modifying the PDF profile function for specific needs,
adding restraints to a fit and the ability to simultaneously refine a
structure to PDF data and data from other sources. This example
demonstrates only the basic configuration.
"""

import multiprocessing as mp
from pathlib import Path

from scipy.optimize import leastsq

from diffpy.srfit.fitbase import (
    FitContribution,
    FitRecipe,
    FitResults,
    Profile,
)
from diffpy.srfit.pdf import PDFGenerator, PDFParser
from diffpy.structure import Structure


def make_recipe(ciffile, datname):
    """Create a fitting recipe for crystalline PDF data."""
    # The Profile
    # This will be used to store the observed and calculated PDF profile.
    profile = Profile()

    # Load data and add it to the Profile. Unlike in other examples, we use a
    # class (PDFParser) to help us load the data. This class will read the data
    # and relevant metadata from a two- to four-column data file generated
    # with PDFGetX2 or PDFGetN. The metadata will be passed to the PDFGenerator
    # when they are associated in the FitContribution, which saves some
    # configuration steps.
    parser = PDFParser()
    parser.parse_file(datname)
    profile.load_parsed_data(parser)
    profile.set_calculation_range(xmax=20)

    # The ProfileGenerator
    # The PDFGenerator is for configuring and calculating a PDF profile. Here,
    # we want to refine a Structure object from diffpy.structure. We tell the
    # PDFGenerator that with the 'setStructure' method. All other configuration
    # options will be inferred from the metadata that is read by the PDFParser.
    # In particular, this will set the scattering type (x-ray or neutron), the
    # Qmax value, as well as initial values for the non-structural Parameters.
    generator = PDFGenerator("G")
    stru = Structure()
    stru.read(ciffile)
    generator.setStructure(stru)

    # The FitContribution
    # Here we associate the Profile and ProfileGenerator, as has been done
    # before.
    contribution = FitContribution("nickel")
    contribution.add_profile_generator(generator)
    contribution.set_profile(profile, xname="r")

    # Make the FitRecipe and add the FitContribution.
    recipe = FitRecipe()
    recipe.add_contribution(contribution)

    # Configure the fit variables

    # The PDFGenerator class holds the ParameterSet associated with the
    # Structure passed above in a data member named "phase". (We could have
    # given the ParameterSet a name other than "phase" when we added it to the
    # PDFGenerator.) The ParameterSet in this case is a StructureParameterSet,
    # the documentation for which is found in the
    # diffpy.srfit.structure.diffpystructure module.
    phase = generator.phase

    # We start by constraining the phase to the known space group. We could do
    # this by hand, but there is a method in diffpy.srfit.structure named
    # 'constrain_as_space_group' for this purpose. The constraints will by
    # default be applied to the sites, the lattice and to the ADPs.
    # See the method documentation for more details.
    # The 'constrain_as_space_group' method may create new
    # Parameters, which it returns in a SpaceGroupParameters object.
    from diffpy.srfit.structure import constrain_as_space_group

    sgpars = constrain_as_space_group(phase, "Fm-3m")

    # The SpaceGroupParameters object returned by
    # 'constrain_as_space_group' holds the free Parameters allowed by
    # the space group constraints. Once a structure is constrained,
    # we need (should) only use the Parameters
    # provided in the SpaceGroupParameters, as the relevant structure
    # Parameters are constrained to these.
    #
    # We know that the space group does not allow for any free sites because
    # each atom is on a special position. There is one free (cubic) lattice
    # parameter and one free (isotropic) ADP. We can access these Parameters in
    # the xyzpars, latpars, and adppars members of the SpaceGroupParameters
    # object.
    for par in sgpars.latpars:
        recipe.add_variable(par)
    for par in sgpars.adppars:
        recipe.add_variable(par, 0.005)

    # We now select non-structural parameters to refine.
    # This controls the scaling of the PDF.
    recipe.add_variable(generator.scale, 1)
    # This is a peak-damping resolution term.
    recipe.add_variable(generator.qdamp, 0.01)
    # This is a vibrational correlation term that sharpens peaks at low-r.
    recipe.add_variable(generator.delta2, 5)

    # Give the recipe away so it can be used!
    return recipe


def refine_recipe(recipe):
    """Helper function."""
    leastsq(recipe.residual, recipe.get_values())
    return recipe


if __name__ == "__main__":

    # Make the data and the recipe
    ciffile = str(Path(__file__).parent / "data/ni.cif")
    data = Path(__file__).parent / "data/ni-q27r100-neutron.gr"

    # sanity check
    print("==== Rw before refinements ====")
    recipe = make_recipe(ciffile, data)
    recipe.clear_fit_hooks()
    res = FitResults(recipe)
    print(res.rw)
    # Make the recipe
    recipe_list = []
    for i in range(5):
        recipe = make_recipe(ciffile, data)
        recipe.clear_fit_hooks()
        recipe_list.append(recipe)

    # Optimize
    print("==== Rw: Sequential refinements ====")
    for recipe in recipe_list:
        recipe = refine_recipe(recipe)
        res = FitResults(recipe)
        print(res.rw)

    # Make the recipe
    recipe_list = []
    for i in range(5):
        recipe = make_recipe(ciffile, data)
        recipe.clear_fit_hooks()
        recipe_list.append(recipe)

    # Optimize
    print("==== Rw: Parallel refinements ====")
    n_process = 4
    with mp.Pool(n_process) as p:
        rv = p.map(refine_recipe, recipe_list)
    for recipe in rv:
        res = FitResults(recipe)
        print(res.rw)
