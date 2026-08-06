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
"""Example of a PDF refinement of two-phase structure.

Like the ones before, this example uses PDFGenerator to refine a
structure to PDF data. However, for a multi-phase structure one must use
multiple PDFGenerators. This example refines a physical mixture of
nickel and silicon to find the structures and phase fractions.
"""

from gaussianrecipe import scipyOptimize
from pyobjcryst import loadCrystal

from diffpy.srfit.fitbase import (
    FitContribution,
    FitRecipe,
    FitResults,
    Profile,
)
from diffpy.srfit.pdf import PDFGenerator, PDFParser

######
#  Example Code


def makeRecipe(niciffile, siciffile, datname):
    """Create a fitting recipe for crystalline PDF data."""
    # The Profile
    profile = Profile()

    # Load data and add it to the profile
    parser = PDFParser()
    parser.parse_file(datname)
    profile.load_parsed_data(parser)
    profile.set_calculation_range(xmax=20)

    # The ProfileGenerator
    # In order to fit two phases simultaneously, we must use two PDFGenerators.
    # PDFGenerator is designed to take care of as little information as it
    # must. (Don't do too much, and do it well.) A PDFGenerator can generate
    # the signal from only a single phase at a time. So, we will create one
    # PDFGenerator for each phase and compose them within the same
    # FitContribution. Note that both generators will be associated with the
    # same Profile within the FitContribution, so they will both be
    # automatically configured according to the metadata.
    #
    # The generator for the nickel phase. We call it "G_ni" and will use this
    # name later when we set the fitting equation in the FitContribution.
    generator_ni = PDFGenerator("G_ni")
    stru = loadCrystal(niciffile)
    generator_ni.setStructure(stru)
    # The generator for the silicon phase. We call it "G_si".
    generator_si = PDFGenerator("G_si")
    stru = loadCrystal(siciffile)
    generator_si.setStructure(stru)

    # The FitContribution
    # Add both generators to the FitContribution. Add the Profile. This will
    # send the metadata to the generators.
    contribution = FitContribution("nisi")
    contribution.add_profile_generator(generator_ni)
    contribution.add_profile_generator(generator_si)
    contribution.set_profile(profile, xname="r")

    # Write the fitting equation. We want to sum the PDFs from each phase and
    # multiply it by a scaling factor. We also want a certain phase scaling
    # relationship between the PDFs which we will enforce with constraints in
    # the FitRecipe.
    contribution.set_equation("scale * (G_ni +  G_si)")

    # Make the FitRecipe and add the FitContribution.
    recipe = FitRecipe()
    recipe.add_contribution(contribution)

    # Configure the fit variables
    # Start by configuring the scale factor and resolution factors.
    # We want the sum of the phase scale factors to be 1.
    recipe.create_new_variable("scale_ni", 0.1)
    recipe.add_constraint(generator_ni.scale, "scale_ni")
    recipe.add_constraint(generator_si.scale, "1 - scale_ni")
    # We also want the resolution factor to be the same on each.
    recipe.create_new_variable("qdamp", 0.03)
    recipe.add_constraint(generator_ni.qdamp, "qdamp")
    recipe.add_constraint(generator_si.qdamp, "qdamp")

    # Vary the global scale as well.
    recipe.add_variable(contribution.scale, 1)

    # Now we can configure the structural parameters. Since we're using
    # ObjCrystCrystalParSets, the space group constraints are automatically
    # applied to each phase. We must selectively vary the free parameters.
    #
    # First the nickel parameters
    phase_ni = generator_ni.phase
    for par in phase_ni.sgpars:
        recipe.add_variable(par, name=par.name + "_ni")
    recipe.add_variable(generator_ni.delta2, name="delta2_ni")
    # Next the silicon parameters
    phase_si = generator_si.phase
    for par in phase_si.sgpars:
        recipe.add_variable(par, name=par.name + "_si")
    recipe.add_variable(generator_si.delta2, name="delta2_si")

    # We have prior information from the earlier examples so we'll use it here
    # in the form of restraints.
    #
    # The nickel lattice parameter was measured to be 3.527. The uncertainty
    # values are invalid for that measurement, since the data from which it is
    # derived has no uncertainty. Thus, we will tell the recipe to scale the
    # residual, which means that it will be weighted as much as the average
    # data point during the fit.
    recipe.add_soft_bounds(
        "a_ni", lower_bound=3.527, upper_bound=3.527, scaled=True
    )
    # Now we do the same with the delta2 and Biso parameters (remember that
    # Biso = 8*pi**2*Uiso)
    recipe.add_soft_bounds(
        "delta2_ni", lower_bound=2.22, upper_bound=2.22, scaled=True
    )
    recipe.add_soft_bounds(
        "Biso_0_ni", lower_bound=0.454, upper_bound=0.454, scaled=True
    )
    #
    # We can do the same with the silicon values. We haven't done a thorough
    # job of measuring the uncertainties in the results, so we'll scale these
    # as well.
    recipe.add_soft_bounds(
        "a_si", lower_bound=5.430, upper_bound=5.430, scaled=True
    )
    recipe.add_soft_bounds(
        "delta2_si", lower_bound=3.54, upper_bound=3.54, scaled=True
    )
    recipe.add_soft_bounds(
        "Biso_0_si", lower_bound=0.645, upper_bound=0.645, scaled=True
    )

    # Give the recipe away so it can be used!
    return recipe


def plot_results(recipe):
    """Plot the results contained within a refined FitRecipe."""
    recipe.plot_recipe(
        data_color="b",
        fit_color="r",
        diff_color="g",
        data_label="G(r) Data",
        fit_label="G(r) Fit",
        diff_label="G(r) diff",
        xlabel=r"$r (\AA)$",
        ylabel=r"$G (\AA^{-2})$",
    )
    return


if __name__ == "__main__":

    # Make the data and the recipe
    niciffile = "data/ni.cif"
    siciffile = "data/si.cif"
    data = "data/si90ni10-q27r60-xray.gr"

    # Make the recipe
    recipe = makeRecipe(niciffile, siciffile, data)

    # Optimize
    scipyOptimize(recipe)

    # Generate and print the FitResults
    res = FitResults(recipe)
    res.print_results()

    # Plot!
    plot_results(recipe)

# End of file
