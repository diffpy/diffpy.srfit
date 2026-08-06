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

This example uses PDFGenerator to refine a the two phase nickel-silicon
structure to all the available data.
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


def makeProfile(datafile):
    """Make an place data within a Profile."""
    profile = Profile()
    parser = PDFParser()
    parser.parse_file(datafile)
    profile.load_parsed_data(parser)
    profile.set_calculation_range(xmax=20)
    return profile


def makeContribution(name, generator, profile):
    """Make a FitContribution and add a generator and profile."""
    contribution = FitContribution(name)
    contribution.add_profile_generator(generator)
    contribution.set_profile(profile, xname="r")
    return contribution


def makeRecipe(
    ciffile_ni, ciffile_si, xdata_ni, ndata_ni, xdata_si, xdata_sini
):
    """Create a fitting recipe for crystalline PDF data."""
    # The Profiles
    # We need a profile for each data set.
    xprofile_ni = makeProfile(xdata_ni)
    xprofile_si = makeProfile(xdata_si)
    nprofile_ni = makeProfile(ndata_ni)
    xprofile_sini = makeProfile(xdata_sini)

    # The ProfileGenerators
    # We create one for each phase and share the phases.
    xgenerator_ni = PDFGenerator("xG_ni")
    stru = loadCrystal(ciffile_ni)
    xgenerator_ni.setStructure(stru)
    phase_ni = xgenerator_ni.phase

    xgenerator_si = PDFGenerator("xG_si")
    stru = loadCrystal(ciffile_si)
    xgenerator_si.setStructure(stru)
    phase_si = xgenerator_si.phase

    ngenerator_ni = PDFGenerator("nG_ni")
    ngenerator_ni.setPhase(phase_ni)

    xgenerator_sini_ni = PDFGenerator("xG_sini_ni")
    xgenerator_sini_ni.setPhase(phase_ni)

    xgenerator_sini_si = PDFGenerator("xG_sini_si")
    xgenerator_sini_si.setPhase(phase_si)

    # The FitContributions
    # We one of these for each data set.
    xcontribution_ni = makeContribution("xnickel", xgenerator_ni, xprofile_ni)
    xcontribution_si = makeContribution("xsilicon", xgenerator_si, xprofile_si)
    ncontribution_ni = makeContribution("nnickel", ngenerator_ni, nprofile_ni)
    xcontribution_sini = makeContribution(
        "xsini", xgenerator_sini_ni, xprofile_sini
    )
    xcontribution_sini.add_profile_generator(xgenerator_sini_si)
    xcontribution_sini.set_equation("scale * (xG_sini_ni +  xG_sini_si)")

    # As explained in another example, we want to minimize using Rw^2.
    xcontribution_ni.set_residual_equation("resv")
    xcontribution_si.set_residual_equation("resv")
    ncontribution_ni.set_residual_equation("resv")
    xcontribution_sini.set_residual_equation("resv")

    # Make the FitRecipe and add the FitContributions.
    recipe = FitRecipe()
    recipe.add_contribution(xcontribution_ni)
    recipe.add_contribution(xcontribution_si)
    recipe.add_contribution(ncontribution_ni)
    recipe.add_contribution(xcontribution_sini)

    # Now we vary and constrain Parameters as before.
    for par in phase_ni.sgpars:
        recipe.add_variable(par, name=par.name + "_ni")
    delta2_ni = recipe.create_new_variable("delta2_ni", 2.5)
    recipe.add_constraint(xgenerator_ni.delta2, delta2_ni)
    recipe.add_constraint(ngenerator_ni.delta2, delta2_ni)
    recipe.add_constraint(xgenerator_sini_ni.delta2, delta2_ni)

    for par in phase_si.sgpars:
        recipe.add_variable(par, name=par.name + "_si")
    delta2_si = recipe.create_new_variable("delta2_si", 2.5)
    recipe.add_constraint(xgenerator_si.delta2, delta2_si)
    recipe.add_constraint(xgenerator_sini_si.delta2, delta2_si)

    # Now the experimental parameters
    recipe.add_variable(xgenerator_ni.scale, name="xscale_ni")
    recipe.add_variable(xgenerator_si.scale, name="xscale_si")
    recipe.add_variable(ngenerator_ni.scale, name="nscale_ni")
    recipe.add_variable(xcontribution_sini.scale, 1.0, "xscale_sini")
    recipe.create_new_variable("pscale_sini_ni", 0.8)
    recipe.add_constraint(xgenerator_sini_ni.scale, "pscale_sini_ni")
    recipe.add_constraint(xgenerator_sini_si.scale, "1 - pscale_sini_ni")

    # The qdamp parameters are too correlated to vary so we fix them based on
    # previous measurements.
    xgenerator_ni.qdamp.value = 0.055
    xgenerator_si.qdamp.value = 0.051
    ngenerator_ni.qdamp.value = 0.030
    xgenerator_sini_ni.qdamp.value = 0.052
    xgenerator_sini_si.qdamp.value = 0.052

    # Give the recipe away so it can be used!
    return recipe


def plot_results(recipe):
    """Plot the results contained within a refined FitRecipe.

    The recipe has four contributions ("xnickel", "xsilicon", "nnickel",
    "xsini"), so plot_recipe produces one figure per contribution.
    """
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
    ciffile_ni = "data/ni.cif"
    ciffile_si = "data/si.cif"
    xdata_ni = "data/ni-q27r60-xray.gr"
    ndata_ni = "data/ni-q27r100-neutron.gr"
    xdata_si = "data/si-q27r60-xray.gr"
    xdata_sini = "data/si90ni10-q27r60-xray.gr"

    # Make the recipe
    recipe = makeRecipe(
        ciffile_ni, ciffile_si, xdata_ni, ndata_ni, xdata_si, xdata_sini
    )

    # Optimize
    scipyOptimize(recipe)

    # Generate and print the FitResults
    res = FitResults(recipe)
    res.print_results()

    # Plot!
    plot_results(recipe)

# End of file
