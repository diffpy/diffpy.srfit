#!/usr/bin/env python
########################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2010 The Trustees of Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Chris Farrow
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
########################################################################
"""Example of a refinement of SAS I(Q) data to an ellipsoidal model."""

import matplotlib.pyplot as plt
from gaussianrecipe import scipyOptimize

from diffpy.srfit.fitbase import (
    FitContribution,
    FitRecipe,
    FitResults,
    Profile,
)
from diffpy.srfit.sas import SASGenerator, SASParser

######
#  Example Code


def makeRecipe(datname):
    """Create a fitting recipe for ellipsoidal SAS data."""
    # The Profile
    # This will be used to store the observed and calculated I(Q) data.
    profile = Profile()

    # Load data and add it to the Profile. We use a SASParser to load the data
    # properly and pass the metadata along.
    parser = SASParser()
    parser.parse_file(datname)
    profile.load_parsed_data(parser)

    # The ProfileGenerator
    # The SASGenerator is for configuring and calculating a SAS profile. We use
    # a sas model to configure and serve as the calculation engine of the
    # generator. This allows us to use the full sas model creation
    # capabilities, and tie this into SrFit when we want to fit a model to
    # data. The documentation for the various sas models can be found at
    # http://www.sasview.org.
    from sas.models.EllipsoidModel import EllipsoidModel

    model = EllipsoidModel()
    generator = SASGenerator("generator", model)

    # The FitContribution
    # Here we associate the Profile and ProfileGenerator, as has been done
    # before.
    contribution = FitContribution("ellipsoid")
    contribution.add_profile_generator(generator)
    contribution.set_profile(profile, xname="q")

    # We want to fit the log of the signal to the log of the data so that the
    # higher-Q information remains significant. There are no I(Q) uncertainty
    # values with the data, so we do not need to worry about the effect this
    # will have on the estimated parameter uncertainties.
    contribution.set_residual_equation("log(eq) - log(y)")

    # Make the FitRecipe and add the FitContribution.
    recipe = FitRecipe()
    recipe.add_contribution(contribution)

    # Configure the fit variables
    # The SASGenerator uses the parameters from the params and dispersion
    # attributes of the model. These vary from model to model, but are adopted
    # as SrFit Parameters within the generator. Whereas the dispersion
    # parameters are accessible as, e.g. "radius.width", within the
    # SASGenerator these are named like "radius_width".
    #
    # We want to fit the scale factor, radii and background factors.
    recipe.add_variable(generator.scale, 1)
    recipe.add_variable(generator.radius_a, 50)
    recipe.add_variable(generator.radius_b, 500)
    recipe.add_variable(generator.background, 0)

    # Give the recipe away so it can be used!
    return recipe


def plot_results(recipe):
    """Plot the results contained within a refined FitRecipe."""
    # I(Q) SAS data is best viewed on a log-log scale, so we prepare the
    # axes ourselves and hand them to plot_recipe. The residual difference
    # curve is not shown since it can go negative and is not meaningful on
    # a log scale.
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xscale("log")
    ax.set_yscale("log")
    recipe.plot_recipe(
        ax=ax,
        show=False,
        show_diff=False,
        data_color="b",
        fit_color="r",
        data_label="I(Q) Data",
        fit_label="I(Q) Fit",
        xlabel=r"$Q (\AA^{-1})$",
        ylabel="$I (arb. units)$",
    )
    plt.show()
    return


if __name__ == "__main__":

    # Make the data and the recipe
    data = "data/sas_ellipsoid_testdata.txt"

    # Make the recipe
    recipe = makeRecipe(data)

    # Optimize
    scipyOptimize(recipe)

    # Generate and print the FitResults
    res = FitResults(recipe)
    res.print_results()

    # Plot!
    plot_results(recipe)

# End of file
