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
"""Example of combining PDF and SAS nanoparticles data.

This is an example of using both PDF and SAS data in the same fit. This
fits a crystal model to the PDF while fitting a shape model to both the
SAS profile and the PDF data. Using the same shape for the PDF and SAS
provides a feedback mechanism into the fit that allows the PDF and SAS
portions of the fit to guide one another, and in the end gives the shape
of the nanoparticle that agrees best with both the PDF and SAS data.
"""

import matplotlib.pyplot as plt
from gaussianrecipe import scipyOptimize
from pyobjcryst import loadCrystal

from diffpy.srfit.fitbase import (
    FitContribution,
    FitRecipe,
    FitResults,
    Profile,
)
from diffpy.srfit.pdf import PDFGenerator, PDFParser
from diffpy.srfit.pdf.characteristicfunctions import SASCF
from diffpy.srfit.sas import SASGenerator, SASParser


def makeRecipe(ciffile, grdata, iqdata):
    """Make complex-modeling recipe where I(q) and G(r) are fit
    simultaneously.

    The fit I(q) is fed into the calculation of G(r), which provides
    feedback for the fit parameters of both.
    """
    # Create a PDF contribution as before
    pdfprofile = Profile()
    pdfparser = PDFParser()
    pdfparser.parse_file(grdata)
    pdfprofile.load_parsed_data(pdfparser)
    pdfprofile.set_calculation_range(xmin=0.1, xmax=20)

    pdfcontribution = FitContribution("pdf")
    pdfcontribution.set_profile(pdfprofile, xname="r")

    pdfgenerator = PDFGenerator("G")
    pdfgenerator.setQmax(30.0)
    stru = loadCrystal(ciffile)
    pdfgenerator.setStructure(stru)
    pdfcontribution.add_profile_generator(pdfgenerator)
    pdfcontribution.set_residual_equation("resv")

    # Create a SAS contribution as well. We assume the nanoparticle is roughly
    # elliptical.
    sasprofile = Profile()
    sasparser = SASParser()
    sasparser.parse_file(iqdata)
    sasprofile.load_parsed_data(sasparser)
    if all(sasprofile.dy == 0):
        sasprofile.dy[:] = 1

    sascontribution = FitContribution("sas")
    sascontribution.set_profile(sasprofile)

    from sas.models.EllipsoidModel import EllipsoidModel

    model = EllipsoidModel()
    sasgenerator = SASGenerator("generator", model)
    sascontribution.add_profile_generator(sasgenerator)
    sascontribution.set_residual_equation("resv")

    # Now we set up a characteristic function calculator that depends on the
    # sas model.
    cfcalculator = SASCF("f", model)

    # Register the calculator with the pdf contribution and define the fitting
    # equation.
    pdfcontribution.register_calculator(cfcalculator)
    # The PDF for a nanoscale crystalline is approximated by
    # Gnano = f * Gcryst
    pdfcontribution.set_equation("f * G")

    # Moving on
    recipe = FitRecipe()
    recipe.add_contribution(pdfcontribution)
    recipe.add_contribution(sascontribution)

    # PDF
    phase = pdfgenerator.phase
    for par in phase.sgpars:
        recipe.add_variable(par)

    recipe.add_variable(pdfgenerator.scale, 1)
    recipe.add_variable(pdfgenerator.delta2, 0)

    # SAS
    recipe.add_variable(sasgenerator.scale, 1, name="iqscale")
    recipe.add_variable(sasgenerator.radius_a, 10)
    recipe.add_variable(sasgenerator.radius_b, 10)

    # Even though the cfcalculator and sasgenerator depend on the same sas
    # model, we must still constrain the cfcalculator Parameters so that it is
    # informed of changes in the refined parameters.
    recipe.add_constraint(cfcalculator.radius_a, "radius_a")
    recipe.add_constraint(cfcalculator.radius_b, "radius_b")

    return recipe


def fitRecipe(recipe):
    """We refine in stages to help the refinement converge."""
    # Tune SAS.
    recipe.set_weight(recipe.pdf, 0)
    recipe.fix("all")
    recipe.free("radius_a", "radius_b", iqscale=1e8)
    recipe.add_constraint("radius_b", "radius_a")
    scipyOptimize(recipe)
    recipe.remove_constraint("radius_b")

    # Tune PDF
    recipe.set_weight(recipe.pdf, 1)
    recipe.set_weight(recipe.sas, 0)
    recipe.fix("all")
    recipe.free("a", "Biso_0", "scale", "delta2")
    scipyOptimize(recipe)

    # Tune all
    recipe.set_weight(recipe.pdf, 1)
    recipe.set_weight(recipe.sas, 1)
    recipe.free("all")
    scipyOptimize(recipe)

    return


def plot_results(recipe):
    """Plot the results contained within a refined FitRecipe.

    The recipe has two contributions ("pdf" and "sas"), so plot_recipe
    produces one figure per contribution. The G(r) crystal and shape
    curves are not part of the standard observed/fit/diff plot, so they
    are overlaid on the "pdf" figure afterwards.
    """
    r = recipe.pdf.profile.x
    g = recipe.pdf.profile.y

    gcryst = recipe.pdf.evaluate_equation("G")
    gcryst /= recipe.scale.value

    fr = recipe.pdf.evaluate_equation("f")
    fr *= max(g) / fr[0]

    figs, axes = recipe.plot_recipe(
        show=False,
        return_fig=True,
        data_color="b",
        fit_color="r",
        diff_color="g",
        data_label="G(r) Data",
        fit_label="G(r) Fit",
        diff_label="G(r) diff",
        xlabel=r"$r (\AA)$",
        ylabel=r"$G (\AA^{-2})$",
    )
    # "pdf" was added to the recipe first, so its axes come first.
    ax = axes[0]
    ax.plot(r, gcryst, "y--", label="G(r) Crystal")
    ax.plot(r, fr, "k--", label="f(r) calculated (scaled)")
    ax.legend(loc=1)

    plt.show()
    return


if __name__ == "__main__":

    ciffile = "data/pb.cif"
    grdata = "data/pb_100_qmin1.gr"
    iqdata = "data/pb_100_qmax1.iq"

    recipe = makeRecipe(ciffile, grdata, iqdata)
    recipe.fithooks[0].verbose = 3
    fitRecipe(recipe)

    res = FitResults(recipe)
    res.print_results()

    plot_results(recipe)

# End of file
