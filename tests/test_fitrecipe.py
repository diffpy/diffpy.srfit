#!/usr/bin/env python
##############################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2010 The Trustees of Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
##############################################################################
"""Tests for refinableobj module."""

import unittest

import matplotlib
import matplotlib.pyplot as plt
import pytest
from numpy import array_equal, dot, linspace, pi, sin
from scipy.optimize import leastsq

from diffpy.srfit.fitbase.fitcontribution import FitContribution
from diffpy.srfit.fitbase.fitrecipe import FitRecipe
from diffpy.srfit.fitbase.parameter import Parameter
from diffpy.srfit.fitbase.profile import Profile

matplotlib.use("Agg")


class TestFitRecipe(unittest.TestCase):

    def setUp(self):
        self.recipe = FitRecipe("recipe")
        self.recipe.fithooks[0].verbose = 0

        # Set up the Profile
        self.profile = Profile()
        x = linspace(0, pi, 10)
        y = sin(x)
        self.profile.setObservedProfile(x, y)

        # Set up the FitContribution
        self.fitcontribution = FitContribution("cont")
        self.fitcontribution.setProfile(self.profile)
        self.fitcontribution.setEquation("A*sin(k*x + c)")
        self.fitcontribution.A.setValue(1)
        self.fitcontribution.k.setValue(1)
        self.fitcontribution.c.setValue(0)

        self.recipe.addContribution(self.fitcontribution)
        return

    def testFixFree(self):
        recipe = self.recipe
        con = self.fitcontribution

        recipe.addVar(con.A, 2, tag="tagA")
        recipe.addVar(con.k, 1, tag="tagk")
        recipe.addVar(con.c, 0)
        recipe.newVar("B", 0)

        self.assertTrue(recipe.isFree(recipe.A))
        recipe.fix("tagA")
        self.assertFalse(recipe.isFree(recipe.A))
        recipe.free("tagA")
        self.assertTrue(recipe.isFree(recipe.A))
        recipe.fix("A")
        self.assertFalse(recipe.isFree(recipe.A))
        recipe.free("A")
        self.assertTrue(recipe.isFree(recipe.A))
        recipe.fix(recipe.A)
        self.assertFalse(recipe.isFree(recipe.A))
        recipe.free(recipe.A)
        self.assertTrue(recipe.isFree(recipe.A))
        recipe.fix(recipe.A)
        self.assertFalse(recipe.isFree(recipe.A))
        recipe.free("all")
        self.assertTrue(recipe.isFree(recipe.A))
        self.assertTrue(recipe.isFree(recipe.k))
        self.assertTrue(recipe.isFree(recipe.c))
        self.assertTrue(recipe.isFree(recipe.B))
        recipe.fix(recipe.A, "tagk", c=3)
        self.assertFalse(recipe.isFree(recipe.A))
        self.assertFalse(recipe.isFree(recipe.k))
        self.assertFalse(recipe.isFree(recipe.c))
        self.assertTrue(recipe.isFree(recipe.B))
        self.assertEqual(3, recipe.c.value)
        recipe.fix("all")
        self.assertFalse(recipe.isFree(recipe.A))
        self.assertFalse(recipe.isFree(recipe.k))
        self.assertFalse(recipe.isFree(recipe.c))
        self.assertFalse(recipe.isFree(recipe.B))

        self.assertRaises(ValueError, recipe.free, "junk")
        self.assertRaises(ValueError, recipe.fix, tagA=1)
        self.assertRaises(ValueError, recipe.fix, "junk")
        return

    def testVars(self):
        """Test to see if variables are added and removed properly."""
        recipe = self.recipe
        con = self.fitcontribution

        recipe.addVar(con.A, 2)
        recipe.addVar(con.k, 1)
        recipe.addVar(con.c, 0)
        recipe.newVar("B", 0)

        names = recipe.getNames()
        self.assertEqual(names, ["A", "k", "c", "B"])
        values = recipe.getValues()
        self.assertTrue((values == [2, 1, 0, 0]).all())

        # Constrain a parameter to the B-variable to give it a value
        p = Parameter("Bpar", -1)
        recipe.constrain(recipe.B, p)
        values = recipe.getValues()
        self.assertTrue((values == [2, 1, 0]).all())
        recipe.delVar(recipe.B)

        recipe.fix(recipe.k)

        names = recipe.getNames()
        self.assertEqual(names, ["A", "c"])
        values = recipe.getValues()
        self.assertTrue((values == [2, 0]).all())

        recipe.fix("all")
        names = recipe.getNames()
        self.assertEqual(names, [])
        values = recipe.getValues()
        self.assertTrue((values == []).all())

        recipe.free("all")
        names = recipe.getNames()
        self.assertEqual(3, len(names))
        self.assertTrue("A" in names)
        self.assertTrue("k" in names)
        self.assertTrue("c" in names)
        values = recipe.getValues()
        self.assertEqual(3, len(values))
        self.assertTrue(0 in values)
        self.assertTrue(1 in values)
        self.assertTrue(2 in values)
        return

    def testResidual(self):
        """Test the residual and everything that can change it."""

        # With thing set up as they are, the residual should be 0
        res = self.recipe.residual()
        self.assertAlmostEqual(0, dot(res, res))

        # Change the c value to 1 so that the equation evaluates as sin(x+1)
        x = self.profile.x
        y = sin(x + 1)
        self.recipe.cont.c.setValue(1)
        res = self.recipe.residual()
        self.assertTrue(array_equal(y - self.profile.y, res))

        # Try some constraints
        # Make c = 2*A, A = Avar
        var = self.recipe.newVar("Avar")
        self.recipe.constrain(
            self.fitcontribution.c, "2*A", {"A": self.fitcontribution.A}
        )
        self.assertEqual(2, self.fitcontribution.c.value)
        self.recipe.constrain(self.fitcontribution.A, var)
        self.assertEqual(1, var.getValue())
        self.assertEqual(self.recipe.cont.A.getValue(), var.getValue())
        # c is constrained to a constrained parameter.
        self.assertEqual(2, self.fitcontribution.c.value)
        # The equation should evaluate to sin(x+2)
        x = self.profile.x
        y = sin(x + 2)
        res = self.recipe.residual()
        self.assertTrue(array_equal(y - self.profile.y, res))

        # Now try some restraints. We want c to be exactly zero. It should give
        # a penalty of (c-0)**2, which is 4 in this case
        r1 = self.recipe.restrain(self.fitcontribution.c, 0, 0, 1)
        self.recipe._ready = False
        res = self.recipe.residual()
        chi2 = 4 + dot(y - self.profile.y, y - self.profile.y)
        self.assertAlmostEqual(chi2, dot(res, res))

        # Clear the constraint and restore the value of c to 0. This should
        # give us chi2 = 0 again.
        self.recipe.unconstrain(self.fitcontribution.c)
        self.fitcontribution.c.setValue(0)
        res = self.recipe.residual([self.recipe.cont.A.getValue()])
        chi2 = 0
        self.assertAlmostEqual(chi2, dot(res, res))

        # Remove the restraint and variable
        self.recipe.unrestrain(r1)
        self.recipe.delVar(self.recipe.Avar)
        self.recipe._ready = False
        res = self.recipe.residual()
        chi2 = 0
        self.assertAlmostEqual(chi2, dot(res, res))

        # Add constraints at the fitcontribution level.
        self.fitcontribution.constrain(self.fitcontribution.c, "2*A")
        # This should evaluate to sin(x+2)
        x = self.profile.x
        y = sin(x + 2)
        res = self.recipe.residual()
        self.assertTrue(array_equal(y - self.profile.y, res))

        # Add a restraint at the fitcontribution level.
        r1 = self.fitcontribution.restrain(self.fitcontribution.c, 0, 0, 1)
        self.recipe._ready = False
        # The chi2 is the same as above, plus 4
        res = self.recipe.residual()
        x = self.profile.x
        y = sin(x + 2)
        chi2 = 4 + dot(y - self.profile.y, y - self.profile.y)
        self.assertAlmostEqual(chi2, dot(res, res))

        # Remove those
        self.fitcontribution.unrestrain(r1)
        self.recipe._ready = False
        self.fitcontribution.unconstrain(self.fitcontribution.c)
        self.fitcontribution.c.setValue(0)
        res = self.recipe.residual()
        chi2 = 0
        self.assertAlmostEqual(chi2, dot(res, res))

        # Now try to use the observed profile inside of the equation
        # Set the equation equal to the data
        self.fitcontribution.setEquation("y")
        res = self.recipe.residual()
        self.assertAlmostEqual(0, dot(res, res))

        # Now add the uncertainty. This should give dy/dy = 1 for the residual
        self.fitcontribution.setEquation("y+dy")
        res = self.recipe.residual()
        self.assertAlmostEqual(len(res), dot(res, res))

        return


# End of class TestFitRecipe


# ----------------------------------------------------------------------------
def testPrintFitHook(capturestdout):
    "check output from default PrintFitHook."
    recipe = FitRecipe("recipe")
    recipe.fithooks[0].verbose = 0

    # Set up the Profile
    profile = Profile()
    x = linspace(0, pi, 10)
    y = sin(x)
    profile.setObservedProfile(x, y)

    # Set up the FitContribution
    fitcontribution = FitContribution("cont")
    fitcontribution.setProfile(profile)
    fitcontribution.setEquation("A*sin(k*x + c)")
    fitcontribution.A.setValue(1)
    fitcontribution.k.setValue(1)
    fitcontribution.c.setValue(0)

    recipe.addContribution(fitcontribution)

    recipe.addVar(fitcontribution.c)
    recipe.restrain("c", lb=5)
    (pfh,) = recipe.getFitHooks()
    out = capturestdout(recipe.scalarResidual)
    assert "" == out
    pfh.verbose = 1
    out = capturestdout(recipe.scalarResidual)
    assert out.strip().isdigit()
    assert "\nRestraints:" not in out
    pfh.verbose = 2
    out = capturestdout(recipe.scalarResidual)
    assert "\nResidual:" in out
    assert "\nRestraints:" in out
    assert "\nVariables" not in out
    pfh.verbose = 3
    out = capturestdout(recipe.scalarResidual)
    assert "\nVariables" in out
    assert "c = " in out
    return


def build_recipe_one_contribution():
    "helper to build a simple recipe"
    profile = Profile()
    x = linspace(0, pi, 10)
    y = sin(x)
    profile.setObservedProfile(x, y)
    contribution = FitContribution("c1")
    contribution.setProfile(profile)
    contribution.setEquation("A*sin(k*x + c)")
    recipe = FitRecipe()
    recipe.addContribution(contribution)
    recipe.addVar(contribution.A, 1)
    recipe.addVar(contribution.k, 1)
    recipe.addVar(contribution.c, 1)
    return recipe


def build_recipe_two_contributions():
    "helper to build a recipe with two contributions"
    profile1 = Profile()
    x = linspace(0, pi, 10)
    y1 = sin(x)
    profile1.setObservedProfile(x, y1)
    contribution1 = FitContribution("c1")
    contribution1.setProfile(profile1)
    contribution1.setEquation("A*sin(k*x + c)")

    profile2 = Profile()
    y2 = 0.5 * sin(2 * x)
    profile2.setObservedProfile(x, y2)
    contribution2 = FitContribution("c2")
    contribution2.setProfile(profile2)
    contribution2.setEquation("B*sin(m*x + d)")
    recipe = FitRecipe()
    recipe.addContribution(contribution1)
    recipe.addContribution(contribution2)
    recipe.addVar(contribution1.A, 1)
    recipe.addVar(contribution1.k, 1)
    recipe.addVar(contribution1.c, 1)
    recipe.addVar(contribution2.B, 0.5)
    recipe.addVar(contribution2.m, 2)
    recipe.addVar(contribution2.d, 0)

    return recipe


def optimize_recipe(recipe):
    recipe.fithooks[0].verbose = 0
    residuals = recipe.residual
    values = recipe.values
    leastsq(residuals, values)


def test_plot_recipe_bad_display():
    recipe = build_recipe_one_contribution()
    # Case: All plots are disabled
    # expected: raised ValueError with message
    plt.close("all")
    msg = "At least one of show_data, show_fit, or show_diff must be True"
    with pytest.raises(ValueError, match=msg):
        recipe.plot_recipe(
            show_data=False,
            show_diff=False,
            show_fit=False,
        )


def test_plot_recipe_no_contribution():
    recipe = FitRecipe()
    # Case: No contributions in recipe
    # expected: raised ValueError with message
    plt.close("all")
    msg = (
        "No contributions found in recipe. "
        "Add contributions before plotting."
    )
    with pytest.raises(ValueError, match=msg):
        recipe.plot_recipe()


def test_plot_recipe_before_refinement(capsys):
    # Case: User tries to plot recipe before refinement
    # expected: Data plotted without fit line or difference curve
    #          and warning message printed
    recipe = build_recipe_one_contribution()
    plt.close("all")
    before = set(plt.get_fignums())
    recipe.plot_recipe(show=False)
    after = set(plt.get_fignums())
    new_figs = after - before
    captured = capsys.readouterr()
    actual = captured.out.strip()
    expected = (
        "Contribution 'c1' has no calculated values (ycalc is None). "
        "Only observed data will be plotted."
    )
    assert len(new_figs) == 1
    assert actual == expected


def test_plot_recipe_after_refinement():
    # Case: User refines recipe and then plots
    # expected: Plot generates with no problem
    recipe = build_recipe_one_contribution()
    optimize_recipe(recipe)
    plt.close("all")
    before = set(plt.get_fignums())
    recipe.plot_recipe(show=False)
    after = set(plt.get_fignums())
    new_figs = after - before
    assert len(new_figs) == 1


def test_plot_recipe_two_contributions():
    # Case: Two contributions in recipe
    # expected: two figures created
    recipe = build_recipe_two_contributions()
    optimize_recipe(recipe)
    plt.close("all")
    before = set(plt.get_fignums())
    recipe.plot_recipe(show=False)
    after = set(plt.get_fignums())
    new_figs = after - before
    assert len(new_figs) == 2


def test_plot_recipe_on_existing_plot():
    # Case: User passes axes to plot_recipe to plot on existing figure
    # expected: User modifications are present in the final figure
    recipe = build_recipe_one_contribution()
    optimize_recipe(recipe)
    plt.close("all")
    fig, ax = plt.subplots()
    ax.set_title("User Title")
    recipe.plot_recipe(ax=ax, show=False)
    assert ax.get_title() == "User Title"


def test_plot_recipe_add_new_data():
    # Case: User wants to add data to figure generated by plot_recipe
    # Expected: New data is added to existing figure (check with labels)
    recipe = build_recipe_one_contribution()
    optimize_recipe(recipe)
    plt.close("all")
    before = set(plt.get_fignums())
    figure, ax = recipe.plot_recipe(return_fig=True, show=False)
    after = set(plt.get_fignums())
    new_figs = after - before
    # add new data to existing plot
    ax.plot([0, pi], [0, 0], label="New Data")
    ax.legend()
    legend = ax.get_legend()
    # get sorted list of legend labels for comparison
    actual_labels = sorted([t.get_text() for t in legend.get_texts()])
    expected_labels = sorted(
        ["Observed", "Calculated", "Difference", "New Data"]
    )
    assert len(new_figs) == 1
    assert actual_labels == expected_labels


def test_plot_recipe_add_new_data_two_figs():
    # Case: User wants to add data to figure generated by plot_recipe
    #       with two contributions
    # Expected: New data is added to existing figure (check with labels)
    recipe = build_recipe_two_contributions()
    optimize_recipe(recipe)
    plt.close("all")
    before = set(plt.get_fignums())
    figure, axes = recipe.plot_recipe(return_fig=True, show=False)
    after = set(plt.get_fignums())
    new_figs = after - before
    # add new data to existing plots
    for ax in axes:
        ax.plot([0, pi], [0, 0], label="New Data")
        ax.legend()
        legend = ax.get_legend()
        # get sorted list of legend labels for comparison
        actual_labels = sorted([t.get_text() for t in legend.get_texts()])
        expected_labels = sorted(
            ["Observed", "Calculated", "Difference", "New Data"]
        )
        assert actual_labels == expected_labels
    assert len(new_figs) == 1


def test_plot_recipe_set_title():
    # Case: User sets title via plot_recipe
    # Expected: Title is set correctly
    recipe = build_recipe_one_contribution()
    optimize_recipe(recipe)
    plt.close("all")
    expected_title = "Custom Recipe Title"
    figure, ax = recipe.plot_recipe(
        title=expected_title, return_fig=True, show=False
    )
    actual_title = ax.get_title()
    assert actual_title == expected_title


if __name__ == "__main__":
    unittest.main()
