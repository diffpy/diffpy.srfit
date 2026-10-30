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
"""Tests for the fitresults module."""

import pytest
from numpy import array, dot, inf, isfinite, linspace, pi, sin, sqrt
from scipy.optimize import leastsq

from diffpy.srfit import FitResults

# The start of the fit results of the recipe fixture in conftest.py
EXPECTED_RESULTS_HEADER = """\
My Custom header
Some quantities invalid due to missing profile uncertainty
Overall (Chi2 and Reduced Chi2 invalid)
------------------------------------------------------------------------------
Residual       0.00000000
Contributions  0.00000000
Restraints     0.00000000
Chi2           0.00000000
Reduced Chi2   0.00000000
Rw             0.00000010

Variables (Uncertainties invalid)
------------------------------------------------------------------------------
"""
EXPECTED_REFINED_VARIABLES = ["amplitude", "wave_number", "phase_shift"]


def _optimize_recipe(recipe):
    recipe.fithooks[0].verbose = 0
    leastsq(recipe.residual, recipe.values)


def _get_results_output(results, method_name, capsys, tmp_path):
    """Return the results text given by get_results_string, printed by
    print_results or saved by save_results."""
    if method_name == "get_results_string":
        return results.get_results_string(header="My Custom header")
    if method_name == "print_results":
        results.print_results(header="My Custom header")
        return capsys.readouterr().out
    results_file = tmp_path / "fit_results.txt"
    results.save_results(results_file, header="My Custom header")
    return results_file.read_text()


# Test the results text of a refined recipe, whether returned, printed or
# saved.
@pytest.mark.parametrize(
    "input_method_name",
    [
        # C1: The results text is returned by get_results_string.
        # Expected: It contains the header, the overall metrics and the
        # refined variables.
        "get_results_string",
        # C2: The results text is printed by print_results.
        # Expected: It contains the header, the overall metrics and the
        # refined variables.
        "print_results",
        # C3: The results text is saved to a file by save_results.
        # Expected: It contains the header, the overall metrics and the
        # refined variables.
        "save_results",
    ],
)
def test_results_output(
    build_recipes_one_contribution, capsys, tmp_path, input_method_name
):
    recipe, _ = build_recipes_one_contribution
    _optimize_recipe(recipe)
    results = FitResults(recipe)
    capsys.readouterr()
    actual_results_output = _get_results_output(
        results, input_method_name, capsys, tmp_path
    ).strip()
    # Refinements vary slightly, so only the header and the variable names
    # are checked.
    assert EXPECTED_RESULTS_HEADER.strip() in actual_results_output
    actual_variables_present = [
        name in actual_results_output for name in EXPECTED_REFINED_VARIABLES
    ]
    expected_variables_present = [True] * len(EXPECTED_REFINED_VARIABLES)
    assert actual_variables_present == expected_variables_present


def test_get_results_dictionary(build_recipes_one_contribution):
    # C1: The results dictionary is taken after the recipe is refined to
    # the noiseless sine it is fit to.
    # Expected: It holds the refined value of each variable and overall
    # metrics that are all zero within 1e-5.
    recipe, _ = build_recipes_one_contribution
    _optimize_recipe(recipe)
    results = FitResults(recipe)
    actual_results_dict = results.get_results_dictionary()
    expected_results_dict = {
        "amplitude": recipe.amplitude.value,
        "wave_number": recipe.wave_number.value,
        "phase_shift": recipe.phase_shift.value,
        "Residual": 0,
        "Contributions": 0,
        "Restraints": 0,
        "Chi2": 0,
        "Reduced Chi2": 0,
        "Rw": 0,
    }
    assert actual_results_dict == pytest.approx(
        expected_results_dict, abs=1e-5
    )
    actual_amplitude = actual_results_dict["amplitude"]
    expected_amplitude = 1
    assert actual_amplitude == pytest.approx(expected_amplitude, abs=1e-5)


# Points whose uncertainties are zero, negative, infinite or not a number
# are left out of Chi2, Reduced Chi2 and Rw, as if they weren't in the
# data, and the results say how many there are. The recipe compares
# A*sin(k*x + c), with A = 0.8, k = 1.2 and c = 0.1, to sin(x) on 21 points
# in [0, pi], without refining. The three variables are counted in Reduced
# Chi2.
@pytest.mark.parametrize(
    "input_dyobs, expected_dy, expected_messages",
    [
        # C1: A user fits data whose uncertainties are all valid.
        # Expected: Every point counts, and there are no messages.
        ([0.1] * 21, [0.1] * 21, []),
        # C2: A user fits data whose last three uncertainties are zero, as
        # at the end of a detector bank.
        # Expected: Only the other 18 points count, and a message says how
        # many points have no weight and how to set their uncertainties.
        (
            [0.1] * 18 + [0.0] * 3,
            [0.1] * 18 + [inf] * 3,
            [
                "3 of the 21 points of contribution 'c1' have no weight, "
                "because their uncertainties are zero, negative, infinite "
                "or not a number. To include them in the fit, please set "
                "valid uncertainties with Profile.set_observed_profile."
            ],
        ),
        # C3: A user fits data whose uncertainties are all zero.
        # Expected: Every point counts with an uncertainty of one, and the
        # message for missing uncertainties is given.
        (
            [0.0] * 21,
            [1.0] * 21,
            ["Some quantities invalid due to missing profile uncertainty"],
        ),
    ],
)
def test_results_with_invalid_uncertainties(
    build_recipe_with_uncertainty, input_dyobs, expected_dy, expected_messages
):
    # Compute the results at the starting values
    recipe, profile = build_recipe_with_uncertainty(input_dyobs)
    results = FitResults(recipe)
    results.get_results_string()
    actual_messages = results.messages
    actual_metrics = [results.chi2, results.rchi2, results.rw]
    # Compute the expected metrics from the points with weight only
    x = linspace(0, pi, 21)
    y = sin(x)
    dy = array(expected_dy)
    is_weighted = isfinite(dy)
    chiv = (0.8 * sin(1.2 * x + 0.1) - y)[is_weighted] / dy[is_weighted]
    weighted_y = y[is_weighted] / dy[is_weighted]
    expected_chi2 = dot(chiv, chiv)
    expected_rchi2 = expected_chi2 / (is_weighted.sum() - 3)
    expected_rw = sqrt(expected_chi2 / dot(weighted_y, weighted_y))
    expected_metrics = [expected_chi2, expected_rchi2, expected_rw]
    assert actual_messages == expected_messages
    assert actual_metrics == pytest.approx(expected_metrics)
