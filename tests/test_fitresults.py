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

import re

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
# The tables printed for the two_sample_decay_recipe fixture in conftest.py,
# with a header and a footer
EXPECTED_TWO_SAMPLE_TABLES = """\
My Custom header
           Overall
┌───────────────┬───────────┐
│ Residual      │   172.418 │
│ Contributions │   172.418 │
│ Restraints    │         0 │
│ Chi2          │   176.886 │
│ Reduced Chi2  │   2.26777 │
│ Rw            │ 0.0911824 │
└───────────────┴───────────┘
                      Contributions
┏━━━━━━━━━━━━━━┳━━━━━━━━┳━━━━━━━━━━┳━━━━━━━━━┳━━━━━━━━━━━┓
┃ Contribution ┃ Weight ┃ Residual ┃    Chi2 ┃        Rw ┃
┡━━━━━━━━━━━━━━╇━━━━━━━━╇━━━━━━━━━━╇━━━━━━━━━╇━━━━━━━━━━━┩
│ sample_A     │      1 │   167.95 │  167.95 │   0.10299 │
│ sample_B     │    0.5 │  17.8726 │ 17.8726 │ 0.0405264 │
└──────────────┴────────┴──────────┴─────────┴───────────┘
                   Variables
┏━━━━━━━━━━┳━━━━━━━━━━┳━━━━━━━━━━━━┳━━━━━━━━━━━┓
┃ Tag      ┃ Variable ┃ Value(unc) ┃ Rel. unc. ┃
┡━━━━━━━━━━╇━━━━━━━━━━╇━━━━━━━━━━━━╇━━━━━━━━━━━┩
│ shared   │ k        │ 0.3000(70) │     2.34% │
├──────────┼──────────┼────────────┼───────────┤
│ sample A │ A_A      │  2.800(34) │     1.20% │
│          │ b_A      │  0.250(13) │     5.03% │
├──────────┼──────────┼────────────┼───────────┤
│ untagged │ A_B      │  1.400(53) │     3.77% │
└──────────┴──────────┴────────────┴───────────┘
  Fixed Variables
┏━━━━━━━━━━┳━━━━━━━┓
┃ Variable ┃ Value ┃
┡━━━━━━━━━━╇━━━━━━━┩
│ b_B      │   0.5 │
└──────────┴───────┘
 Variable Correlations greater than 25%
┏━━━━━━━━━━━━┳━━━━━━━━━━━━┳━━━━━━━━━━━━━┓
┃ Variable 1 ┃ Variable 2 ┃ Correlation ┃
┡━━━━━━━━━━━━╇━━━━━━━━━━━━╇━━━━━━━━━━━━━┩
│ k          │ b_A        │      0.6548 │
│ k          │ A_A        │      0.3832 │
│ A_B        │ k          │      0.2663 │
└────────────┴────────────┴─────────────┘
My Custom footer"""
# The tables printed for a sine recipe whose only variable is the
# untagged amplitude
EXPECTED_UNTAGGED_SINE_TABLES = """\
Some quantities invalid due to missing profile uncertainty
        Overall
┌───────────────┬──────┐
│ Residual      │ 0.18 │
│ Contributions │ 0.18 │
│ Restraints    │    0 │
│ Chi2          │ 0.18 │
│ Reduced Chi2  │ 0.02 │
│ Rw            │  0.2 │
└───────────────┴──────┘
              Variables
┏━━━━━━━━━━┳━━━━━━━━━━━━┳━━━━━━━━━━━┓
┃ Variable ┃ Value(unc) ┃ Rel. unc. ┃
┡━━━━━━━━━━╇━━━━━━━━━━━━╇━━━━━━━━━━━┩
│ A        │   1.20(47) │    39.28% │
└──────────┴────────────┴───────────┘
No correlations greater than 25%"""
# The ANSI codes rich prints around red and yellow text in a 16-color
# terminal, and after styled text
ANSI_RED = "\x1b[31m"
ANSI_YELLOW = "\x1b[33m"
ANSI_RESET = "\x1b[0m"


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


def _read_printed_tables(capsys):
    """Return the printed text without the spaces that pad its lines."""
    printed_lines = capsys.readouterr().out.splitlines()
    stripped_lines = [line.rstrip() for line in printed_lines]
    return "\n".join(stripped_lines)


def _find_colored_text(printed_text, color_code):
    """Return each piece of the printed text shown in the color with the
    ANSI code ``color_code``, without the spaces that pad it."""
    pattern = re.escape(color_code) + "(.*?)" + re.escape(ANSI_RESET)
    colored_pieces = re.findall(pattern, printed_text)
    return [piece.strip() for piece in colored_pieces]


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


def test_print_results_table_with_tagged_two_sample_fit(
    two_sample_decay_recipe, capsys, plain_terminal
):
    # Test the tables printed for a fit to two data sets.
    # C1: A user fits the decays of two samples that share a rate, adds
    # sample B's amplitude first without a tag, tags the rate and sample A's
    # variables, fixes sample B's baseline, and prints the tables with a
    # header and a footer.
    # Expected: The header comes first, then the overall metrics, the
    # metrics of each sample, the tagged variables grouped by tag in the
    # order they were added, sample B's amplitude under "untagged" after
    # them, the fixed baseline, the correlations greater than 25% strongest
    # first, and the footer.
    # Compute the results at the starting values
    results = FitResults(two_sample_decay_recipe)
    capsys.readouterr()
    # Print the tables
    results.print_results_table(
        header="My Custom header", footer="My Custom footer"
    )
    actual_printed_tables = _read_printed_tables(capsys)
    expected_printed_tables = EXPECTED_TWO_SAMPLE_TABLES
    assert actual_printed_tables == expected_printed_tables


def test_print_results_table_without_tags(
    build_sine_recipe, capsys, plain_terminal
):
    # Test the tables printed when no variable is tagged.
    # C1: A user fits sin(x) with only the amplitude free, starting at 1.2,
    # without tagging it or giving the data uncertainties, and prints the
    # tables.
    # Expected: The message about the missing uncertainties comes first,
    # the variables table has no Tag column, and a line says there are
    # no correlations greater than 25%.
    # Free the amplitude and compute the results at its starting value
    recipe = build_sine_recipe()
    recipe.add_variable(recipe.cont.A, 1.2)
    results = FitResults(recipe)
    capsys.readouterr()
    # Print the tables
    results.print_results_table()
    actual_printed_tables = _read_printed_tables(capsys)
    expected_printed_tables = EXPECTED_UNTAGGED_SINE_TABLES
    assert actual_printed_tables == expected_printed_tables


def test_print_results_table_without_show_fixed(
    two_sample_decay_recipe, capsys
):
    # Test the tables printed when fixed variables are hidden.
    # C1: A user fits the decays of two samples with sample B's baseline
    # fixed, turns off showfixed, and prints the tables.
    # Expected: The fixed baseline b_B is not listed.
    # Compute the results with the fixed variables hidden
    results = FitResults(two_sample_decay_recipe, showfixed=False)
    capsys.readouterr()
    # Print the tables
    results.print_results_table()
    actual_printed_tables = _read_printed_tables(capsys)
    actual_fixed_baseline_listed = "b_B" in actual_printed_tables
    expected_fixed_baseline_listed = False
    assert actual_fixed_baseline_listed == expected_fixed_baseline_listed


def test_print_results_table_in_color_terminal(
    short_decay_recipe, capsys, color_terminal
):
    # Test the colors that point out the results to look at.
    # C1: A user measures a decay only until it is about to level off, so
    # the rate and the baseline can't be told apart, sets the uncertainties
    # of the last three points to zero by mistake, and prints the tables in
    # a terminal that shows colors.
    # Expected: The message about the points with no weight is yellow.
    # Only the relative uncertainty greater than 10%, of b, and the
    # correlation greater than 80%, of k with b, are red.
    # Compute the results at the starting values
    results = FitResults(short_decay_recipe)
    capsys.readouterr()
    # Print the tables and find the yellow and the red text
    results.print_results_table()
    printed_text = capsys.readouterr().out
    actual_yellow_text = _find_colored_text(printed_text, ANSI_YELLOW)
    actual_red_text = _find_colored_text(printed_text, ANSI_RED)
    expected_yellow_text = [
        "3 of the 41 points of contribution 'decay' have no weight, "
        "because their uncertainties are zero, negative, infinite or not "
        "a number. To include them in the fit, please set valid "
        "uncertainties with Profile.set_observed_profile."
    ]
    expected_red_text = ["13.19%", "0.9157"]
    assert actual_yellow_text == expected_yellow_text
    assert actual_red_text == expected_red_text


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
