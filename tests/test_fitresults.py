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
"""Tests for fitresults module."""

import unittest

import pytest
from scipy.optimize import leastsq

from diffpy.srfit.fitbase.fitrecipe import FitRecipe
from diffpy.srfit.fitbase.fitresults import (
    FitResults,
    initializeRecipe,
    resultsDictionary,
)

# The fit results from the recipe fixture in conftest.py
expected_fitresults = """\
My Custom header
Some quantities invalid due to missing profile uncertainty
Overall (Chi2 and Reduced Chi2 invalid)
------------------------------------------------------------------------------
Residual       0.00000000
Contributions  0.00000000
Restraints     0.00000000
Chi2           0.00000000
Reduced Chi2   0.00000000
Rw             0.00000000

Variables (Uncertainties invalid)
------------------------------------------------------------------------------
"""
expected_refined_variables = ["amplitude", "wave_number", "phase_shift"]


def optimize_recipe(recipe):
    recipe.fithooks[0].verbose = 0
    residuals = recipe.residual
    values = recipe.values
    leastsq(residuals, values)


def test_formatResults(build_recipe_one_contribution):
    recipe = build_recipe_one_contribution
    optimize_recipe(recipe)
    results = FitResults(recipe)
    actual_results_string = results.formatResults(header="My Custom header")
    # Because slight variations in refinement, just check
    # that the header of the results are the same.
    assert expected_fitresults.strip() in actual_results_string.strip()
    # check if the refined variables are in the results
    for expected_var in expected_refined_variables:
        assert expected_var in actual_results_string.strip()


def test_get_results_string(build_recipe_one_contribution):
    recipe = build_recipe_one_contribution
    optimize_recipe(recipe)
    results = FitResults(recipe)
    actual_results_string = results.get_results_string(
        header="My Custom header"
    )
    # Because slight variations in refinement, just check
    # that the header of the results are the same.
    assert expected_fitresults.strip() in actual_results_string.strip()
    # check if the refined variables are in the results
    for expected_var in expected_refined_variables:
        assert expected_var in actual_results_string.strip()


def test_printResults(build_recipe_one_contribution, capsys):
    recipe = build_recipe_one_contribution
    optimize_recipe(recipe)
    results = FitResults(recipe)
    results.printResults(header="My Custom header")
    actual_results = capsys.readouterr().out
    # Because slight variations in refinement, just check
    # that the header of the results are the same.
    assert expected_fitresults.strip() in actual_results.strip()
    # check if the refined variables are in the results
    for expected_var in expected_refined_variables:
        assert expected_var in actual_results.strip()


def test_print_results(build_recipe_one_contribution, capsys):
    recipe = build_recipe_one_contribution
    optimize_recipe(recipe)
    results = FitResults(recipe)
    results.print_results(header="My Custom header")
    actual_results = capsys.readouterr().out
    # Because slight variations in refinement, just check
    # that the header of the results are the same.
    assert expected_fitresults.strip() in actual_results.strip()
    # check if the refined variables are in the results
    for expected_var in expected_refined_variables:
        assert expected_var in actual_results.strip()


def test_saveResults(build_recipe_one_contribution, tmp_path):
    recipe = build_recipe_one_contribution
    optimize_recipe(recipe)
    results = FitResults(recipe)
    actual_results_file = tmp_path / "fit_results.txt"
    results.saveResults(actual_results_file, header="My Custom header")
    assert actual_results_file.exists()
    with open(actual_results_file, "r") as res_file:
        actual_results = res_file.read()
    # Because slight variations in refinement, just check
    # that the header of the results are the same.
    assert expected_fitresults.strip() in actual_results.strip()
    # check if the refined variables are in the results
    for expected_var in expected_refined_variables:
        assert expected_var in actual_results.strip()


def test_save_results(build_recipe_one_contribution, tmp_path):
    recipe = build_recipe_one_contribution
    optimize_recipe(recipe)
    results = FitResults(recipe)
    actual_results_file = tmp_path / "fit_results.txt"
    results.save_results(actual_results_file, header="My Custom header")
    assert actual_results_file.exists()
    with open(actual_results_file, "r") as res_file:
        actual_results = res_file.read()
    # Because slight variations in refinement, just check
    # that the header of the results are the same.
    assert expected_fitresults.strip() in actual_results.strip()
    # check if the refined variables are in the results
    for expected_var in expected_refined_variables:
        assert expected_var in actual_results.strip()


def test_get_results_dictionary(build_recipe_one_contribution):
    recipe = build_recipe_one_contribution
    optimize_recipe(recipe)
    results = FitResults(recipe)
    results_dict = results.get_results_dictionary()
    expected_results_dict = {
        "Residual": 0.0,
        "Contributions": 0.0,
        "Restraints": 0.0,
        "Chi2": 0.0,
        "Reduced Chi2": 0.0,
        "Rw": 0.0,
        "amplitude": (1.0, 4.82804000e-01),
        "phase_shift": (-1.61291146e-18, 1.00000000e00),
        "wave_number": (1.00000000e00, 2.17496687e-01),
    }
    # iterate through and assert actual_key == expected_key
    # and actual_value == expected_value (with approximation)
    for expected_key, expected_value in expected_results_dict.items():
        assert expected_key in results_dict
        actual_value = results_dict.get(expected_key)
        if isinstance(expected_value, tuple):
            # If the expected value is a tuple, check each element with approx
            assert len(actual_value) == len(expected_value)
            for actual_val, expected_val in zip(actual_value, expected_value):
                assert actual_val == pytest.approx(expected_val, rel=1e-3)
        else:
            # If the expected value is not a tuple, check with approx
            assert actual_value == pytest.approx(expected_value, rel=1e-3)


def test_resultsDictionary(build_recipe_one_contribution, tmp_path):
    recipe = build_recipe_one_contribution
    optimize_recipe(recipe)
    results = FitResults(recipe)
    results.save_results(
        tmp_path / "fit_results.txt", header="My Custom header"
    )
    actual_results_dict = resultsDictionary(str(tmp_path / "fit_results.txt"))
    expected_results_dict_keys = [
        "than",  # bad behavior: shouldn't be here
        "Feb",  # bad behavior: shouldn't be here
        "Residual",
        # "Reduced Chi2", # bad behavior: should be here but is not
        "Contributions",
        "Restraints",
        "Chi2",
        "Rw",
        "amplitude",
        "phase_shift",
        "wave_number",
    ]
    # get list of keys from actual_results_dict and check
    # that all expected keys are in it
    assert sorted(expected_results_dict_keys) == sorted(
        list(actual_results_dict.keys())
    )


def testInitializeFromFileName(datafile):
    recipe = FitRecipe("recipe")
    recipe.create_new_variable("A", 0)
    recipe.create_new_variable("sig", 0)
    recipe.create_new_variable("x0", 0)
    filename = datafile("results.res")
    Aval = 5.77619823e-01
    sigval = -9.22758690e-01
    x0val = 6.12422115e00

    assert 0 == recipe.A.value
    assert 0 == recipe.sig.value
    assert 0 == recipe.x0.value
    initializeRecipe(recipe, filename)
    assert Aval == pytest.approx(recipe.A.value)
    assert sigval == pytest.approx(recipe.sig.value)
    assert x0val == pytest.approx(recipe.x0.value)
    return


def testInitializeFromFileObj(datafile):
    recipe = FitRecipe("recipe")
    recipe.create_new_variable("A", 0)
    recipe.create_new_variable("sig", 0)
    recipe.create_new_variable("x0", 0)
    filename = datafile("results.res")
    Aval = 5.77619823e-01
    sigval = -9.22758690e-01
    x0val = 6.12422115e00

    assert 0 == recipe.A.value
    assert 0 == recipe.sig.value
    assert 0 == recipe.x0.value
    infile = open(filename, "r")
    initializeRecipe(recipe, infile)
    assert not infile.closed
    infile.close()
    assert Aval == pytest.approx(recipe.A.value)
    assert sigval == pytest.approx(recipe.sig.value)
    assert x0val == pytest.approx(recipe.x0.value)
    return


def testInitializeFromString(datafile):
    recipe = FitRecipe("recipe")
    recipe.create_new_variable("A", 0)
    recipe.create_new_variable("sig", 0)
    recipe.create_new_variable("x0", 0)
    filename = datafile("results.res")
    Aval = 5.77619823e-01
    sigval = -9.22758690e-01
    x0val = 6.12422115e00

    assert 0 == recipe.A.value
    assert 0 == recipe.sig.value
    assert 0 == recipe.x0.value
    infile = open(filename, "r")
    resstr = infile.read()
    infile.close()
    initializeRecipe(recipe, resstr)
    assert Aval == pytest.approx(recipe.A.value)
    assert sigval == pytest.approx(recipe.sig.value)
    assert x0val == pytest.approx(recipe.x0.value)
    return


if __name__ == "__main__":

    unittest.main()
