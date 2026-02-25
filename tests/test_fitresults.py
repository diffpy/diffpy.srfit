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

import numpy as np
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
    results.print_results()
    actual_results_dict = results.get_results_dictionary()
    expected_results_dict = {
        "amplitude": 1.000000000060171,
        "wave_number": 1.00000000012548,
        "phase_shift": -1.6129114631049646e-18,
        "Residual": 3.3284672708760557e-19,
        "Contributions": 3.3284672708760557e-19,
        "Restraints": 0,
        "Chi2": 3.3284672708760557e-19,
        "Reduced Chi2": 4.7549532441086507e-20,
        "Rw": 2.7196679825449506e-10,
    }
    actual_values = np.round(np.array(list(actual_results_dict.values())), 5)
    actual_keys = set(actual_results_dict.keys())
    expected_values = np.round(
        np.array(list(expected_results_dict.values())), 5
    )
    expected_keys = set(expected_results_dict.keys())
    assert expected_keys == actual_keys
    assert list(expected_values == list(actual_values))


def test_resultsDictionary(temp_data_files):
    actual_results_dict = resultsDictionary(
        temp_data_files / "fit_results.res"
    )
    # bad behavior: values are stored as strings
    expected_results_dict = {
        "than": "25",  # bad behavior: shouldn't be here
        "wave_number": "1.00000000e+00",
        "phase_shift": "-1.61291146e-18",
        "amplitude": "1.00000000e+00",
        "Rw": "0.00000000",
        "Chi2": "0.00000000",
        "Restraints": "0.00000000",
        "Contributions": "0.00000000",
        "Residual": "0.00000000",
        "Feb": "25",  # bad behavior: shouldn't be here
    }
    # convert values to float for comparison (with rounding)
    for key in expected_results_dict:
        expected_results_dict[key] = float(expected_results_dict[key])
    for key in actual_results_dict:
        actual_results_dict[key] = float(actual_results_dict[key])

    actual_keys = set(actual_results_dict.keys())
    actual_values = np.round(np.array(list(actual_results_dict.values())), 5)
    expected_keys = set(expected_results_dict.keys())
    expected_values = np.round(
        np.array(list(expected_results_dict.values())), 5
    )
    assert expected_keys == actual_keys
    assert list(expected_values == list(actual_values))


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
