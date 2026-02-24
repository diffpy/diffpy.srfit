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
from diffpy.srfit.fitbase.fitresults import FitResults, initializeRecipe

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
