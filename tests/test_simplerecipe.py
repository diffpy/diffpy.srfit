#!/usr/bin/env python
##############################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2026 The Trustees of Columbia University
#                   in the City of New York.  All rights reserved.
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
##############################################################################
"""Tests for the SimpleRecipe class."""

import pytest
from numpy import arange, array_equal, linspace, ones_like

from diffpy.srfit import (
    FitContribution,
    FitResults,
    Profile,
    ProfileParser,
    SimpleRecipe,
)

XOBS = linspace(0, 1, 5)
YOBS = 2 * XOBS + 1


class _CustomContribution(FitContribution):
    pass


# A SimpleRecipe is built with its own Profile, FitContribution and
# FitResults, linked together.
@pytest.mark.parametrize(
    "input_kwargs, expected_name, expected_contribution_class",
    [
        # C1: No arguments are given.
        # Expected: The recipe is named "fit" and holds a FitContribution.
        ({}, "fit", FitContribution),
        # C2: A custom name and contribution class are given.
        # Expected: The recipe has that name and holds a contribution of
        # that class.
        (
            {"name": "custom", "conclass": _CustomContribution},
            "custom",
            _CustomContribution,
        ),
    ],
)
def test_init(input_kwargs, expected_name, expected_contribution_class):
    recipe = SimpleRecipe(**input_kwargs)
    actual_name = recipe.name
    assert actual_name == expected_name
    actual_contribution_class = type(recipe.contribution)
    assert actual_contribution_class == expected_contribution_class
    actual_contribution_name = recipe.contribution.name
    expected_contribution_name = "contribution"
    assert actual_contribution_name == expected_contribution_name
    actual_contributions = [recipe.contribution]
    expected_contributions = list(recipe._contributions.values())
    assert actual_contributions == expected_contributions
    actual_profile = recipe.contribution.profile
    expected_profile = recipe.profile
    assert actual_profile == expected_profile
    actual_profile_class = type(recipe.profile)
    expected_profile_class = Profile
    assert actual_profile_class == expected_profile_class
    actual_results_class = type(recipe.results)
    expected_results_class = FitResults
    assert actual_results_class == expected_results_class
    actual_results_recipe = recipe.results.recipe
    expected_results_recipe = recipe
    assert actual_results_recipe == expected_results_recipe
    actual_verbose = recipe.fithooks[0].verbose
    expected_verbose = 3
    assert actual_verbose == expected_verbose


# Public FitContribution methods that SimpleRecipe does not define itself
# are made available on the recipe.
@pytest.mark.parametrize(
    "input_method_name, expected_owner_name",
    [
        # C1: set_residual_equation is defined only on FitContribution.
        # Expected: The recipe method is the contribution's method.
        ("set_residual_equation", "contribution"),
        # C2: add_profile_generator is defined only on FitContribution.
        # Expected: The recipe method is the contribution's method.
        ("add_profile_generator", "contribution"),
        # C3: set_equation is defined on both classes.
        # Expected: The recipe keeps its own method.
        ("set_equation", "recipe"),
        # C4: residual is defined on both classes.
        # Expected: The recipe keeps its own method.
        ("residual", "recipe"),
    ],
)
def test_contribution_methods_are_adopted(
    input_method_name, expected_owner_name
):
    recipe = SimpleRecipe()
    owners = {"recipe": recipe, "contribution": recipe.contribution}
    actual_owner = getattr(recipe, input_method_name).__self__
    expected_owner = owners[expected_owner_name]
    assert actual_owner == expected_owner


# set_observed_profile stores the observed data in the built-in profile.
@pytest.mark.parametrize(
    "input_dyobs, expected_dyobs_is_set, expected_dy",
    [
        # C1: Uncertainties are given.
        # Expected: dyobs is set and dy is the given uncertainties.
        (0.5 * ones_like(XOBS), True, 0.5 * ones_like(XOBS)),
        # C2: No uncertainties are given.
        # Expected: dyobs stays None and dy is one at every point.
        (None, False, ones_like(XOBS)),
    ],
)
def test_set_observed_profile(input_dyobs, expected_dyobs_is_set, expected_dy):
    recipe = SimpleRecipe()
    recipe.set_observed_profile(XOBS, YOBS, input_dyobs)
    actual_xobs = recipe.profile.xobs
    expected_xobs = XOBS
    assert array_equal(actual_xobs, expected_xobs)
    actual_yobs = recipe.profile.yobs
    expected_yobs = YOBS
    assert array_equal(actual_yobs, expected_yobs)
    actual_dyobs_is_set = recipe.profile.dyobs is not None
    assert actual_dyobs_is_set == expected_dyobs_is_set
    actual_dy = recipe.profile.dy
    assert array_equal(actual_dy, expected_dy)


# set_observed_profile rejects arrays of different lengths.
@pytest.mark.parametrize(
    "input_yobs, input_dyobs",
    [
        # C1: yobs is shorter than xobs.
        # Expected: Raises a ValueError.
        (YOBS[:-1], None),
        # C2: dyobs is shorter than xobs.
        # Expected: Raises a ValueError.
        (YOBS, ones_like(XOBS)[:-1]),
    ],
)
def test_set_observed_profile_bad(input_yobs, input_dyobs):
    recipe = SimpleRecipe()
    with pytest.raises(ValueError):
        recipe.set_observed_profile(XOBS, input_yobs, input_dyobs)


# set_calculation_range and set_calculation_points set the points of the
# built-in profile. The observed points are 0, 0.25, 0.5, 0.75 and 1.
@pytest.mark.parametrize(
    "input_method_name, input_kwargs, expected_x",
    [
        # C1: The range is narrowed to [0.25, 0.75].
        # Expected: The observed points within that range.
        (
            "set_calculation_range",
            {"xmin": 0.25, "xmax": 0.75},
            [0.25, 0.5, 0.75],
        ),
        # C2: The range is set beyond the observed data.
        # Expected: The range is clipped to the observed points.
        (
            "set_calculation_range",
            {"xmin": -1, "xmax": 2},
            [0.0, 0.25, 0.5, 0.75, 1.0],
        ),
        # C3: The calculation points 0.1, 0.3 and 0.5 are given.
        # Expected: Those points are used.
        (
            "set_calculation_points",
            {"x": arange(0.1, 0.6, 0.2)},
            [0.1, 0.3, 0.5],
        ),
    ],
)
def test_set_calculation_points(input_method_name, input_kwargs, expected_x):
    recipe = SimpleRecipe()
    recipe.set_observed_profile(XOBS, YOBS)
    getattr(recipe, input_method_name)(**input_kwargs)
    actual_x = recipe.profile.x
    assert actual_x == pytest.approx(expected_x)


# loadtxt reads x, y and optionally dy from a text file into the built-in
# profile.
@pytest.mark.parametrize(
    "input_lines, expected_dyobs_is_set, expected_dy",
    [
        # C1: The file has the columns x, y and dy.
        # Expected: dyobs is set and dy is the third column.
        (
            ["1.0 2.0 0.1", "1.1 2.1 0.2", "1.2 2.2 0.3"],
            True,
            [0.1, 0.2, 0.3],
        ),
        # C2: The file has only the columns x and y.
        # Expected: dyobs stays None and dy is one at every point.
        (["1.0 2.0", "1.1 2.1", "1.2 2.2"], False, [1.0, 1.0, 1.0]),
    ],
)
def test_loadtxt(tmp_path, input_lines, expected_dyobs_is_set, expected_dy):
    data_file = tmp_path / "data.txt"
    data_file.write_text("\n".join(input_lines) + "\n")
    recipe = SimpleRecipe()
    recipe.loadtxt(data_file)
    actual_xobs = recipe.profile.xobs.tolist()
    expected_xobs = [1.0, 1.1, 1.2]
    assert actual_xobs == expected_xobs
    actual_yobs = recipe.profile.yobs.tolist()
    expected_yobs = [2.0, 2.1, 2.2]
    assert actual_yobs == expected_yobs
    actual_dyobs_is_set = recipe.profile.dyobs is not None
    assert actual_dyobs_is_set == expected_dyobs_is_set
    actual_dy = recipe.profile.dy.tolist()
    assert actual_dy == expected_dy


def test_load_parsed_data(temp_data_files):
    # C1: A two-column file is parsed by a ProfileParser and loaded.
    # Expected: The profile holds the parsed x and y and the parser's
    # metadata.
    parser = ProfileParser()
    parser.parse_file(temp_data_files / "gr_file.gr")
    recipe = SimpleRecipe()
    recipe.load_parsed_data(parser)
    actual_xobs = recipe.profile.xobs.tolist()
    expected_xobs = [1.0, 1.1, 1.2]
    assert actual_xobs == expected_xobs
    actual_yobs = recipe.profile.yobs.tolist()
    expected_yobs = [2.0, 2.1, 2.2]
    assert actual_yobs == expected_yobs
    actual_metadata = recipe.profile.meta
    expected_metadata = parser.get_metadata()
    assert actual_metadata == expected_metadata


# set_equation sets the contribution equation and makes each of its new
# parameters, other than x, y and dy, a variable of the recipe.
@pytest.mark.parametrize(
    "input_equation, expected_variables",
    [
        # C1: An equation with the new parameters A and B.
        # Expected: A and B become variables with the value 0.
        ("A*x + B", {"A": 0, "B": 0}),
        # C2: An equation that uses only x and y.
        # Expected: No variables are created.
        ("2*x + y", {}),
    ],
)
def test_set_equation(input_equation, expected_variables):
    recipe = SimpleRecipe()
    recipe.set_observed_profile(XOBS, YOBS)
    recipe.set_equation(input_equation)
    actual_variables = dict(zip(recipe.names, recipe.values))
    assert actual_variables == expected_variables


def test_set_equation_keeps_existing_variables(linear_recipe):
    # C1: The equation of a recipe with the variables A=2 and B=1 is set
    # again to one that also uses the new parameter C.
    # Expected: A and B keep their values, and C is added with the
    # value 0.
    linear_recipe.set_equation("A*x + B + C")
    actual_variables = dict(zip(linear_recipe.names, linear_recipe.values))
    expected_variables = {"A": 2.0, "B": 1.0, "C": 0}
    assert actual_variables == expected_variables


def test_call(linear_recipe):
    # C1: The recipe for A*x + B with A=2 and B=1 is called.
    # Expected: It returns 2*x + 1 on the calculation points and stores
    # it in profile.ycalc.
    actual_value = linear_recipe()
    expected_value = 2 * linear_recipe.profile.x + 1
    assert array_equal(actual_value, expected_value)
    actual_ycalc = linear_recipe.profile.ycalc
    assert array_equal(actual_ycalc, expected_value)


EXPECTED_VARIABLE_LINES = [
    "A  2.00000000e+00 +/- 1.26491106e+00",
    "B  1.00000000e+00 +/- 7.74596671e-01",
]


def test_print_results(linear_recipe, capsys):
    # C1: The results of the recipe with A=2 and B=1 are printed with a
    # header and a footer.
    # Expected: The output starts with the header, ends with the footer
    # and lists A and B with their values and uncertainties.
    linear_recipe.print_results(header="HEADER", footer="FOOTER")
    output_lines = capsys.readouterr().out.splitlines()
    actual_first_line = output_lines[0]
    expected_first_line = "HEADER"
    assert actual_first_line == expected_first_line
    actual_last_line = output_lines[-1]
    expected_last_line = "FOOTER"
    assert actual_last_line == expected_last_line
    actual_variable_lines = sorted(
        line for line in output_lines if line.startswith(("A ", "B "))
    )
    expected_variable_lines = EXPECTED_VARIABLE_LINES
    assert actual_variable_lines == expected_variable_lines


def test_save_results(linear_recipe, tmp_path):
    # C1: The results of the recipe with A=2 and B=1 are saved to a file
    # with a header and a footer.
    # Expected: The header follows the two lines that record when and by
    # whom the file was written, the file ends with the footer, and A and B
    # are listed with their values and uncertainties.
    results_file = tmp_path / "results.res"
    linear_recipe.save_results(results_file, header="HEADER", footer="FOOTER")
    saved_lines = results_file.read_text().splitlines()
    actual_header_line = saved_lines[2]
    expected_header_line = "HEADER"
    assert actual_header_line == expected_header_line
    actual_last_line = saved_lines[-1]
    expected_last_line = "FOOTER"
    assert actual_last_line == expected_last_line
    actual_variable_lines = sorted(
        line for line in saved_lines if line.startswith(("A ", "B "))
    )
    expected_variable_lines = EXPECTED_VARIABLE_LINES
    assert actual_variable_lines == expected_variable_lines
