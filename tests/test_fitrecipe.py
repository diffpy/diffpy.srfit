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
"""Tests for the FitRecipe class."""

import re

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pytest
from bg_mpl_stylesheets.styles import all_styles
from numpy import linspace, ones, pi, sin
from scipy.optimize import leastsq

from diffpy.srfit import (
    FitContribution,
    FitRecipe,
    FitResults,
    Profile,
    ProfileParser,
)
from diffpy.srfit.fitbase.parameter import Parameter
from diffpy.srfit.fitbase.parameterset import ParameterSet

matplotlib.use("Agg")


# Fix and then free a single variable, selecting it by tag, by name or by
# the variable itself. A variable is free when it is added.
@pytest.mark.parametrize(
    "input_selection",
    [
        # C1: Variable A is selected by its tag "tagA".
        # Expected: A is free, then fixed, then free again.
        "tag",
        # C2: Variable A is selected by its name "A".
        # Expected: A is free, then fixed, then free again.
        "name",
        # C3: Variable A is selected by passing the variable itself.
        # Expected: A is free, then fixed, then free again.
        "variable",
    ],
)
def test_fix_and_free_one_variable(
    sine_recipe_with_variables, input_selection
):
    recipe = sine_recipe_with_variables
    identifiers = {"tag": "tagA", "name": "A", "variable": recipe.A}
    identifier = identifiers[input_selection]
    actual_is_free_states = [recipe.is_free(recipe.A)]
    recipe.fix(identifier)
    actual_is_free_states.append(recipe.is_free(recipe.A))
    recipe.free(identifier)
    actual_is_free_states.append(recipe.is_free(recipe.A))
    expected_is_free_states = [True, False, True]
    assert actual_is_free_states == expected_is_free_states


# Fix and then free several variables at once. The recipe has the
# variables A, k, c and B, which are all free to begin with.
@pytest.mark.parametrize(
    "input_fix_args, input_fix_kwargs, input_free_args, "
    "expected_free_variables, expected_c_value",
    [
        # C1: A is fixed, then "all" variables are freed.
        # Expected: All variables are free and c keeps its value of 0.
        (["A"], {}, ["all"], ["A", "k", "c", "B"], 0),
        # C2: A is fixed by name, k by tag and c by name as a keyword with
        # the value 3.
        # Expected: Only B is free and c takes the value 3.
        (["A", "tagk"], {"c": 3}, [], ["B"], 3),
        # C3: "all" variables are fixed.
        # Expected: No variable is free and c keeps its value of 0.
        (["all"], {}, [], [], 0),
    ],
)
def test_fix_and_free_several_variables(
    sine_recipe_with_variables,
    input_fix_args,
    input_fix_kwargs,
    input_free_args,
    expected_free_variables,
    expected_c_value,
):
    recipe = sine_recipe_with_variables
    recipe.fix(*input_fix_args, **input_fix_kwargs)
    recipe.free(*input_free_args)
    actual_free_variables = [
        name
        for name in ["A", "k", "c", "B"]
        if recipe.is_free(getattr(recipe, name))
    ]
    actual_c_value = recipe.c.value
    assert actual_free_variables == expected_free_variables
    assert actual_c_value == expected_c_value


# Fix or free with an identifier that does not select any variable.
@pytest.mark.parametrize(
    "input_method_name, input_args, input_kwargs",
    [
        # C1: An unknown name is freed.
        # Expected: Raises a ValueError.
        ("free", ("junk",), {}),
        # C2: An unknown name is fixed.
        # Expected: Raises a ValueError.
        ("fix", ("junk",), {}),
        # C3: A tag is fixed as if it were a variable name in a keyword.
        # Expected: Raises a ValueError.
        ("fix", (), {"tagA": 1}),
    ],
)
def test_fix_and_free_bad(
    sine_recipe_with_variables, input_method_name, input_args, input_kwargs
):
    recipe = sine_recipe_with_variables
    with pytest.raises(ValueError):
        getattr(recipe, input_method_name)(*input_args, **input_kwargs)


def test_get_names_and_values(sine_recipe_with_variables):
    # Get the names and values of the free variables after each change to
    # the variables A=2, k=1, c=0 and B=0, which were added in that order.
    # Each step changes the recipe left by the step before it.
    recipe = sine_recipe_with_variables
    actual_names_and_values = []
    # C1: The variables are left as they were added.
    # Expected: All variables are listed in the order they were added.
    actual_names_and_values.append(
        (recipe.get_names(), recipe.get_values().tolist())
    )
    # C2: B is constrained to a parameter.
    # Expected: B is fixed by the constraint and no longer listed.
    recipe.add_constraint(recipe.B, Parameter("Bpar", -1))
    actual_names_and_values.append(
        (recipe.get_names(), recipe.get_values().tolist())
    )
    # C3: B is deleted and k is fixed.
    # Expected: Only A and c are listed.
    recipe.delete_variable(recipe.B)
    recipe.fix(recipe.k)
    actual_names_and_values.append(
        (recipe.get_names(), recipe.get_values().tolist())
    )
    # C4: All variables are fixed.
    # Expected: Nothing is listed.
    recipe.fix("all")
    actual_names_and_values.append(
        (recipe.get_names(), recipe.get_values().tolist())
    )
    # C5: All variables are freed.
    # Expected: The remaining variables A, k and c are all listed.
    recipe.free("all")
    actual_names_and_values.append(
        (recipe.get_names(), recipe.get_values().tolist())
    )
    expected_names_and_values = [
        (["A", "k", "c", "B"], [2, 1, 0, 0]),
        (["A", "k", "c"], [2, 1, 0]),
        (["A", "c"], [2, 0]),
        ([], []),
        (["A", "k", "c"], [2, 1, 0]),
    ]
    assert actual_names_and_values == expected_names_and_values


def test_add_constraint_to_constrained_parameter(build_sine_recipe):
    # C1: c is constrained to 2*A, then A is constrained to the new
    # variable Avar, which has no value yet.
    # Expected: Avar takes the value of A, which is 1, and c follows the
    # chain of constraints to 2.
    recipe = build_sine_recipe()
    variable = recipe.create_new_variable("Avar")
    recipe.add_constraint(recipe.cont.c, "2*A", {"A": recipe.cont.A})
    recipe.add_constraint(recipe.cont.A, variable)
    actual_values = [
        variable.get_value(),
        recipe.cont.A.get_value(),
        recipe.cont.c.get_value(),
    ]
    expected_values = [1, 1, 2]
    assert actual_values == expected_values


# The residual of the sine recipe is (ycalc - y)/dy, followed by the square
# root of each restraint penalty. The recipe starts with ycalc = sin(x) = y
# and its profile has no uncertainties, so dy = 1.


# Set the phase c of the sine recipe, so that ycalc = sin(x + c).
@pytest.mark.parametrize(
    "input_c_value",
    [
        # C1: c is left at 0.
        # Expected: The residual is sin(x) - y, which is zero.
        0,
        # C2: c is set to 1.
        # Expected: The residual is sin(x + 1) - y.
        1,
    ],
)
def test_residual_with_parameter_value(build_sine_recipe, input_c_value):
    recipe = build_sine_recipe()
    profile = recipe.cont.profile
    recipe.cont.c.set_value(input_c_value)
    actual_residual = recipe.residual()
    expected_residual = sin(profile.x + input_c_value) - profile.y
    assert actual_residual == pytest.approx(expected_residual)


# Replace the equation of an evaluated sine recipe.
@pytest.mark.parametrize(
    "input_equation, expected_residual_value",
    [
        # C1: The equation is set to the observed profile "y".
        # Expected: The residual is zero at every point.
        ("y", 0),
        # C2: The equation is set to "y+dy".
        # Expected: The residual is dy/dy = 1 at every point.
        ("y+dy", 1),
    ],
)
def test_residual_with_new_equation(
    build_sine_recipe, input_equation, expected_residual_value
):
    recipe = build_sine_recipe()
    profile = recipe.cont.profile
    recipe.residual()
    recipe.cont.set_equation(input_equation)
    actual_residual = recipe.residual()
    expected_residual = np.full_like(profile.y, expected_residual_value)
    assert actual_residual == pytest.approx(expected_residual)


def test_residual_with_recipe_constraints_and_restraints(build_sine_recipe):
    # Constrain and restrain c through the recipe, then undo each change.
    # Each step changes the recipe left by the step before it, and the
    # residual is evaluated after each step.
    recipe = build_sine_recipe()
    x, y = recipe.cont.profile.x, recipe.cont.profile.y
    # C1: c is constrained to 2*A, and A to the new variable Avar.
    # Expected: c is 2, so the residual is sin(x + 2) - y.
    variable = recipe.create_new_variable("Avar")
    recipe.add_constraint(recipe.cont.c, "2*A", {"A": recipe.cont.A})
    recipe.add_constraint(recipe.cont.A, variable)
    actual_constrained_residual = recipe.residual()
    expected_constrained_residual = sin(x + 2) - y
    assert actual_constrained_residual == pytest.approx(
        expected_constrained_residual
    )
    # C2: c is also restrained to zero with sig=1.
    # Expected: As C1, followed by a penalty term of |c - 0| = 2.
    restraint = recipe.add_soft_bounds(recipe.cont.c, 0, 0, 1)
    actual_restrained_residual = recipe.residual()
    expected_restrained_residual = np.append(sin(x + 2) - y, 2)
    assert actual_restrained_residual == pytest.approx(
        expected_restrained_residual
    )
    # C3: The constraint on c is removed and c is reset to 0.
    # Expected: The residual is zero, followed by a penalty term of 0 for
    # the restraint, which is now met.
    recipe.remove_constraint(recipe.cont.c)
    recipe.cont.c.set_value(0)
    actual_unconstrained_residual = recipe.residual()
    expected_unconstrained_residual = np.append(sin(x) - y, 0)
    assert actual_unconstrained_residual == pytest.approx(
        expected_unconstrained_residual
    )
    # C4: The restraint and the variable Avar are removed.
    # Expected: The residual is zero, with no penalty term.
    recipe.remove_soft_bounds(restraint)
    recipe.delete_variable(recipe.Avar)
    actual_final_residual = recipe.residual()
    expected_final_residual = sin(x) - y
    assert actual_final_residual == pytest.approx(expected_final_residual)


def test_residual_with_contribution_constraints_and_restraints(
    build_sine_recipe,
):
    # Constrain and restrain c through the contribution of an evaluated
    # recipe, then undo each change. Each step changes the recipe left by
    # the step before it, and the residual is evaluated after each step.
    # A change made on a contribution does not tell the recipe that its
    # configuration changed, so the recipe is reset by hand after adding or
    # removing a restraint.
    recipe = build_sine_recipe()
    x, y = recipe.cont.profile.x, recipe.cont.profile.y
    recipe.residual()
    # C1: The contribution constrains c to 2*A.
    # Expected: c is 2, so the residual is sin(x + 2) - y.
    recipe.cont.add_constraint(recipe.cont.c, "2*A")
    actual_constrained_residual = recipe.residual()
    expected_constrained_residual = sin(x + 2) - y
    assert actual_constrained_residual == pytest.approx(
        expected_constrained_residual
    )
    # C2: The contribution also restrains c to zero with sig=1.
    # Expected: As C1, followed by a penalty term of |c - 0| = 2.
    restraint = recipe.cont.add_soft_bounds(recipe.cont.c, 0, 0, 1)
    recipe._ready = False
    actual_restrained_residual = recipe.residual()
    expected_restrained_residual = np.append(sin(x + 2) - y, 2)
    assert actual_restrained_residual == pytest.approx(
        expected_restrained_residual
    )
    # C3: The restraint and the constraint on c are removed and c is reset
    # to 0.
    # Expected: The residual is zero, with no penalty term.
    recipe.cont.remove_soft_bounds(restraint)
    recipe._ready = False
    recipe.cont.remove_constraint(recipe.cont.c)
    recipe.cont.c.set_value(0)
    actual_final_residual = recipe.residual()
    expected_final_residual = sin(x) - y
    assert actual_final_residual == pytest.approx(expected_final_residual)


def _get_restraint_bounds(restraints):
    """Return the sorted (lower_bound, upper_bound, sig) of each
    restraint."""
    sorted_restraints = sorted(
        (restraint.lower_bound, restraint.upper_bound, restraint.sig)
        for restraint in restraints
    )
    return sorted_restraints


def test_add_bounds_after_recipe_evaluation():
    # Test that soft bounds added to a contribution after the recipe is
    # evaluated are used in the next evaluation.
    profile = Profile()
    profile.set_observed_profile(linspace(0, 1, 5), ones(5))
    contribution = FitContribution("contrib")
    contribution.set_profile(profile)
    contribution.set_equation("m*x + b")
    # C1: m is restrained in the contribution before it is added to a recipe.
    # Expected: The contribution holds the restraint on m.
    contribution.add_soft_bounds("m", lower_bound=0, upper_bound=2, sig=1)
    actual_bounds_before_eval = _get_restraint_bounds(
        contribution._get_restraints()
    )
    expected_bounds_before_eval = [(0, 2, 1)]
    assert actual_bounds_before_eval == expected_bounds_before_eval
    # C2: The recipe is evaluated, then b is restrained in the contribution
    # and the recipe is evaluated again.
    # Expected: The recipe uses the restraints on both m and b.
    recipe = FitRecipe()
    recipe.add_contribution(contribution)
    recipe.add_variable(contribution.m, 1)
    recipe.add_variable(contribution.b, 0)
    # Evaluate the recipe.
    recipe.residual()
    recipe.contrib.add_soft_bounds("b", lower_bound=-1, upper_bound=1, sig=3)
    # Evaluate the recipe again to include the new restraint.
    recipe.residual()
    actual_bounds_after_eval = _get_restraint_bounds(recipe._restraintlist)
    expected_bounds_after_eval = [(-1, 1, 3), (0, 2, 1)]
    assert actual_bounds_after_eval == expected_bounds_after_eval


def test_add_bounds_to_nested_parameter_set_after_recipe_evaluation(
    build_sine_recipe,
):
    # Test that soft bounds added to a ParameterSet nested two levels below
    # the recipe after the recipe is evaluated are used in the next
    # evaluation.
    outer = ParameterSet("outer")
    inner = ParameterSet("inner")
    inner.new_parameter("p", 1)
    inner.new_parameter("q", 1)
    # C1: p is restrained in "inner" before "inner" is added to a recipe.
    # Expected: "inner" holds the restraint on p.
    inner.add_soft_bounds("p", lower_bound=0, upper_bound=2, sig=1)
    actual_bounds_before_eval = _get_restraint_bounds(inner._get_restraints())
    expected_bounds_before_eval = [(0, 2, 1)]
    assert actual_bounds_before_eval == expected_bounds_before_eval
    # C2: The recipe holding "outer", which holds "inner", is evaluated, then
    # q is restrained in "inner" and the recipe is evaluated again.
    # Expected: The recipe uses the restraints on both p and q.
    recipe = build_sine_recipe()
    outer.add_parameter_set(inner)
    recipe.add_parameter_set(outer)
    # Evaluate the recipe.
    recipe.residual()
    recipe.outer.inner.add_soft_bounds(
        "q", lower_bound=-1, upper_bound=1, sig=3
    )
    # Evaluate the recipe again to include the new restraint.
    recipe.residual()
    actual_bounds_after_eval = _get_restraint_bounds(recipe._restraintlist)
    expected_bounds_after_eval = [(-1, 1, 3), (0, 2, 1)]
    assert actual_bounds_after_eval == expected_bounds_after_eval


def test_add_constraint_after_recipe_evaluation(build_sine_recipe):
    # Test that a constraint added to a contribution after the recipe is
    # evaluated is updated in the next evaluations.
    recipe = build_sine_recipe()
    recipe.add_variable(recipe.cont.A)
    # Evaluate the recipe.
    recipe.residual()
    # C1: c is constrained to 2*A in the contribution with A = 1.
    # Expected: c is 2.
    recipe.cont.add_constraint(recipe.cont.c, "2*A")
    actual_c_value = recipe.cont.c.value
    expected_c_value = 2
    assert actual_c_value == expected_c_value
    # C2: The recipe is evaluated again with A = 3.
    # Expected: c follows A and is 6.
    recipe.residual([3])
    actual_c_value = recipe.cont.c.value
    expected_c_value = 6
    assert actual_c_value == expected_c_value


def test_convert_bounds_to_restraints():
    # C1: A variable bounded to [-1, 1] is converted with sig=2 and
    # scaled=True.
    # Expected: The recipe holds one restraint with those bounds, sig and
    # scaling.
    recipe = FitRecipe("recipe")
    recipe.create_new_variable("var1", 1)
    recipe.var1.bounds = (-1, 1)
    recipe.convert_bounds_to_restraints(sig=2, scaled=True)
    actual_restraints = [
        (
            restraint.lower_bound,
            restraint.upper_bound,
            restraint.sig,
            restraint.scaled,
        )
        for restraint in recipe._restraints
    ]
    expected_restraints = [(-1, 1, 2, True)]
    assert actual_restraints == expected_restraints


def test_add_and_remove_parameter_set():
    # C1: A parameter is added to the recipe and then removed.
    # Expected: After adding, the recipe holds the parameter under its
    # name. After removing, the name is gone.
    recipe = FitRecipe("recipe")
    parameter_to_add = Parameter("added_param", 1)
    recipe.add_parameter_set(parameter_to_add)
    actual_added_parameter = recipe.added_param
    actual_added_value = recipe.added_param.value
    expected_added_parameter = parameter_to_add
    expected_added_value = 1
    assert actual_added_parameter == expected_added_parameter
    assert actual_added_value == expected_added_value
    recipe.remove_parameter_set(parameter_to_add)
    actual_has_parameter = hasattr(recipe, "added_param")
    expected_has_parameter = False
    assert actual_has_parameter == expected_has_parameter


def _get_printed_line_labels(output):
    """Return the label that starts each printed line.

    A line that holds only the iteration count is labelled "count".
    """
    line_labels = []
    for line in output.splitlines():
        if line.strip().isdigit():
            line_labels.append("count")
        else:
            line_labels.append(line.split(":")[0].split("=")[0].strip())
    return line_labels


# Print the progress of a refinement with the default PrintFitHook. The
# sine recipe has the variable c and the restraint c >= 5, and each case
# evaluates the scalar residual once.
@pytest.mark.parametrize(
    "input_verbose, expected_line_labels",
    [
        # C1: The verbosity is 0.
        # Expected: Nothing is printed.
        (0, []),
        # C2: The verbosity is 1.
        # Expected: Only the iteration count is printed.
        (1, ["count"]),
        # C3: The verbosity is 2.
        # Expected: The count is followed by the residual, the
        # contributions and the restraints.
        (2, ["count", "Residual", "FitContributions", "Restraints"]),
        # C4: The verbosity is 3.
        # Expected: As C3, followed by the variables and their values.
        (
            3,
            [
                "count",
                "Residual",
                "FitContributions",
                "Restraints",
                "Variables",
                "c",
            ],
        ),
    ],
)
def test_print_fit_hook_output(
    build_sine_recipe, capturestdout, input_verbose, expected_line_labels
):
    recipe = build_sine_recipe()
    recipe.add_variable(recipe.cont.c)
    recipe.add_soft_bounds("c", lower_bound=5)
    (print_fit_hook,) = recipe.get_fit_hooks()
    print_fit_hook.verbose = input_verbose
    actual_output = capturestdout(recipe.scalar_residual)
    actual_line_labels = _get_printed_line_labels(actual_output)
    assert actual_line_labels == expected_line_labels


def optimize_recipe(recipe):
    recipe.fithooks[0].verbose = 0
    residuals = recipe.residual
    values = recipe.values
    leastsq(residuals, values)


def test_initialize_recipe_from_recipe(build_recipes_one_contribution):
    # C1: A FitRecipe is initialized from a previously optimized recipe.
    # Expected: The recipe gets everything: contributions, profiles
    # (contained in contributions), variables, restraints and
    # constraints, and refines to the same variable names and values.
    recipe1, _ = build_recipes_one_contribution
    optimize_recipe(recipe1)
    expected_parameters_dict = recipe1._parameters
    expected_constraints_dict = recipe1._constraints
    expected_restraints_set = recipe1._restraints
    expected_contributions_dict = recipe1._contributions
    expected_profiles_list = []
    for contribution in expected_contributions_dict.values():
        expected_profiles_list.append(contribution.profile)

    recipe2 = FitRecipe()
    actual_recipes_differ = recipe1 != recipe2
    expected_recipes_differ = True
    assert actual_recipes_differ == expected_recipes_differ
    recipe2.initialize_recipe_with_recipe(recipe1)
    actual_parameters_dict = recipe2._parameters
    actual_constraints_dict = recipe2._constraints
    actual_restraints_set = recipe2._restraints
    actual_contributions_dict = recipe2._contributions
    actual_profiles_list = []
    for contribution in actual_contributions_dict.values():
        actual_profiles_list.append(contribution.profile)

    assert actual_parameters_dict == expected_parameters_dict
    assert actual_constraints_dict == expected_constraints_dict
    assert actual_restraints_set == expected_restraints_set
    assert actual_contributions_dict == expected_contributions_dict
    assert actual_profiles_list == expected_profiles_list

    # Check that the refined values and variable names are the same in
    # the results objects for each recipe.
    results1 = FitResults(recipe1)
    # Round to account for small numerical differences.
    expected_values = sorted(np.round(results1.varvals, 7).tolist())
    expected_names = sorted(results1.varnames)

    optimize_recipe(recipe2)
    results2 = FitResults(recipe2)
    # Round to account for small numerical differences.
    actual_values = sorted(np.round(results2.varvals, 7).tolist())
    actual_names = sorted(results2.varnames)

    assert actual_names == expected_names
    assert actual_values == expected_values


def test_initialize_recipe_from_recipe_bad():
    # C1: A FitRecipe is initialized from an integer instead of a recipe.
    # Expected: Raises a ValueError that names the type it got.
    recipe_bad = 12345
    recipe2 = FitRecipe()
    expected_error_message = (
        "The input recipe_object must be a FitRecipe, "
        "but got <class 'int'>."
    )
    with pytest.raises(ValueError, match=re.escape(expected_error_message)):
        recipe2.initialize_recipe_with_recipe(recipe_bad)


def test_initialize_recipe_from_results_object(build_recipes_one_contribution):
    # C1: A FitRecipe with an extra variable is initialized from the
    # FitResults of a previous fit.
    # Expected: The recipe refines to the variables of the previous fit,
    # and keeps its extra variable.
    recipe1, recipe2 = build_recipes_one_contribution
    optimize_recipe(recipe1)
    results1 = FitResults(recipe1)
    expected_values = np.round(results1.varvals, 5)
    expected_names = results1.varnames

    actual_recipes_differ = recipe1 != recipe2
    expected_recipes_differ = True
    assert actual_recipes_differ == expected_recipes_differ
    # Create a new variable that should be kept in the initialized recipe.
    recipe2.create_new_variable("extra_var", 5)
    actual_values_before_init = [val for val in recipe2.get_values()]
    actual_names_before_init = sorted(recipe2.get_names())
    # The three variables plus extra_var.
    expected_names_before_init = sorted(
        [
            "amplitude",
            "extra_var",
            "phase_shift",
            "wave_number",
        ]
    )
    expected_values_before_init = [4, 3, 2, 5]
    assert actual_values_before_init == expected_values_before_init
    assert actual_names_before_init == expected_names_before_init
    recipe2.initialize_recipe_with_results(results1)
    optimize_recipe(recipe2)
    results2 = FitResults(recipe2)
    actual_values = sorted(np.round(results2.varvals, 5).tolist())
    actual_names = sorted(results2.varnames)

    # Add the name and value of the new variable to the expected ones.
    expected_names = sorted(expected_names + ["extra_var"])
    expected_values = sorted(expected_values.tolist() + [5])
    assert actual_names == expected_names
    assert actual_values == expected_values


def test_initialize_recipe_from_results_file(
    build_recipes_one_contribution, temp_data_files
):
    # C1: A FitRecipe is initialized from a FitResults file.
    # Expected: The recipe takes the variable values saved in the file.
    results_file = temp_data_files / "fit_results.res"
    expected_names = ["amplitude", "phase_shift", "wave_number"]
    expected_values = [1, 1, 0]

    recipe, _ = build_recipes_one_contribution
    recipe.initialize_recipe_with_results(results_file)
    results = FitResults(recipe)
    actual_values = np.round(results.varvals, 5).tolist()
    actual_names = sorted(results.varnames)

    assert actual_names == expected_names
    assert actual_values == expected_values


def test_initialize_recipe_from_results_file_bad(
    build_recipes_one_contribution,
):
    # C1: A recipe is initialized from an integer, which is not a path,
    # str or FitResults object.
    # Expected: Raises a ValueError that names the type it got.
    recipe, _ = build_recipes_one_contribution
    bad_input = 12345
    expected_error_message = (
        "The input results must be a FitResults object or a path to a "
        "results file, but got <class 'int'>."
    )
    with pytest.raises(ValueError, match=re.escape(expected_error_message)):
        recipe.initialize_recipe_with_results(bad_input)


def test_initialize_recipe_from_results_file_wrong(
    build_recipe_two_contributions, temp_data_files, capsys
):
    # C1: A FitRecipe is initialized from a results file whose
    # parameters are not in the recipe.
    # Expected: A warning is printed for each unknown parameter and the
    # recipe keeps its own variables.
    results_file_from_single_contrib = temp_data_files / "fit_results.res"
    recipe = build_recipe_two_contributions
    recipe.initialize_recipe_with_results(results_file_from_single_contrib)
    captured = capsys.readouterr()
    actual_print_msg = captured.out

    results_file_param_names = ["amplitude", "phase_shift", "wave_number"]
    expected_print_messages = []
    for param_name in results_file_param_names:
        msg = (
            f"Warning: Parameter '{param_name}' from results not found "
            "in FitRecipe and will be ignored."
        )
        expected_print_messages.append(msg)

    for expected_print_msg in expected_print_messages:
        assert expected_print_msg in actual_print_msg


def get_labels_and_linecount(ax):
    """Helper to get line labels and count from a matplotlib Axes."""
    labels = [
        line.get_label()
        for line in ax.get_lines()
        if not line.get_label().startswith("_")
    ]
    line_count = len(
        [
            line
            for line in ax.get_lines()
            if not line.get_label().startswith("_")
        ]
    )
    return labels, line_count


def build_recipe_from_datafile(datafile):
    """Helper to build a FitRecipe from a datafile using
    ProfileParser."""
    profile = Profile()
    parser = ProfileParser()
    parser.parse_file(str(datafile))
    profile.load_parsed_data(parser)

    contribution = FitContribution("c")
    contribution.set_profile(profile)
    contribution.set_equation("m*x + b")
    recipe = FitRecipe()
    recipe.add_contribution(contribution)
    recipe.add_variable(contribution.m, 1)
    recipe.add_variable(contribution.b, 0)
    return recipe


def test_plot_recipe_bad_display(build_recipes_one_contribution):
    recipe, _ = build_recipes_one_contribution
    # C1: The observed, fit and difference curves are all disabled.
    # Expected: Raises a ValueError with a message.
    plt.close("all")
    msg = "At least one of show_observed, show_fit, or show_diff must be True"
    with pytest.raises(ValueError, match=msg):
        recipe.plot_recipe(
            show_observed=False,
            show_diff=False,
            show_fit=False,
        )


def test_plot_recipe_no_contribution():
    recipe = FitRecipe()
    # C1: The recipe has no contributions.
    # Expected: Raises a ValueError with a message.
    plt.close("all")
    msg = (
        "No contributions found in recipe. "
        "Add contributions before plotting."
    )
    with pytest.raises(ValueError, match=msg):
        recipe.plot_recipe()


def test_plot_recipe_before_refinement(capsys, build_recipes_one_contribution):
    # C1: The recipe is plotted before refinement.
    # Expected: The observed data is plotted without the fit or
    # difference curves, and a warning is printed.
    recipe, _ = build_recipes_one_contribution
    plt.close("all")
    before = set(plt.get_fignums())
    # Include fit_label="nothing" to make sure the fit is not plotted.
    fig, ax = recipe.plot_recipe(
        show=False, data_label="my data", fit_label="nothing", return_fig=True
    )
    after = set(plt.get_fignums())
    new_figs = after - before
    captured = capsys.readouterr()
    actual_message = captured.out.strip()
    expected_message = (
        "Contribution 'c1' has no calculated values (ycalc is None). "
        "Only observed data will be plotted."
    )
    # Get the labels from the plotted line.
    actual_label, actual_line_count = get_labels_and_linecount(ax)
    expected_line_count = 1
    expected_label = ["my data"]
    assert actual_line_count == expected_line_count
    assert actual_label == expected_label
    actual_new_figure_count = len(new_figs)
    expected_new_figure_count = 1
    assert actual_new_figure_count == expected_new_figure_count
    assert actual_message == expected_message


def test_plot_recipe_after_refinement(build_recipes_one_contribution):
    # C1: The recipe is refined and then plotted.
    # Expected: One figure with the observed, calculated and difference
    # curves is created.
    recipe, _ = build_recipes_one_contribution
    optimize_recipe(recipe)
    plt.close("all")
    before = set(plt.get_fignums())
    fig, ax = recipe.plot_recipe(show=False, return_fig=True)
    after = set(plt.get_fignums())
    new_figs = after - before
    actual_label, actual_line_count = get_labels_and_linecount(ax)
    expected_label = ["Observed", "Calculated", "Difference"]
    expected_line_count = 3
    assert actual_line_count == expected_line_count
    assert actual_label == expected_label
    actual_new_figure_count = len(new_figs)
    expected_new_figure_count = 1
    assert actual_new_figure_count == expected_new_figure_count


def test_plot_recipe_styles_only_its_own_figure(
    build_recipes_one_contribution,
):
    # C1: A refined recipe is plotted.
    # Expected: The srfit figure uses the bg-style stylesheet, and the
    # global matplotlib style is the same after plotting as before.
    recipe, _ = build_recipes_one_contribution
    optimize_recipe(recipe)
    plt.close("all")
    expected_global_linewidth = plt.rcParams["axes.linewidth"]
    _, ax = recipe.plot_recipe(show=False, return_fig=True)
    actual_figure_linewidth = ax.spines["left"].get_linewidth()
    expected_figure_linewidth = all_styles["bg-style"]["axes.linewidth"]
    actual_global_linewidth = plt.rcParams["axes.linewidth"]
    assert actual_figure_linewidth == expected_figure_linewidth
    assert actual_global_linewidth == expected_global_linewidth


def test_plot_recipe_two_contributions(build_recipe_two_contributions):
    # C1: The recipe has two contributions.
    # Expected: Two figures are created, each with the observed,
    # calculated and difference curves.
    recipe = build_recipe_two_contributions
    optimize_recipe(recipe)
    plt.close("all")
    before = set(plt.get_fignums())
    figs, axes = recipe.plot_recipe(show=False, return_fig=True)
    for ax in axes:
        actual_label, actual_line_count = get_labels_and_linecount(ax)
        expected_label = ["Observed", "Calculated", "Difference"]
        expected_line_count = 3
        assert actual_line_count == expected_line_count
        assert actual_label == expected_label
    after = set(plt.get_fignums())
    new_figs = after - before
    actual_new_figure_count = len(new_figs)
    expected_new_figure_count = 2
    assert actual_new_figure_count == expected_new_figure_count


def test_plot_recipe_on_existing_plot(build_recipes_one_contribution):
    # C1: The user passes existing axes, with a title and a curve, to
    # plot_recipe.
    # Expected: The user's title and curve are kept in the final figure.
    recipe, _ = build_recipes_one_contribution
    optimize_recipe(recipe)
    plt.close("all")
    fig, ax = plt.subplots()
    ax.set_title("User Title")
    ax.plot([0, 1], [0, 1], label="New Data")
    recipe.plot_recipe(ax=ax, show=False)
    actual_title = ax.get_title()
    expected_title = "User Title"
    actual_labels, actual_line_count = get_labels_and_linecount(ax)
    expected_line_count = 4
    expected_labels = ["New Data", "Observed", "Calculated", "Difference"]
    assert actual_line_count == expected_line_count
    assert actual_labels == expected_labels
    assert actual_title == expected_title


def test_plot_recipe_add_new_data(build_recipes_one_contribution):
    # C1: The user adds a curve to the figure returned by plot_recipe.
    # Expected: The new curve is added to the existing figure.
    recipe, _ = build_recipes_one_contribution
    optimize_recipe(recipe)
    plt.close("all")
    before = set(plt.get_fignums())
    fig, ax = recipe.plot_recipe(return_fig=True, show=False)
    after = set(plt.get_fignums())
    new_figs = after - before
    # Add new data to the existing plot.
    ax.plot([0, pi], [0, 0], label="New Data")
    ax.legend()
    actual_labels, actual_line_count = get_labels_and_linecount(ax)
    expected_labels = ["Observed", "Calculated", "Difference", "New Data"]
    expected_line_count = 4
    actual_new_figure_count = len(new_figs)
    expected_new_figure_count = 1
    assert actual_new_figure_count == expected_new_figure_count
    assert actual_line_count == expected_line_count
    assert actual_labels == expected_labels


def test_plot_recipe_add_new_data_two_figs(build_recipe_two_contributions):
    # C1: The user adds a curve to each figure returned by plot_recipe
    # for a recipe with two contributions.
    # Expected: The new curve is added to each existing figure.
    recipe = build_recipe_two_contributions
    optimize_recipe(recipe)
    plt.close("all")
    before = set(plt.get_fignums())
    figure, axes = recipe.plot_recipe(return_fig=True, show=False)
    after = set(plt.get_fignums())
    new_figs = after - before
    # Add new data to the existing plots.
    for ax in axes:
        ax.plot([0, pi], [0, 0], label="New Data")
        ax.legend()
        actual_labels, actual_line_count = get_labels_and_linecount(ax)
        expected_labels = ["Observed", "Calculated", "Difference", "New Data"]
        expected_line_count = 4
        assert actual_line_count == expected_line_count
        assert actual_labels == expected_labels
    actual_new_figure_count = len(new_figs)
    expected_new_figure_count = 2
    assert actual_new_figure_count == expected_new_figure_count


def test_plot_recipe_set_title(build_recipes_one_contribution):
    # C1: The user sets the title through plot_recipe.
    # Expected: The figure has that title.
    recipe, _ = build_recipes_one_contribution
    optimize_recipe(recipe)
    plt.close("all")
    expected_title = "Custom Recipe Title"
    figure, ax = recipe.plot_recipe(
        title=expected_title, return_fig=True, show=False
    )
    actual_title = ax.get_title()
    assert actual_title == expected_title


def test_plot_recipe_default_title(build_recipes_one_contribution):
    # C1: A single contribution is plotted with no title given.
    # Expected: The figure is titled with the name of the contribution.
    recipe, _ = build_recipes_one_contribution
    optimize_recipe(recipe)
    plt.close("all")
    _, ax = recipe.plot_recipe(show=False, return_fig=True)
    actual_title = ax.get_title()
    expected_title = "c1"
    assert actual_title == expected_title


def test_plot_recipe_titles_two_contributions(build_recipe_two_contributions):
    # C1: Two contributions are plotted on separate figures with no
    # title given.
    # Expected: Each figure is titled with the name of its contribution.
    recipe = build_recipe_two_contributions
    optimize_recipe(recipe)
    plt.close("all")
    _, axes = recipe.plot_recipe(show=False, return_fig=True)
    actual_titles = [ax.get_title() for ax in axes]
    expected_titles = ["c1", "c2"]
    assert actual_titles == expected_titles


def test_plot_recipe_labels_shared_axes(build_recipe_two_contributions):
    # C1: Two contributions are plotted on a single user-supplied axes.
    # Expected: Legend labels are prefixed with the contribution name so
    # that the curves of each contribution can be told apart.
    recipe = build_recipe_two_contributions
    optimize_recipe(recipe)
    plt.close("all")
    _, ax = plt.subplots()
    recipe.plot_recipe(ax=ax, show=False)
    actual_labels, _ = get_labels_and_linecount(ax)
    expected_labels = [
        "c1: Observed",
        "c1: Calculated",
        "c1: Difference",
        "c2: Observed",
        "c2: Calculated",
        "c2: Difference",
    ]
    assert actual_labels == expected_labels


def test_plot_recipe_label_contribution_placeholder(
    build_recipe_two_contributions,
):
    # C1: User passes a label containing the {contribution} placeholder.
    # Expected: The placeholder is replaced by the contribution name and no
    # prefix is added.
    recipe = build_recipe_two_contributions
    optimize_recipe(recipe)
    plt.close("all")
    _, axes = recipe.plot_recipe(
        fit_label="{contribution} calculated",
        show_observed=False,
        show_diff=False,
        show=False,
        return_fig=True,
    )
    actual_labels = [get_labels_and_linecount(ax)[0] for ax in axes]
    expected_labels = [["c1 calculated"], ["c2 calculated"]]
    assert actual_labels == expected_labels


def test_plot_recipe_missing_ycalc_one_contribution(
    build_recipe_two_contributions,
):
    # C1: Only the second contribution has been evaluated, so the ycalc of
    # the first one is None.
    # Expected: The unevaluated contribution shows observed data only while
    # the other one is still plotted in full.
    recipe = build_recipe_two_contributions
    optimize_recipe(recipe)
    plt.close("all")
    # Manually set ycalc to None.
    recipe._contributions["c1"].profile.ycalc = None
    _, axes = recipe.plot_recipe(show=False, return_fig=True)
    actual_labels = [get_labels_and_linecount(ax)[0] for ax in axes]
    expected_labels = [
        ["Observed"],
        ["Observed", "Calculated", "Difference"],
    ]
    assert actual_labels == expected_labels


def test_plot_recipe_set_defaults(build_recipes_one_contribution):
    # C1: The user sets default plot options with set_plot_defaults.
    # Expected: plot_recipe uses the default options when called without
    # arguments.
    recipe, _ = build_recipes_one_contribution
    optimize_recipe(recipe)
    plt.close("all")
    # Set new defaults.
    recipe.set_plot_defaults(
        show_observed=False,
        show_fit=True,
        show_diff=False,
        data_label="Data Default",
        fit_label="Fit Default",
        diff_label="Diff Default",
        title="Default Title",
    )
    # Call plot_recipe without any arguments.
    figure, ax = recipe.plot_recipe(return_fig=True, show=False)
    actual_title = ax.get_title()
    actual_labels, actual_line_count = get_labels_and_linecount(ax)
    expected_title = "Default Title"
    expected_labels = ["Fit Default"]
    expected_line_count = 1
    assert actual_title == expected_title
    assert actual_line_count == expected_line_count
    assert actual_labels == expected_labels


def test_plot_recipe_set_defaults_bad(capsys, build_recipes_one_contribution):
    # C1: The user sets a default that is not a plot_recipe option.
    # Expected: A warning is printed and the recipe is still plotted.
    recipe, _ = build_recipes_one_contribution
    optimize_recipe(recipe)
    plt.close("all")
    recipe.set_plot_defaults(
        invalid_option="blah",
    )
    captured = capsys.readouterr()
    actual_msg = captured.out.strip()
    expected_msg = (
        "Warning: 'invalid_option' is not a valid "
        "plot_recipe option and will be ignored."
    )
    assert actual_msg == expected_msg
    before = set(plt.get_fignums())
    figure, ax = recipe.plot_recipe(return_fig=True, show=False)
    after = set(plt.get_fignums())
    new_figs = after - before
    actual_new_figure_count = len(new_figs)
    expected_new_figure_count = 1
    assert actual_new_figure_count == expected_new_figure_count


# Infer the axis labels of the plot from the extension of the data file.
@pytest.mark.parametrize(
    "input_filename, expected_xlabel, expected_ylabel",
    [
        # C1: A .gr file.
        # Expected: The labels of a PDF are used.
        ("gr_file.gr", r"r ($\mathrm{\AA}$)", r"G ($\mathrm{\AA}^{-2}$)"),
        # C2: A .cgr file.
        # Expected: The labels of a PDF are used.
        ("cgr_file.cgr", r"r ($\mathrm{\AA}$)", r"G ($\mathrm{\AA}^{-2}$)"),
        # C3: A .dat file, whose kind of data is unknown.
        # Expected: The labels are empty.
        ("dat_file.dat", "", ""),
    ],
)
def test_plot_recipe_labels_from_gr_file(
    temp_data_files, input_filename, expected_xlabel, expected_ylabel
):
    data_file = temp_data_files / input_filename
    recipe = build_recipe_from_datafile(data_file)
    optimize_recipe(recipe)
    plt.close("all")
    _, ax = recipe.plot_recipe(return_fig=True, show=False)
    actual_xlabel = ax.get_xlabel()
    actual_ylabel = ax.get_ylabel()
    assert actual_xlabel == expected_xlabel
    assert actual_ylabel == expected_ylabel


def test_plot_recipe_labels_from_gr_file_overwrite(temp_data_files):
    # C1: A recipe loaded from a .gr file is plotted with xlabel and
    # ylabel given.
    # Expected: The given labels replace the ones inferred from the file.
    gr_file = temp_data_files / "gr_file.gr"
    recipe = build_recipe_from_datafile(gr_file)
    optimize_recipe(recipe)
    plt.close("all")
    fig, ax = recipe.plot_recipe(
        return_fig=True, show=False, xlabel="My X", ylabel="My Y"
    )
    actual_xlabel = ax.get_xlabel()
    actual_ylabel = ax.get_ylabel()
    expected_xlabel = "My X"
    expected_ylabel = "My Y"
    assert actual_xlabel == expected_xlabel
    assert actual_ylabel == expected_ylabel


def test_plot_recipe_reset_all_defaults(build_recipes_one_contribution):
    # C1: Every plot option is given a new default with set_plot_defaults.
    # Expected: plot_recipe called without arguments uses every one of the
    # new defaults.
    expected_defaults = {
        "show_observed": True,
        "show_fit": True,
        "show_diff": True,
        "offset_scale": 0.5,
        "xmin": 1,
        "xmax": 10,
        "figsize": (9, 10),
        "data_style": "-",
        "fit_style": "o",
        "diff_style": "o",
        "data_color": "blue",
        "fit_color": "purple",
        "diff_color": "orange",
        "data_label": "my data",
        "fit_label": "my fit",
        "diff_label": "my diff",
        "xlabel": "my x label",
        "ylabel": "my y label",
        "title": "my title",
        "legend": False,
        "legend_loc": "upper right",
        "markersize": 5,
        "linewidth": 3,
        "alpha": 0.5,
        "show": True,
    }

    recipe, _ = build_recipes_one_contribution
    optimize_recipe(recipe)
    plt.close("all")

    recipe.set_plot_defaults(**expected_defaults)
    fig, ax = recipe.plot_recipe(return_fig=True, show=False)

    actual_title = ax.get_title()
    actual_xlabel = ax.get_xlabel()
    actual_ylabel = ax.get_ylabel()

    expected_title = expected_defaults["title"]
    expected_xlabel = expected_defaults["xlabel"]
    expected_ylabel = expected_defaults["ylabel"]

    assert actual_title == expected_title
    assert actual_xlabel == expected_xlabel
    assert actual_ylabel == expected_ylabel

    actual_labels, actual_line_count = get_labels_and_linecount(ax)

    expected_labels = [
        expected_defaults["data_label"],
        expected_defaults["fit_label"],
        expected_defaults["diff_label"],
    ]
    expected_line_count = 3

    assert actual_line_count == expected_line_count
    assert actual_labels == expected_labels

    lines_by_label = {line.get_label(): line for line in ax.get_lines()}

    data_line = lines_by_label[expected_defaults["data_label"]]
    fit_line = lines_by_label[expected_defaults["fit_label"]]
    diff_line = lines_by_label[expected_defaults["diff_label"]]

    actual_line_styles = {
        "data_color": data_line.get_color(),
        "fit_color": fit_line.get_color(),
        "diff_color": diff_line.get_color(),
        "data_style": data_line.get_linestyle(),
        "fit_style": fit_line.get_marker(),
        "diff_style": diff_line.get_marker(),
        "markersize": data_line.get_markersize(),
        "alpha": data_line.get_alpha(),
    }
    expected_line_styles = {
        key: expected_defaults[key] for key in actual_line_styles
    }
    assert actual_line_styles == expected_line_styles

    actual_xlim = ax.get_xlim()
    expected_xlim = (expected_defaults["xmin"], expected_defaults["xmax"])
    assert actual_xlim == expected_xlim

    # No legend.
    actual_legend = ax.get_legend() is not None
    expected_legend = expected_defaults["legend"]

    assert actual_legend == expected_legend


# Observed uncertainties are optional. A profile loaded without an
# uncertainty column keeps dyobs as None and is refined unweighted, with
# dy falling back to one at every calculation point. The cases below
# refine the same noiseless sine profile so that the known solution
# (A=1, k=1, c=0) is recovered no matter how the fit is weighted.


# The profile of build_recipe_with_uncertainty has 21 points.
UNIFORM_DYOBS = ones(21)
VARYING_DYOBS = linspace(0.1, 1.0, 21)


@pytest.mark.parametrize(
    "input_dyobs, expected_dyobs_is_set",
    [
        # C1: No uncertainties are observed, as for a file with no
        # uncertainty column.
        # Expected: The refinement converges and dyobs stays None.
        (None, False),
        # C2: Uniform uncertainties of one are observed explicitly.
        # Expected: The refinement converges to the same values as C1,
        # since an unweighted fit is weighted by one everywhere.
        (UNIFORM_DYOBS, True),
        # C3: Uncertainties vary across the profile, so the refinement
        # is weighted point by point.
        # Expected: The refinement still converges to the known values.
        (VARYING_DYOBS, True),
    ],
)
def test_refine_with_and_without_uncertainty(
    build_recipe_with_uncertainty, input_dyobs, expected_dyobs_is_set
):
    recipe, profile = build_recipe_with_uncertainty(input_dyobs)
    optimize_recipe(recipe)
    actual_dyobs_is_set = profile.dyobs is not None
    actual_names = recipe.get_names()
    actual_values = recipe.get_values()
    expected_names = ["A", "k", "c"]
    expected_values = [1.0, 1.0, 0.0]
    assert actual_dyobs_is_set == expected_dyobs_is_set
    assert actual_names == expected_names
    assert actual_values == pytest.approx(expected_values, abs=1e-5)


# The residual that is optimized is (ycalc - y)/dy, so dy is what
# actually weights a refinement. These cases check that a profile with
# no observed uncertainties is weighted identically to one whose
# uncertainties are all one, and that observed uncertainties are carried
# through to the residual unchanged.


@pytest.mark.parametrize(
    "input_dyobs, expected_dy",
    [
        # C1: No uncertainties are observed.
        # Expected: dy falls back to one everywhere, so the residual is
        # unweighted.
        (None, UNIFORM_DYOBS),
        # C2: Uniform uncertainties of one are observed explicitly.
        # Expected: dy is one everywhere, matching C1.
        (UNIFORM_DYOBS, UNIFORM_DYOBS),
        # C3: Uncertainties vary across the profile.
        # Expected: dy keeps the observed values, so the residual is
        # weighted point by point.
        (VARYING_DYOBS, VARYING_DYOBS),
    ],
)
def test_residual_is_weighted_by_uncertainty(
    build_recipe_with_uncertainty, input_dyobs, expected_dy
):
    recipe, profile = build_recipe_with_uncertainty(input_dyobs)
    actual_residual = recipe.residual(recipe.values)
    actual_dy = profile.dy
    expected_residual = (profile.ycalc - profile.y) / expected_dy
    assert actual_dy == pytest.approx(expected_dy)
    assert actual_residual == pytest.approx(expected_residual)
