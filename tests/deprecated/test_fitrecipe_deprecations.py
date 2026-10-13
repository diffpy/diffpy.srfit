"""Tests for the deprecated camelCase names used with FitRecipe.

Each test calls the deprecated name, checks that it warns, and checks
that it gives the same result as the new name on an identical fresh
object.

Remove this file in version 4.0.0, together with the deprecated names.
"""

import re

import pytest

from diffpy.srfit import FitRecipe, Profile, ProfileParser
from diffpy.srfit.fitbase.fithook import FitHook
from diffpy.srfit.fitbase.parameter import Parameter


def expected_message(base, old_name, new_name):
    return re.escape(
        f"'{base}.{old_name}' is deprecated and will be removed in version "
        f"4.0.0. Please use '{base}.{new_name}' instead."
    )


FIT_RECIPE = "diffpy.srfit.fitbase.FitRecipe"
PROFILE = "diffpy.srfit.fitbase.profile.Profile"
PROFILE_PARSER = "diffpy.srfit.fitbase.profileparser.ProfileParser"
RECIPE_ORGANIZER = "diffpy.srfit.fitbase.recipeorganizer.RecipeOrganizer"


def _add_variables(recipe):
    recipe.add_variable(recipe.cont.A, 2)
    recipe.add_variable(recipe.cont.k, 1)
    recipe.add_variable(recipe.cont.c, 0)


def _describe_variables(recipe):
    return recipe.get_names(), recipe.get_values().tolist()


def _call_add_variable(recipe, method_name):
    method = getattr(recipe, method_name)
    method(recipe.cont.A, 2, tag="tagA")
    recipe.fix("tagA")
    return _describe_variables(recipe), recipe.is_free(recipe.A)


def _call_create_new_variable(recipe, method_name):
    method = getattr(recipe, method_name)
    variable = method("B", 3)
    return _describe_variables(recipe), variable is recipe.B


def _call_delete_variable(recipe, method_name):
    method = getattr(recipe, method_name)
    _add_variables(recipe)
    method(recipe.k)
    return _describe_variables(recipe), hasattr(recipe, "k")


def _call_is_free(recipe, method_name):
    method = getattr(recipe, method_name)
    _add_variables(recipe)
    recipe.fix(recipe.k)
    return method(recipe.A), method(recipe.k)


def _call_get_names(recipe, method_name):
    method = getattr(recipe, method_name)
    _add_variables(recipe)
    return method()


def _call_get_values(recipe, method_name):
    method = getattr(recipe, method_name)
    _add_variables(recipe)
    return method().tolist()


def _call_get_bounds_pairs(recipe, method_name):
    method = getattr(recipe, method_name)
    _add_variables(recipe)
    recipe.A.bounds = (0, 5)
    return method()


def _call_get_bounds_array(recipe, method_name):
    method = getattr(recipe, method_name)
    _add_variables(recipe)
    recipe.A.bounds = (0, 5)
    return [bounds.tolist() for bounds in method()]


def _call_convert_bounds_to_restraints(recipe, method_name):
    method = getattr(recipe, method_name)
    recipe.create_new_variable("var1", 1)
    recipe.var1.bounds = (-1, 1)
    method(sig=2, scaled=True)
    return [
        (
            restraint.lower_bound,
            restraint.upper_bound,
            restraint.sig,
            restraint.scaled,
        )
        for restraint in recipe._restraints
    ]


def _call_add_constraint(recipe, method_name):
    method = getattr(recipe, method_name)
    _add_variables(recipe)
    recipe.create_new_variable("B", 0)
    method(recipe.B, Parameter("Bpar", -1))
    return _describe_variables(recipe), recipe.scalar_residual()


def _call_remove_constraint(recipe, method_name):
    method = getattr(recipe, method_name)
    _add_variables(recipe)
    recipe.add_constraint(recipe.c, "2*A")
    method(recipe.c)
    return (
        _describe_variables(recipe),
        recipe.c.constrained,
        recipe.scalar_residual(),
    )


def _call_add_contribution(recipe, method_name):
    new_recipe = FitRecipe("recipe")
    new_recipe.fithooks[0].verbose = 0
    getattr(new_recipe, method_name)(recipe.cont)
    return list(new_recipe._contributions), new_recipe.scalar_residual()


def _call_set_weight(recipe, method_name):
    method = getattr(recipe, method_name)
    recipe.cont.c.set_value(1)
    method(recipe.cont, 0.5)
    return recipe._weights, recipe.scalar_residual()


def _call_push_fit_hook(recipe, method_name):
    method = getattr(recipe, method_name)
    method(FitHook(), 0)
    return [type(fit_hook).__name__ for fit_hook in recipe.fithooks]


def _call_pop_fit_hook(recipe, method_name):
    method = getattr(recipe, method_name)
    method(recipe.fithooks[0])
    return recipe.fithooks


def _call_get_fit_hooks(recipe, method_name):
    method = getattr(recipe, method_name)
    return [type(fit_hook).__name__ for fit_hook in method()]


def _call_clear_fit_hooks(recipe, method_name):
    method = getattr(recipe, method_name)
    method()
    return recipe.fithooks


def _call_add_parameter_set(recipe, method_name):
    method = getattr(recipe, method_name)
    parameter = Parameter("added_param", 1)
    method(parameter)
    return recipe.added_param is parameter, recipe.added_param.value


def _call_remove_parameter_set(recipe, method_name):
    method = getattr(recipe, method_name)
    parameter = Parameter("added_param", 1)
    recipe.add_parameter_set(parameter)
    method(parameter)
    return hasattr(recipe, "added_param")


def _call_scalar_residual(recipe, method_name):
    method = getattr(recipe, method_name)
    recipe.cont.c.set_value(1)
    return method()


_FIT_RECIPE_CALLS = {
    "add_variable": _call_add_variable,
    "create_new_variable": _call_create_new_variable,
    "delete_variable": _call_delete_variable,
    "is_free": _call_is_free,
    "get_names": _call_get_names,
    "get_values": _call_get_values,
    "get_bounds_pairs": _call_get_bounds_pairs,
    "get_bounds_array": _call_get_bounds_array,
    "convert_bounds_to_restraints": _call_convert_bounds_to_restraints,
    "add_constraint": _call_add_constraint,
    "remove_constraint": _call_remove_constraint,
    "add_contribution": _call_add_contribution,
    "set_weight": _call_set_weight,
    "push_fit_hook": _call_push_fit_hook,
    "pop_fit_hook": _call_pop_fit_hook,
    "get_fit_hooks": _call_get_fit_hooks,
    "clear_fit_hooks": _call_clear_fit_hooks,
    "add_parameter_set": _call_add_parameter_set,
    "remove_parameter_set": _call_remove_parameter_set,
    "scalar_residual": _call_scalar_residual,
}


def _call_fit_recipe_method(build_sine_recipe, method_name, new_name):
    """Call a method of a fresh sine recipe the way its new name is
    tested, and return the result or the resulting state."""
    recipe = build_sine_recipe()
    return _FIT_RECIPE_CALLS[new_name](recipe, method_name)


# FitRecipe: each deprecated method does the same as its new name on a
# fresh sine recipe (see build_sine_recipe in conftest.py). See
# test_fitrecipe.py.
@pytest.mark.parametrize(
    "input_old_name, input_new_name",
    [
        # C1: addVar adds the variable A = 2 with the tag "tagA", which
        # is then fixed.
        # Expected: Warns and gives the same result as add_variable.
        ("addVar", "add_variable"),
        # C2: newVar creates the variable B = 3.
        # Expected: Warns and gives the same result as
        # create_new_variable.
        ("newVar", "create_new_variable"),
        # C3: delVar deletes the variable k.
        # Expected: Warns and gives the same result as delete_variable.
        ("delVar", "delete_variable"),
        # C4: isFree checks the free variable A and the fixed variable k.
        # Expected: Warns and gives the same result as is_free.
        ("isFree", "is_free"),
        # C5: getNames gets the names of the variables A, k and c.
        # Expected: Warns and gives the same result as get_names.
        ("getNames", "get_names"),
        # C6: getValues gets the values of the variables A, k and c.
        # Expected: Warns and gives the same result as get_values.
        ("getValues", "get_values"),
        # C7: getBounds gets the bounds of the variables, where A is
        # bounded to [0, 5].
        # Expected: Warns and gives the same result as get_bounds_pairs.
        ("getBounds", "get_bounds_pairs"),
        # C8: getBounds2 gets the bounds of the variables as two arrays,
        # where A is bounded to [0, 5].
        # Expected: Warns and gives the same result as get_bounds_array.
        ("getBounds2", "get_bounds_array"),
        # C9: boundsToRestraints converts the bounds [-1, 1] of a
        # variable with sig=2 and scaled=True.
        # Expected: Warns and gives the same result as
        # convert_bounds_to_restraints.
        ("boundsToRestraints", "convert_bounds_to_restraints"),
        # C10: constrain constrains the variable B to a parameter.
        # Expected: Warns and gives the same result as add_constraint.
        ("constrain", "add_constraint"),
        # C11: unconstrain removes the constraint c = 2*A.
        # Expected: Warns and gives the same result as remove_constraint.
        ("unconstrain", "remove_constraint"),
        # C12: addContribution adds the sine contribution to an empty
        # recipe.
        # Expected: Warns and gives the same result as add_contribution.
        ("addContribution", "add_contribution"),
        # C13: setWeight sets the weight of the contribution to 0.5.
        # Expected: Warns and gives the same result as set_weight.
        ("setWeight", "set_weight"),
        # C14: pushFitHook inserts a FitHook at the start.
        # Expected: Warns and gives the same result as push_fit_hook.
        ("pushFitHook", "push_fit_hook"),
        # C15: popFitHook removes the default PrintFitHook.
        # Expected: Warns and gives the same result as pop_fit_hook.
        ("popFitHook", "pop_fit_hook"),
        # C16: getFitHooks gets the default PrintFitHook.
        # Expected: Warns and gives the same result as get_fit_hooks.
        ("getFitHooks", "get_fit_hooks"),
        # C17: clearFitHooks removes all fit hooks.
        # Expected: Warns and gives the same result as clear_fit_hooks.
        ("clearFitHooks", "clear_fit_hooks"),
        # C18: addParameterSet adds a parameter.
        # Expected: Warns and gives the same result as add_parameter_set.
        ("addParameterSet", "add_parameter_set"),
        # C19: removeParameterSet removes an added parameter.
        # Expected: Warns and gives the same result as
        # remove_parameter_set.
        ("removeParameterSet", "remove_parameter_set"),
        # C20: scalarResidual calculates the residual with c = 1.
        # Expected: Warns and gives the same result as scalar_residual.
        ("scalarResidual", "scalar_residual"),
    ],
)
def test_fit_recipe_deprecations(
    build_sine_recipe, input_old_name, input_new_name
):
    with pytest.warns(
        DeprecationWarning,
        match=expected_message(FIT_RECIPE, input_old_name, input_new_name),
    ):
        actual_old_result = _call_fit_recipe_method(
            build_sine_recipe, input_old_name, input_new_name
        )
    actual_new_result = _call_fit_recipe_method(
        build_sine_recipe, input_new_name, input_new_name
    )
    assert actual_old_result == actual_new_result


def _parse_gr_file(parse_method_name, file_path, as_list):
    parser = ProfileParser()
    getattr(parser, parse_method_name)(str(file_path))
    return [as_list(data) for data in parser.get_data()], parser.get_metadata()


# ProfileParser: see test_fitrecipe.py::test_plot_recipe_labels_from_gr_file.
def test_parseFile(temp_data_files, as_list):
    # C1: parseFile parses a .gr file.
    # Expected: Warns, and the parser holds the same data and metadata as
    # after parse_file.
    file_path = temp_data_files / "gr_file.gr"
    with pytest.warns(
        DeprecationWarning,
        match=expected_message(PROFILE_PARSER, "parseFile", "parse_file"),
    ):
        actual_old_result = _parse_gr_file("parseFile", file_path, as_list)
    actual_new_result = _parse_gr_file("parse_file", file_path, as_list)
    assert actual_old_result == actual_new_result


def _load_gr_file(load_method_name, file_path, as_list):
    parser = ProfileParser()
    parser.parse_file(str(file_path))
    profile = Profile()
    getattr(profile, load_method_name)(parser)
    return (
        as_list(profile.xobs),
        as_list(profile.yobs),
        as_list(profile.dyobs),
        profile.meta,
    )


# Profile: see test_fitrecipe.py::test_plot_recipe_labels_from_gr_file.
def test_loadParsedData(temp_data_files, as_list):
    # C1: loadParsedData loads a parsed .gr file into a profile.
    # Expected: Warns, and the profile holds the same observed data and
    # metadata as after load_parsed_data.
    file_path = temp_data_files / "gr_file.gr"
    with pytest.warns(
        DeprecationWarning,
        match=expected_message(PROFILE, "loadParsedData", "load_parsed_data"),
    ):
        actual_old_result = _load_gr_file("loadParsedData", file_path, as_list)
    actual_new_result = _load_gr_file("load_parsed_data", file_path, as_list)
    assert actual_old_result == actual_new_result


def _unconstrain_c_with_contribution(build_sine_recipe, method_name):
    recipe = build_sine_recipe()
    recipe.cont.add_constraint(recipe.cont.c, "2*A")
    getattr(recipe.cont, method_name)(recipe.cont.c)
    recipe.cont.c.set_value(0)
    return recipe.cont.c.constrained, recipe.residual().tolist()


# RecipeOrganizer, through FitContribution: see
# test_fitrecipe.py::test_residual.
def test_unconstrain_with_contribution(build_sine_recipe):
    # C1: The contribution of the sine recipe constrains c = 2*A, then
    # unconstrain removes the constraint and c is reset to 0.
    # Expected: Warns, and c and the residual are the same as after
    # remove_constraint.
    with pytest.warns(
        DeprecationWarning,
        match=expected_message(
            RECIPE_ORGANIZER, "unconstrain", "remove_constraint"
        ),
    ):
        actual_old_result = _unconstrain_c_with_contribution(
            build_sine_recipe, "unconstrain"
        )
    actual_new_result = _unconstrain_c_with_contribution(
        build_sine_recipe, "remove_constraint"
    )
    assert actual_old_result == actual_new_result
