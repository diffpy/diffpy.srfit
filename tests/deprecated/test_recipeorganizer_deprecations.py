"""Tests for the deprecated camelCase names in
diffpy.srfit.fitbase.recipeorganizer.

Each test calls a deprecated name, checks that it warns, and checks that
it gives the same result as its new name on an identical fresh object.

Remove this file in version 4.0.0, together with the deprecated names.
"""

import re
import types

import pytest
from numpy import inf

from diffpy.srfit.equation import Equation
from diffpy.srfit.equation.builder import EquationFactory
from diffpy.srfit.fitbase.calculator import Calculator
from diffpy.srfit.fitbase.parameter import Parameter
from diffpy.srfit.fitbase.recipeorganizer import (
    RecipeContainer,
    RecipeOrganizer,
    equationFromString,
    get_equation_from_string,
)
from diffpy.srfit.fitbase.restraint import Restraint


def expected_message(base, old_name, new_name):
    return re.escape(
        f"'{base}.{old_name}' is deprecated and will be removed in version "
        f"4.0.0. Please use '{base}.{new_name}' instead."
    )


RECIPEORGANIZER = "diffpy.srfit.fitbase.recipeorganizer"
RECIPE_CONTAINER = "diffpy.srfit.fitbase.recipeorganizer.RecipeContainer"
RECIPE_ORGANIZER = "diffpy.srfit.fitbase.recipeorganizer.RecipeOrganizer"


def _describe(result):
    """Return a comparable description of what a method returned.

    Parameters are described by their name and value, Equations by their
    value and Restraints by their penalty.
    """
    if isinstance(result, Parameter):
        return result.name, result.value
    if isinstance(result, Equation):
        return result()
    if isinstance(result, Restraint):
        return result.penalty()
    if isinstance(result, (list, types.GeneratorType, type({}.keys()))):
        return [_describe(item) for item in result]
    return result


def _make_container():
    """Return the container root holding x = 1 and y = 2, and the
    container child, which holds x = 3."""
    root = RecipeContainer("root")
    root._containers = {}
    root._manage(root._containers)
    root._add_object(Parameter("x", 1), root._parameters)
    root._add_object(Parameter("y", 2), root._parameters)
    child = RecipeContainer("child")
    child._add_object(Parameter("x", 3), child._parameters)
    root._add_object(child, root._containers)
    return root


# RecipeContainer: each deprecated method does the same as its new name on
# a container holding x = 1, y = 2 and a child container holding x = 3.
# See test_recipeorganizer.py.
@pytest.mark.parametrize(
    "input_old_name, input_new_name, input_args",
    [
        # C1: getNames gets the names of the container's Parameters.
        # Expected: Warns and gives the same names as get_names.
        ("getNames", "get_names", ()),
        # C2: getValues gets the values of the container's Parameters.
        # Expected: Warns and gives the same values as get_values.
        ("getValues", "get_values", ()),
        # C3: iterPars iterates over the Parameters named x in the
        # container and its child.
        # Expected: Warns and gives the same Parameters as
        # iterate_over_parameters.
        ("iterPars", "iterate_over_parameters", ("^x$",)),
    ],
)
def test_recipe_container_deprecations(
    input_old_name, input_new_name, input_args
):
    with pytest.warns(
        DeprecationWarning,
        match=expected_message(
            RECIPE_CONTAINER, input_old_name, input_new_name
        ),
    ):
        actual_old_result = _describe(
            getattr(_make_container(), input_old_name)(*input_args)
        )
    actual_new_result = _describe(
        getattr(_make_container(), input_new_name)(*input_args)
    )
    assert actual_old_result == actual_new_result


class _ScaleCalculator(Calculator):
    """Calculator of scale * x."""

    def __init__(self, name):
        Calculator.__init__(self, name)
        self.new_parameter("scale", 2.0)

    def __call__(self, x):
        return self.scale.get_value() * x


def _double(x):
    return 2 * x


def _make_organizer():
    """Return an organizer holding x = 1, y = 2 and z = 3, with the
    constraint y = 2*x and a soft bound of z <= 2."""
    organizer = RecipeOrganizer("organizer")
    for name, value in (("x", 1), ("y", 2), ("z", 3)):
        organizer._new_parameter(name, value)
    organizer.add_constraint("y", "2*x")
    organizer.add_soft_bounds("z", upper_bound=2)
    return organizer


def _call_organizer_method(organizer, method_name, args):
    """Call a method on an organizer from _make_organizer.

    Return what the method returned, then the Parameters, the
    constrained Parameters and the total restraint penalty of the
    organizer afterward.
    """
    result = getattr(organizer, method_name)(*args)
    constrained_names = sorted(
        parameter.name for parameter in organizer._constraints
    )
    total_penalty = sum(
        restraint.penalty() for restraint in organizer._restraints
    )
    return (
        _describe(result),
        _describe(list(organizer._parameters.values())),
        constrained_names,
        total_penalty,
    )


# RecipeOrganizer: each deprecated method does the same as its new name on
# an organizer holding x = 1, y = 2 and z = 3, with the constraint y = 2*x
# and a soft bound of z <= 2. See test_recipeorganizer.py.
@pytest.mark.parametrize(
    "input_old_name, input_new_name, input_args",
    [
        # C1: registerFunction registers the function 2 * x.
        # Expected: Warns and gives the same result as register_function.
        ("registerFunction", "register_function", (_double,)),
        # C2: registerStringFunction registers f = x + y.
        # Expected: Warns and gives the same result as
        # register_string_function.
        (
            "registerStringFunction",
            "register_string_function",
            ("x + y", "f"),
        ),
        # C3: evaluateEquation evaluates "3 * x".
        # Expected: Warns and gives the same result as evaluate_equation.
        ("evaluateEquation", "evaluate_equation", ("3 * x",)),
        # C4: constrain constrains z to "x + y".
        # Expected: Warns and gives the same result as add_constraint.
        ("constrain", "add_constraint", ("z", "x + y")),
        # C5: isConstrained checks whether y is constrained.
        # Expected: Warns and gives the same result as is_constrained.
        ("isConstrained", "is_constrained", ("y",)),
        # C6: unconstrain removes the constraint on y.
        # Expected: Warns and gives the same result as remove_constraint.
        ("unconstrain", "remove_constraint", ("y",)),
        # C7: getConstrainedPars gets the constrained Parameters.
        # Expected: Warns and gives the same result as
        # get_constrained_parameters.
        ("getConstrainedPars", "get_constrained_parameters", ()),
        # C8: clearConstraints removes all constraints.
        # Expected: Warns and gives the same result as
        # clear_all_constraints.
        ("clearConstraints", "clear_all_constraints", ()),
        # C9: restrain adds the soft bound x <= 0.
        # Expected: Warns and gives the same result as add_soft_bounds.
        ("restrain", "add_soft_bounds", ("x", -inf, 0)),
        # C10: clearRestraints removes all soft bounds.
        # Expected: Warns and gives the same result as
        # clear_all_soft_bounds.
        ("clearRestraints", "clear_all_soft_bounds", ()),
    ],
)
def test_recipe_organizer_deprecations(
    input_old_name, input_new_name, input_args
):
    with pytest.warns(
        DeprecationWarning,
        match=expected_message(
            RECIPE_ORGANIZER, input_old_name, input_new_name
        ),
    ):
        actual_old_result = _call_organizer_method(
            _make_organizer(), input_old_name, input_args
        )
    actual_new_result = _call_organizer_method(
        _make_organizer(), input_new_name, input_args
    )
    assert actual_old_result == actual_new_result


# RecipeOrganizer: the deprecated methods below take an object that is
# made for each organizer. The organizer is the same as in
# test_recipe_organizer_deprecations.


def test_registerCalculator():
    # C1: registerCalculator registers the calculator g = 2 * x.
    # Expected: Warns and gives the same result as register_calculator.
    with pytest.warns(
        DeprecationWarning,
        match=expected_message(
            RECIPE_ORGANIZER, "registerCalculator", "register_calculator"
        ),
    ):
        actual_old_result = _call_organizer_method(
            _make_organizer(), "registerCalculator", (_ScaleCalculator("g"),)
        )
    actual_new_result = _call_organizer_method(
        _make_organizer(), "register_calculator", (_ScaleCalculator("g"),)
    )
    assert actual_old_result == actual_new_result


def test_addRestraint():
    # C1: addRestraint adds a Restraint of x <= 0.
    # Expected: Warns and gives the same result as register_soft_bounds.
    old_organizer = _make_organizer()
    old_restraint = Restraint(Equation(root=old_organizer.x), -inf, 0)
    with pytest.warns(
        DeprecationWarning,
        match=expected_message(
            RECIPE_ORGANIZER, "addRestraint", "register_soft_bounds"
        ),
    ):
        actual_old_result = _call_organizer_method(
            old_organizer, "addRestraint", (old_restraint,)
        )
    new_organizer = _make_organizer()
    new_restraint = Restraint(Equation(root=new_organizer.x), -inf, 0)
    actual_new_result = _call_organizer_method(
        new_organizer, "register_soft_bounds", (new_restraint,)
    )
    assert actual_old_result == actual_new_result


def test_unrestrain():
    # C1: unrestrain removes the organizer's soft bound on z.
    # Expected: Warns and gives the same result as remove_soft_bounds.
    old_organizer = _make_organizer()
    with pytest.warns(
        DeprecationWarning,
        match=expected_message(
            RECIPE_ORGANIZER, "unrestrain", "remove_soft_bounds"
        ),
    ):
        actual_old_result = _call_organizer_method(
            old_organizer, "unrestrain", tuple(old_organizer._restraints)
        )
    new_organizer = _make_organizer()
    actual_new_result = _call_organizer_method(
        new_organizer, "remove_soft_bounds", tuple(new_organizer._restraints)
    )
    assert actual_old_result == actual_new_result


def _build_equation(function):
    """Build "p1+p2+p3" with a factory that has p1 = 1 and p2 = 2
    registered, passing p3 = 3 in the namespace dictionary.

    Return the argument names and value of the equation, and whether p3
    is still registered with the factory afterward.
    """
    factory = EquationFactory()
    factory.register_argument("p1", Parameter("p1", 1))
    factory.register_argument("p2", Parameter("p2", 2))
    equation = function("p1+p2+p3", factory, {"p3": Parameter("p3", 3)})
    return (
        sorted(arg.name for arg in equation.args),
        equation(),
        "p3" in factory.builders,
    )


# recipeorganizer module: see
# test_recipeorganizer.py::test_get_equation_from_string.
def test_equationFromString():
    # C1: "p1+p2+p3" is built with p3 passed in the namespace dictionary.
    # Expected: Warns and gives the same equation as
    # get_equation_from_string.
    with pytest.warns(
        DeprecationWarning,
        match=expected_message(
            RECIPEORGANIZER, "equationFromString", "get_equation_from_string"
        ),
    ):
        actual_old_result = _build_equation(equationFromString)
    actual_new_result = _build_equation(get_equation_from_string)
    assert actual_old_result == actual_new_result
