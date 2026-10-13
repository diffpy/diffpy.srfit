"""Tests for the deprecated camelCase names in
diffpy.srfit.fitbase.constraint.

Each test calls the deprecated name, checks that it warns, and checks
that it gives the same result as the new name on an identical object.

Remove this file in version 4.0.0, together with the deprecated names.
"""

import re

import pytest

from diffpy.srfit.fitbase.constraint import Constraint
from diffpy.srfit.fitbase.recipeorganizer import get_equation_from_string


def expected_message(base, old_name, new_name):
    return re.escape(
        f"'{base}.{old_name}' is deprecated and will be removed in version "
        f"4.0.0. Please use '{base}.{new_name}' instead."
    )


CONSTRAINT = "diffpy.srfit.fitbase.constraint.Constraint"


def _call_constraint_method(make_parameters_with_factory, method_name):
    """Constrain p1 = 2*p2 with the given method, or remove that
    constraint with it, where p1 = 1 and p2 = 2.

    Return whether the Constraint holds p1 and the equation, the
    constrained flags of p1 and p2, and the value of p1 afterward.
    """
    p1, p2, factory = make_parameters_with_factory()
    equation = get_equation_from_string("2*p2", factory)
    constraint = Constraint()
    if method_name in ("constrain", "add_constraint"):
        getattr(constraint, method_name)(p1, equation)
    else:
        constraint.add_constraint(p1, equation)
        getattr(constraint, method_name)()
    return (
        constraint.par is p1,
        constraint.eq is equation,
        p1.constrained,
        p2.constrained,
        p1.get_value(),
    )


# Constraint: each deprecated method does the same as its new name. See
# test_constraint.py::test_add_constraint.
@pytest.mark.parametrize(
    "input_old_name, input_new_name",
    [
        # C1: constrain constrains p1 to 2*p2.
        # Expected: Warns and gives the same result as add_constraint.
        ("constrain", "add_constraint"),
        # C2: unconstrain removes the constraint p1 = 2*p2.
        # Expected: Warns and gives the same result as remove_constraint.
        ("unconstrain", "remove_constraint"),
    ],
)
def test_constraint_deprecations(
    make_parameters_with_factory, input_old_name, input_new_name
):
    with pytest.warns(
        DeprecationWarning,
        match=expected_message(CONSTRAINT, input_old_name, input_new_name),
    ):
        actual_old_result = _call_constraint_method(
            make_parameters_with_factory, input_old_name
        )
    actual_new_result = _call_constraint_method(
        make_parameters_with_factory, input_new_name
    )
    assert actual_old_result == actual_new_result


# Constraint: see test_constraint.py::test_add_constraint_bad.
@pytest.mark.parametrize(
    "input_parameter_state",
    [
        # C1: p1 is already constrained by another Constraint.
        # Expected: Warns and raises a ValueError.
        "constrained",
        # C2: p1 is constant.
        # Expected: Warns and raises a ValueError.
        "constant",
    ],
)
def test_constrain_bad(make_parameters_with_factory, input_parameter_state):
    p1, _, factory = make_parameters_with_factory()
    if input_parameter_state == "constrained":
        Constraint().add_constraint(
            p1, get_equation_from_string("2*p2", factory)
        )
    else:
        p1.set_constant()
    equation = get_equation_from_string("2*p2+1", factory)
    with (
        pytest.warns(
            DeprecationWarning,
            match=expected_message(CONSTRAINT, "constrain", "add_constraint"),
        ),
        pytest.raises(ValueError),
    ):
        Constraint().constrain(p1, equation)
