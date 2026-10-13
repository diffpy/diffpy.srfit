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
"""Tests for the constraint module."""

import pytest

from diffpy.srfit.fitbase.constraint import Constraint
from diffpy.srfit.fitbase.recipeorganizer import get_equation_from_string


# Test that add_constraint marks only the constrained Parameter and sets its
# value from the equation.
def test_add_constraint(make_parameters_with_factory):
    # C1: p1 = 1 is constrained to 2*p2, where p2 = 2.
    # Expected: p1 is constrained, p2 is not, and p1 takes the value 4.
    p1, p2, factory = make_parameters_with_factory()
    constraint = Constraint()
    constraint.add_constraint(p1, get_equation_from_string("2*p2", factory))
    actual_constrained = [p1.constrained, p2.constrained]
    expected_constrained = [True, False]
    assert actual_constrained == expected_constrained
    actual_value = p1.get_value()
    expected_value = 4.0
    assert actual_value == expected_value


# Test that add_constraint refuses a Parameter that cannot be constrained.
@pytest.mark.parametrize(
    "input_parameter_state",
    [
        # C1: p1 is already constrained by another Constraint.
        # Expected: Raises a ValueError.
        "constrained",
        # C2: p1 is constant.
        # Expected: Raises a ValueError.
        "constant",
    ],
)
def test_add_constraint_bad(
    make_parameters_with_factory, input_parameter_state
):
    p1, _, factory = make_parameters_with_factory()
    if input_parameter_state == "constrained":
        Constraint().add_constraint(
            p1, get_equation_from_string("2*p2", factory)
        )
    else:
        p1.set_constant()
    equation = get_equation_from_string("2*p2+1", factory)
    with pytest.raises(ValueError):
        Constraint().add_constraint(p1, equation)


# Test that update sets the constrained Parameter from the current values in
# the equation, and that the Parameter is not changed until update is called.
@pytest.mark.parametrize(
    "input_p2_value, input_call_update, expected_p1_value",
    [
        # C1: p2 is set to 2.5 and update is called, for p1 = 2*p2.
        # Expected: p1 becomes 5.0.
        (2.5, True, 5.0),
        # C2: p2 is set to 8.1 and update is called, for p1 = 2*p2.
        # Expected: p1 becomes 16.2.
        (8.1, True, 16.2),
        # C3: p2 is set to 8.1 and update is not called, for p1 = 2*p2.
        # Expected: p1 keeps the value 4.0 it had when it was constrained.
        (8.1, False, 4.0),
    ],
)
def test_update(
    make_parameters_with_factory,
    input_p2_value,
    input_call_update,
    expected_p1_value,
):
    p1, p2, factory = make_parameters_with_factory()
    constraint = Constraint()
    constraint.add_constraint(p1, get_equation_from_string("2*p2", factory))
    p2.set_value(input_p2_value)
    if input_call_update:
        constraint.update()
    actual_p1_value = p1.get_value()
    assert actual_p1_value == expected_p1_value
