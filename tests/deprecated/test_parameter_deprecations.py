"""Tests for the deprecated camelCase names in
diffpy.srfit.fitbase.parameter.

Each test calls the deprecated name, checks that it warns, and checks
that it gives the same result as the new name on an identical object.

Remove this file in version 4.0.0, together with the deprecated names.
"""

import re

import numpy as np
import pytest

from diffpy.srfit.fitbase.parameter import (
    Parameter,
    ParameterAdapter,
    ParameterProxy,
)


def expected_message(base, old_name, new_name):
    return re.escape(
        f"'{base}.{old_name}' is deprecated and will be removed in version "
        f"4.0.0. Please use '{base}.{new_name}' instead."
    )


PARAMETER = "diffpy.srfit.fitbase.Parameter"


class _ValueHolder:
    """A plain object for ParameterAdapter to wrap."""

    def __init__(self, value):
        self.value = value


def _make_parameter(parameter_kind):
    """Return a Parameter of the given kind whose value is 3.14."""
    if parameter_kind == "parameter":
        return Parameter("l", 3.14)
    if parameter_kind == "proxy":
        return ParameterProxy("l_proxy", Parameter("l", 3.14))
    return ParameterAdapter("l_adapted", _ValueHolder(3.14), attr="value")


# Parameter: getValue forwards to get_value for Parameter and both of its
# wrapping subclasses.
@pytest.mark.parametrize(
    "input_parameter_kind",
    [
        # C1: A Parameter stores the value itself.
        # Expected: Warns and returns the same value as get_value.
        "parameter",
        # C2: A ParameterProxy defers to the Parameter it proxies.
        # Expected: Warns and returns the same value as get_value.
        "proxy",
        # C3: A ParameterAdapter defers to the attribute it wraps.
        # Expected: Warns and returns the same value as get_value.
        "adapter",
    ],
)
def test_getValue(input_parameter_kind):
    with pytest.warns(
        DeprecationWarning,
        match=expected_message(PARAMETER, "getValue", "get_value"),
    ):
        actual_old_value = _make_parameter(input_parameter_kind).getValue()
    actual_new_value = _make_parameter(input_parameter_kind).get_value()
    assert actual_old_value == actual_new_value


def _call_parameter_method(method_name, args):
    """Call a method of a Parameter whose value is 3.14.

    Return whether the method returned the Parameter, and the value,
    constant flag and bounds the Parameter has afterward.
    """
    parameter = Parameter("l", 3.14)
    result = getattr(parameter, method_name)(*args)
    return (
        result is parameter,
        parameter.get_value(),
        parameter.const,
        list(parameter.bounds),
    )


# Parameter: setValue and setConst do the same as their new names. See
# test_parameter.py::test_set_value and test_parameter.py::test_set_constant.
@pytest.mark.parametrize(
    "input_old_name, input_new_name, input_args",
    [
        # C1: setValue sets the value to 2.5.
        # Expected: Warns and gives the same result as set_value.
        ("setValue", "set_value", (2.5,)),
        # C2: setConst marks the Parameter as constant.
        # Expected: Warns and gives the same result as set_constant.
        ("setConst", "set_constant", ()),
    ],
)
def test_parameter_deprecations(input_old_name, input_new_name, input_args):
    with pytest.warns(
        DeprecationWarning,
        match=expected_message(PARAMETER, input_old_name, input_new_name),
    ):
        actual_old_result = _call_parameter_method(input_old_name, input_args)
    actual_new_result = _call_parameter_method(input_new_name, input_args)
    assert actual_old_result == actual_new_result


def _make_bounded_parameter(initial_bounds):
    parameter = Parameter("a", value=5)
    parameter.bounds = list(initial_bounds)
    return parameter


# Parameter: see test_parameter.py::test_bound_range.
@pytest.mark.parametrize(
    "input_initial_bounds, input_lower_bound, input_upper_bound",
    [
        # C1: Both bounds are given for a Parameter with no bounds.
        # Expected: Warns and sets the same bounds as bound_range.
        ([-np.inf, np.inf], 1, 10),
        # C2: Only a lower bound is given for a Parameter with no bounds.
        # Expected: Warns and sets the same bounds as bound_range.
        ([-np.inf, np.inf], 2, None),
        # C3: Only an upper bound is given for a Parameter with no bounds.
        # Expected: Warns and sets the same bounds as bound_range.
        ([-np.inf, np.inf], None, 8),
        # C4: Both bounds are given for a Parameter with bounds [0, 10].
        # Expected: Warns and sets the same bounds as bound_range.
        ([0, 10], 2, 6),
    ],
)
def test_boundRange(
    input_initial_bounds, input_lower_bound, input_upper_bound
):
    old_parameter = _make_bounded_parameter(input_initial_bounds)
    with pytest.warns(
        DeprecationWarning,
        match=expected_message(PARAMETER, "boundRange", "bound_range"),
    ):
        old_result = old_parameter.boundRange(
            lb=input_lower_bound, ub=input_upper_bound
        )
    new_parameter = _make_bounded_parameter(input_initial_bounds)
    new_result = new_parameter.bound_range(
        lower_bound=input_lower_bound, upper_bound=input_upper_bound
    )
    actual_old_result = (old_result is old_parameter, old_parameter.bounds)
    actual_new_result = (new_result is new_parameter, new_parameter.bounds)
    assert actual_old_result == actual_new_result


# Parameter: see test_parameter.py::test_bound_window.
@pytest.mark.parametrize(
    "input_value, input_lower_radius, input_upper_radius",
    [
        # C1: Only a lower radius of 2 is given, around the value 10.
        # Expected: Warns and sets the same bounds as bound_window.
        (10, 2, None),
        # C2: A lower radius of 3 and an upper radius of 5 are given,
        # around the value 10.
        # Expected: Warns and sets the same bounds as bound_window.
        (10, 3, 5),
        # C3: A radius of 0 is given, around the value 4.
        # Expected: Warns and sets the same bounds as bound_window.
        (4, 0, None),
        # C4: The value is changed from 5 to 20 before a radius of 2 is
        # given.
        # Expected: Warns and sets the same bounds as bound_window.
        (20, 2, None),
    ],
)
def test_boundWindow(input_value, input_lower_radius, input_upper_radius):
    old_parameter = Parameter("a", value=5).set_value(input_value)
    with pytest.warns(
        DeprecationWarning,
        match=expected_message(PARAMETER, "boundWindow", "bound_window"),
    ):
        old_result = old_parameter.boundWindow(
            lr=input_lower_radius, ur=input_upper_radius
        )
    new_parameter = Parameter("a", value=5).set_value(input_value)
    new_result = new_parameter.bound_window(
        lower_radius=input_lower_radius, upper_radius=input_upper_radius
    )
    actual_old_result = (old_result is old_parameter, old_parameter.bounds)
    actual_new_result = (new_result is new_parameter, new_parameter.bounds)
    assert actual_old_result == actual_new_result
