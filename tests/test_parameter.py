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
"""Tests for the parameter module."""

import numpy as np
import pytest

from diffpy.srfit.fitbase.parameter import (
    Parameter,
    ParameterAdapter,
    ParameterProxy,
)

FINE_ARRAY = np.arange(0, 10, 0.1)
COARSE_ARRAY = np.arange(0, 10, 0.5)


# Test that a value set on a Parameter, with set_value or through the value
# property, is the very object get_value and value return afterward.
@pytest.mark.parametrize(
    "input_initial_value, input_setter, input_value",
    [
        # C1: A scalar is set with set_value on a Parameter with no value.
        # Expected: get_value and value return that scalar.
        (None, "set_value", 3.14),
        # C2: An array is set with set_value on a Parameter holding a
        # scalar.
        # Expected: get_value and value return that same array object.
        (3.14, "set_value", FINE_ARRAY),
        # C3: An array is assigned to the value property of a Parameter
        # holding another array.
        # Expected: get_value and value return the new array object.
        (FINE_ARRAY, "value", COARSE_ARRAY),
        # C4: A scalar is set with set_value on a Parameter holding an
        # array.
        # Expected: get_value and value return that scalar.
        (FINE_ARRAY, "set_value", 1.01),
    ],
)
def test_set_value(input_initial_value, input_setter, input_value):
    parameter = Parameter("l", input_initial_value)
    if input_setter == "set_value":
        parameter.set_value(input_value)
    else:
        parameter.value = input_value
    # The value is stored as given, without being copied.
    actual_is_input_object = [
        parameter.get_value() is input_value,
        parameter.value is input_value,
    ]
    expected_is_input_object = [True, True]
    assert actual_is_input_object == expected_is_input_object


# Test that set_constant sets the constant flag of a Parameter, and sets its
# value only when a value is given.
@pytest.mark.parametrize(
    "input_initial_const, input_args, expected_const, expected_value",
    [
        # C1: A variable Parameter with value 3.14 is set constant with no
        # arguments.
        # Expected: The Parameter is constant and keeps the value 3.14.
        (False, (), True, 3.14),
        # C2: A variable Parameter with value 3.14 is set constant with the
        # value 2.5.
        # Expected: The Parameter is constant and has the value 2.5.
        (False, (True, 2.5), True, 2.5),
        # C3: A constant Parameter with value 3.14 is set not constant.
        # Expected: The Parameter is variable and keeps the value 3.14.
        (True, (False,), False, 3.14),
    ],
)
def test_set_constant(
    input_initial_const, input_args, expected_const, expected_value
):
    parameter = Parameter("l", 3.14, const=input_initial_const)
    parameter.set_constant(*input_args)
    actual_const = parameter.const
    assert actual_const == expected_const
    actual_value = parameter.get_value()
    assert actual_value == expected_value


# Test that a ParameterProxy has its own name but shares its value with the
# Parameter it proxies, whichever of the two is changed.
@pytest.mark.parametrize(
    "input_changed_object, input_new_value, expected_value",
    [
        # C1: Neither the Parameter nor the proxy is changed.
        # Expected: Both have the initial value 3.14.
        (None, None, 3.14),
        # C2: The value of the Parameter is changed to 2.3.
        # Expected: Both have the value 2.3.
        ("parameter", 2.3, 2.3),
        # C3: The value of the proxy is changed to 3.2.
        # Expected: Both have the value 3.2.
        ("proxy", 3.2, 3.2),
    ],
)
def test_parameter_proxy(
    input_changed_object, input_new_value, expected_value
):
    parameter = Parameter("l", 3.14)
    proxy = ParameterProxy("l2", parameter)
    changeable_objects = {"parameter": parameter, "proxy": proxy}
    if input_changed_object is not None:
        changeable_objects[input_changed_object].value = input_new_value
    actual_name = proxy.name
    expected_name = "l2"
    assert actual_name == expected_name
    actual_values = [
        parameter.get_value(),
        parameter.value,
        proxy.get_value(),
        proxy.value,
    ]
    expected_values = [expected_value] * 4
    assert actual_values == expected_values


ACCESSOR_KWARGS = {
    "getter": Parameter.get_value,
    "setter": Parameter.set_value,
}
ATTRIBUTE_KWARGS = {"attr": "value"}


# Test that a ParameterAdapter wrapping a Parameter, through accessor
# functions or through an attribute name, reads and writes the value of the
# wrapped Parameter.
@pytest.mark.parametrize(
    "input_adapter_kwargs, input_changed_object, input_new_value, "
    "expected_attr, expected_value",
    [
        # C1: The adapter uses get_value and set_value, and nothing is
        # changed.
        # Expected: The adapter has no attr and the initial value 3.14.
        (ACCESSOR_KWARGS, None, None, None, 3.14),
        # C2: The adapter uses get_value and set_value, and the Parameter is
        # set to 2.3.
        # Expected: The adapter reads the new value 2.3.
        (ACCESSOR_KWARGS, "parameter", 2.3, None, 2.3),
        # C3: The adapter uses get_value and set_value, and the adapter is
        # set to 3.2.
        # Expected: The Parameter is set to 3.2.
        (ACCESSOR_KWARGS, "adapter", 3.2, None, 3.2),
        # C4: The adapter uses the attribute "value", and nothing is
        # changed.
        # Expected: The adapter has attr "value" and the initial value 3.14.
        (ATTRIBUTE_KWARGS, None, None, "value", 3.14),
        # C5: The adapter uses the attribute "value", and the Parameter is
        # set to 2.3.
        # Expected: The adapter reads the new value 2.3.
        (ATTRIBUTE_KWARGS, "parameter", 2.3, "value", 2.3),
        # C6: The adapter uses the attribute "value", and the adapter is set
        # to 3.2.
        # Expected: The Parameter is set to 3.2.
        (ATTRIBUTE_KWARGS, "adapter", 3.2, "value", 3.2),
    ],
)
def test_parameter_adapter(
    input_adapter_kwargs,
    input_changed_object,
    input_new_value,
    expected_attr,
    expected_value,
):
    parameter = Parameter("l", 3.14)
    adapter = ParameterAdapter("l", parameter, **input_adapter_kwargs)
    changeable_objects = {"parameter": parameter, "adapter": adapter}
    if input_changed_object is not None:
        changeable_objects[input_changed_object].set_value(input_new_value)
    actual_name = adapter.name
    expected_name = parameter.name
    assert actual_name == expected_name
    actual_attr = adapter.attr
    assert actual_attr == expected_attr
    actual_values = [parameter.get_value(), adapter.get_value()]
    expected_values = [expected_value, expected_value]
    assert actual_values == expected_values


# Test that bound_range sets the bounds it is given and leaves a bound
# unchanged when it is given None for it.
@pytest.mark.parametrize(
    "input_initial_bounds, input_lower_bound, input_upper_bound, "
    "expected_bounds",
    [
        # C1: Both bounds are given for a Parameter with no bounds.
        # Expected: Both bounds are set.
        ([-np.inf, np.inf], 1, 10, [1, 10]),
        # C2: Only a lower bound is given for a Parameter with no bounds.
        # Expected: The lower bound is set and the upper bound stays inf.
        ([-np.inf, np.inf], 2, None, [2, np.inf]),
        # C3: Only an upper bound is given for a Parameter with no bounds.
        # Expected: The upper bound is set and the lower bound stays -inf.
        ([-np.inf, np.inf], None, 8, [-np.inf, 8]),
        # C4: Both bounds are given for a Parameter with bounds [0, 10].
        # Expected: Both existing bounds are overwritten.
        ([0, 10], 2, 6, [2, 6]),
    ],
)
def test_bound_range(
    input_initial_bounds, input_lower_bound, input_upper_bound, expected_bounds
):
    parameter = Parameter("a", value=5)
    parameter.bounds = list(input_initial_bounds)
    parameter.bound_range(
        lower_bound=input_lower_bound, upper_bound=input_upper_bound
    )
    actual_bounds = parameter.bounds
    assert actual_bounds == expected_bounds


# Test that bound_window sets bounds around the current value of the
# Parameter, using the lower radius for the upper one when no upper radius
# is given.
@pytest.mark.parametrize(
    "input_initial_value, input_current_value, input_lower_radius, "
    "input_upper_radius, expected_bounds",
    [
        # C1: Only a lower radius of 2 is given, around the value 10.
        # Expected: The window is symmetric, [8, 12].
        (10, 10, 2, None, [8, 12]),
        # C2: A lower radius of 3 and an upper radius of 5 are given, around
        # the value 10.
        # Expected: The window is asymmetric, [7, 15].
        (10, 10, 3, 5, [7, 15]),
        # C3: A radius of 0 is given, around the value 4.
        # Expected: Both bounds equal the value, [4, 4].
        (4, 4, 0, None, [4, 4]),
        # C4: The value is changed from 5 to 20 before a radius of 2 is
        # given.
        # Expected: The window is around the current value, [18, 22].
        (5, 20, 2, None, [18, 22]),
    ],
)
def test_bound_window(
    input_initial_value,
    input_current_value,
    input_lower_radius,
    input_upper_radius,
    expected_bounds,
):
    parameter = Parameter("a", value=input_initial_value)
    parameter.set_value(input_current_value)
    parameter.bound_window(
        lower_radius=input_lower_radius, upper_radius=input_upper_radius
    )
    actual_bounds = parameter.bounds
    assert actual_bounds == expected_bounds
