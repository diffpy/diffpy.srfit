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
"""Tests for the profilegenerator module."""

import pickle

import pytest
from numpy import arange, array_equal

from diffpy.srfit import Profile


def test_operation(generator_and_profile):
    # C1: operation is called, and the generator is called with 2*x.
    # Expected: The base ProfileGenerator returns its input x unchanged,
    # so operation returns the profile's calculation points x and the
    # call returns 2*x.
    generator, profile = generator_and_profile
    actual_operation_value = generator.operation()
    expected_operation_value = profile.x
    assert array_equal(actual_operation_value, expected_operation_value)
    actual_call_value = generator(2 * profile.x)
    expected_call_value = 2 * profile.x
    assert array_equal(actual_call_value, expected_call_value)


# A change in the profile flushes the cached value of the generator, so its
# next value is computed on the new calculation points.
@pytest.mark.parametrize(
    "input_update, input_new_x",
    [
        # C1: The profile's calculation points are set with
        # set_calculation_points.
        # Expected: The cached value is cleared, and the value is the new x.
        ("set_calculation_points", arange(0, 9, 0.1)),
        # C2: The profile's x is assigned directly.
        # Expected: The cached value is cleared, and the value is the new x.
        ("assign_x", arange(0, 9, 0.1) + 3),
        # C3: A new profile is assigned with set_profile.
        # Expected: The cached value is cleared, and the value is the new
        # profile's x.
        ("set_profile", arange(0, 8, 0.1)),
    ],
)
def test_update(generator_and_profile, input_update, input_new_x):
    generator, profile = generator_and_profile
    # Evaluate once so there is a cached value to flush.
    actual_initial_value = generator.value
    expected_initial_value = profile.x
    assert array_equal(actual_initial_value, expected_initial_value)
    if input_update == "set_calculation_points":
        profile.set_calculation_points(input_new_x)
    elif input_update == "assign_x":
        profile.x = input_new_x
    else:
        new_profile = Profile()
        new_profile.set_calculation_points(input_new_x)
        generator.set_profile(new_profile)
    actual_cached_value = generator._value
    expected_cached_value = None
    assert actual_cached_value == expected_cached_value
    actual_value = generator.value
    expected_value = input_new_x
    assert array_equal(actual_value, expected_value)


def test_pickling(generator_and_profile):
    # C1: A generator with a profile is pickled and unpickled.
    # Expected: The copy keeps the name and evaluates the same profile.
    generator, profile = generator_and_profile
    unpickled_generator = pickle.loads(pickle.dumps(generator))
    actual_name = unpickled_generator.name
    expected_name = "test"
    assert actual_name == expected_name
    actual_operation_value = unpickled_generator.operation()
    expected_operation_value = profile.x
    assert array_equal(actual_operation_value, expected_operation_value)
    actual_call_value = unpickled_generator(3 * profile.x)
    expected_call_value = 3 * profile.x
    assert array_equal(actual_call_value, expected_call_value)
