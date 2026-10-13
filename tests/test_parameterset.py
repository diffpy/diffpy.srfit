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
"""Tests for the parameterset module."""

import pytest

from diffpy.srfit.fitbase.parameter import Parameter
from diffpy.srfit.fitbase.parameterset import ParameterSet


# Test that a stored ParameterSet, and a Parameter added to it afterward,
# can be reached as attributes of the parent set.
def test_add_parameter_set():
    # C1: The ParameterSet "child" is stored in "test", then the Parameter
    # "p1" is added to "child".
    # Expected: test.child is the stored set and test.child.p1 is p1.
    parset = ParameterSet("test")
    child = ParameterSet("child")
    parset.add_parameter_set(child)
    p1 = Parameter("p1", 1)
    child.add_parameter(p1)
    actual_child = parset.child
    expected_child = child
    assert actual_child == expected_child
    actual_p1 = parset.child.p1
    expected_p1 = p1
    assert actual_p1 == expected_p1


# Test that a ParameterSet refuses to store an object under a name that is
# already taken.
@pytest.mark.parametrize(
    "input_object",
    [
        # C1: A Parameter is named like the stored ParameterSet "child".
        # Expected: Raises a ValueError.
        Parameter("child", 1),
        # C2: Another ParameterSet is named like the stored ParameterSet
        # "child".
        # Expected: Raises a ValueError.
        ParameterSet("child"),
    ],
)
def test_add_parameter_set_name_taken(input_object):
    parset = ParameterSet("test")
    parset.add_parameter_set(ParameterSet("child"))
    with pytest.raises(ValueError):
        parset.add_parameter_set(input_object)


# Test that remove_parameter_set removes a stored ParameterSet.
def test_remove_parameter_set():
    # C1: The stored ParameterSet "child" is removed from "test".
    # Expected: "test" no longer has a "child" attribute.
    parset = ParameterSet("test")
    child = ParameterSet("child")
    parset.add_parameter_set(child)
    parset.remove_parameter_set(child)
    actual_has_child = hasattr(parset, "child")
    expected_has_child = False
    assert actual_has_child == expected_has_child


# Test that set_constant sets the constant flag of every Parameter in the
# set, including those in nested sets.
@pytest.mark.parametrize(
    "input_initial_const, input_args, expected_consts",
    [
        # C1: A set of variable Parameters is set constant with no
        # arguments.
        # Expected: Both p1 and the nested p2 are constant.
        (False, (), [True, True]),
        # C2: A set of constant Parameters is set not constant.
        # Expected: Both p1 and the nested p2 are variable.
        (True, (False,), [False, False]),
    ],
)
def test_set_constant(input_initial_const, input_args, expected_consts):
    parset = ParameterSet("test")
    parset.add_parameter(Parameter("p1", 1, const=input_initial_const))
    child = ParameterSet("child")
    child.add_parameter(Parameter("p2", 2, const=input_initial_const))
    parset.add_parameter_set(child)
    parset.set_constant(*input_args)
    actual_consts = [parset.p1.const, parset.child.p2.const]
    assert actual_consts == expected_consts


# Test that add_parameter stores a Parameter under its name.
def test_add_parameter():
    # C1: The Parameter "p1" is added to an empty set.
    # Expected: The set holds p1 under the name "p1".
    parset = ParameterSet("test")
    p1 = Parameter("p1", 1)
    parset.add_parameter(p1)
    actual_names = parset.get_names()
    expected_names = ["p1"]
    assert actual_names == expected_names
    actual_p1 = parset.p1
    expected_p1 = p1
    assert actual_p1 == expected_p1


# Test that new_parameter creates, stores and returns a Parameter.
def test_new_parameter():
    # C1: The Parameter "p1" with value 1 is created in an empty set.
    # Expected: The set holds a Parameter "p1" with value 1, and it is the
    # Parameter that was returned.
    parset = ParameterSet("test")
    p1 = parset.new_parameter("p1", 1)
    actual_names = parset.get_names()
    expected_names = ["p1"]
    assert actual_names == expected_names
    actual_values = parset.get_values()
    expected_values = [1]
    assert actual_values == expected_values
    actual_p1 = parset.p1
    expected_p1 = p1
    assert actual_p1 == expected_p1


# Test that remove_parameter removes a stored Parameter from the set.
def test_remove_parameter():
    # C1: The stored Parameter "p1" is removed.
    # Expected: The set holds no Parameters.
    parset = ParameterSet("test")
    p1 = parset.new_parameter("p1", 1)
    parset.remove_parameter(p1)
    actual_names = parset.get_names()
    expected_names = []
    assert actual_names == expected_names
