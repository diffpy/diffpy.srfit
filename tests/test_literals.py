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
"""Tests for the diffpy.srfit.equation.literals module."""

import numpy
import pytest

import diffpy.srfit.equation.literals as literals
import diffpy.srfit.equation.literals.abcs as abcs
from diffpy.srfit.equation.equationmod import Equation


def _make_add_operator():
    return literals.make_operator(
        name="add", symbol="+", operation=numpy.add, nin=2, nout=1
    )


def _make_sum(literal_kind):
    """Return a literal whose value is 1.5 + 2.0 = 3.5."""
    if literal_kind == "argument":
        return literals.Argument(name="a", value=3.5)
    operator = literals.AdditionOperator()
    operator.add_literal(literals.Argument(name="a", value=1.5))
    operator.add_literal(literals.Argument(name="b", value=2.0))
    if literal_kind == "operator":
        return operator
    return Equation(name="eq", root=operator)


def test_argument_init():
    # C1: An Argument is created without any arguments.
    # Expected: It has no value and no name, and it is not constant.
    argument = literals.Argument()
    actual_attributes = (argument._value, argument.const, argument.name)
    expected_attributes = (None, False, None)
    assert actual_attributes == expected_attributes


def test_operator_init():
    # C1: An add operator is made with make_operator.
    # Expected: It stores the given symbol, operation and argument counts,
    # and has no value and no literals yet.
    operator = _make_add_operator()
    actual_attributes = (
        operator.symbol,
        operator.operation,
        operator.nin,
        operator.nout,
        operator._value,
        operator.args,
    )
    expected_attributes = ("+", numpy.add, 2, 1, None, [])
    assert actual_attributes == expected_attributes


# Each literal class is registered with its abstract base class.
@pytest.mark.parametrize(
    "input_class, input_literal, expected_abc",
    [
        # C1: An Argument.
        # Expected: Argument and its instances count as ArgumentABC.
        (literals.Argument, literals.Argument(), abcs.ArgumentABC),
        # C2: An Operator.
        # Expected: Operator and its instances count as OperatorABC.
        (literals.Operator, _make_add_operator(), abcs.OperatorABC),
    ],
)
def test_literal_abc(input_class, input_literal, expected_abc):
    actual_is_abc = (
        issubclass(input_class, expected_abc),
        isinstance(input_literal, expected_abc),
    )
    expected_is_abc = (True, True)
    assert actual_is_abc == expected_is_abc


# get_value and value return the last value given to set_value.
@pytest.mark.parametrize(
    "input_values, expected_value",
    [
        # C1: No value is set.
        # Expected: None.
        ([], None),
        # C2: One value is set.
        # Expected: That value.
        ([3.14], 3.14),
        # C3: Two values are set one after the other.
        # Expected: The second value.
        ([1.0, 3.14], 3.14),
    ],
)
def test_argument_value(input_values, expected_value):
    argument = literals.Argument()
    for value in input_values:
        argument.set_value(value)
    actual_values = (argument.get_value(), argument.value)
    expected_values = (expected_value, expected_value)
    assert actual_values == expected_values


# get_value of every literal kind returns that literal's own value.
@pytest.mark.parametrize(
    "input_literal_kind",
    [
        # C1: An Argument holds the value itself.
        # Expected: 3.5.
        "argument",
        # C2: An Operator computes the value from its literals.
        # Expected: 3.5.
        "operator",
        # C3: An Equation evaluates the operator at its root.
        # Expected: 3.5.
        "equation",
    ],
)
def test_get_value(input_literal_kind):
    literal = _make_sum(input_literal_kind)
    actual_value = literal.get_value()
    expected_value = 3.5
    assert actual_value == expected_value


# An add operator over a=0 and b=0 caches its value. Changing an argument
# clears the cache, and the next value is the new sum.
@pytest.mark.parametrize(
    "input_change, expected_value",
    [
        # C1: a is changed to 4 with set_value.
        # Expected: The cache is cleared and the value is 4.
        ("set_value", 4),
        # C2: b is changed to 2 by assigning to value.
        # Expected: The cache is cleared and the value is 2.
        ("assign_value", 2),
    ],
)
def test_operator_value_after_argument_change(input_change, expected_value):
    operator = _make_add_operator()
    a = literals.Argument(value=0)
    b = literals.Argument(value=0)
    operator.add_literal(a)
    operator.add_literal(b)
    actual_initial_value = operator.value
    expected_initial_value = 0
    assert actual_initial_value == expected_initial_value
    if input_change == "set_value":
        a.set_value(4)
    else:
        b.value = 2
    actual_cached_value = operator._value
    expected_cached_value = None
    assert actual_cached_value == expected_cached_value
    actual_value = operator.value
    assert actual_value == expected_value


# An operator that has fewer literals than it needs has no value.
@pytest.mark.parametrize(
    "input_literal_count",
    [
        # C1: The add operator has no literals.
        # Expected: get_value raises a TypeError.
        0,
        # C2: The add operator has one of its two literals.
        # Expected: get_value raises a TypeError.
        1,
    ],
)
def test_operator_get_value_missing_literals_bad(input_literal_count):
    operator = _make_add_operator()
    for _ in range(input_literal_count):
        operator.add_literal(literals.Argument(value=1))
    with pytest.raises(TypeError):
        operator.get_value()


def test_operator_get_value_returns_cached_value():
    # C1: An operator without literals already has a cached value.
    # Expected: get_value returns the cached value.
    operator = _make_add_operator()
    operator._value = 1
    actual_value = operator.get_value()
    expected_value = 1
    assert actual_value == expected_value


# An operator may not contain itself, directly or through another operator.
@pytest.mark.parametrize(
    "input_cycle",
    [
        # C1: The add operator is added to itself.
        # Expected: A ValueError is raised.
        "self",
        # C2: A subtract operator that contains the add operator is added
        # to the add operator.
        # Expected: A ValueError is raised.
        "containing_operator",
    ],
)
def test_operator_add_literal_cycle_bad(input_cycle):
    add = _make_add_operator()
    add.add_literal(literals.Argument(name="a", value=1))
    subtract = literals.make_operator(
        name="sub", symbol="-", operation=numpy.subtract, nin=2, nout=1
    )
    subtract.add_literal(add)
    literal_to_add = {"self": add, "containing_operator": subtract}
    with pytest.raises(ValueError):
        add.add_literal(literal_to_add[input_cycle])


def test_convolution_operator():
    # C1: Two Gaussians, centred at 4.5 and 2.5 with widths 0.1 and 0.4,
    # are convolved.
    # Expected: A Gaussian centred at the first one's position, with width
    # sqrt(0.1**2 + 0.4**2), scaled to the first one's area.
    x = numpy.linspace(0, 10, 1000)
    g1 = numpy.exp(-0.5 * ((x - 4.5) / 0.1) ** 2)
    g2 = numpy.exp(-0.5 * ((x - 2.5) / 0.4) ** 2)
    operator = literals.ConvolutionOperator()
    operator.add_literal(literals.Argument(name="g1", value=g1))
    operator.add_literal(literals.Argument(name="g2", value=g2))
    actual_convolution = operator.value
    width = (0.1**2 + 0.4**2) ** 0.5
    expected_convolution = numpy.exp(-0.5 * ((x - 4.5) / width) ** 2)
    expected_convolution *= sum(g1) / sum(expected_convolution)
    assert sum(actual_convolution) == pytest.approx(
        sum(expected_convolution), abs=1e-7
    )
    actual_squared_difference = sum(
        (expected_convolution - actual_convolution) ** 2
    )
    assert actual_squared_difference == pytest.approx(0, abs=1e-7)


# ArrayOperator collects the values of its literals x=1, y=2 and z=3 into
# an array.
@pytest.mark.parametrize(
    "input_literal_count, input_new_z_value, expected_array",
    [
        # C1: No literals.
        # Expected: An empty array.
        (0, None, []),
        # C2: Only x.
        # Expected: [1].
        (1, None, [1]),
        # C3: x, y and z.
        # Expected: [1, 2, 3].
        (3, None, [1, 2, 3]),
        # C4: x, y and z, then z is changed to 7.
        # Expected: The array follows the change, [1, 2, 7].
        (3, 7, [1, 2, 7]),
    ],
)
def test_array_operator(
    input_literal_count, input_new_z_value, expected_array
):
    arguments = [
        literals.Argument("x", 1.0),
        literals.Argument("y", 2.0),
        literals.Argument("z", 3.0),
    ]
    operator = literals.ArrayOperator()
    for argument in arguments[:input_literal_count]:
        operator.add_literal(argument)
    if input_new_z_value is not None:
        operator.get_value()
        arguments[2].value = input_new_z_value
    actual_array = operator.value
    assert isinstance(actual_array, numpy.ndarray)
    assert actual_array.tolist() == expected_array
