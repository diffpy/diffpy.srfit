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
"""Tests for the equationmod module."""

import pytest

import diffpy.srfit.equation.literals as literals
from diffpy.srfit.equation import Equation

# Every test builds the equation c*(v1+v3)*(v4-v2) with v1=1, v2=2, v3=3,
# v4=4 and the constant c=2.5, so it evaluates to 2.5*(1+3)*(4-2) = 20.
#
# Every test runs for both kinds of root, given as input_root_kind:
# - "operator": the operator tree is the equation's root.
# - "equation": the operator tree is wrapped in its own Equation, which is
#   then the root (an embedded equation).
ROOT_KINDS = ["operator", "equation"]


def _build_equation(make_args, root_kind):
    v1, v2, v3, v4, c = make_args(5)
    c.name = "c"
    c.const = True
    c.set_value(2.5)
    plus = literals.AdditionOperator()
    plus.add_literal(v1)
    plus.add_literal(v3)
    minus = literals.SubtractionOperator()
    minus.add_literal(v4)
    minus.add_literal(v2)
    mult = literals.MultiplicationOperator()
    mult.add_literal(plus)
    mult.add_literal(minus)
    top = literals.MultiplicationOperator()
    top.add_literal(mult)
    top.add_literal(c)
    root = top if root_kind == "operator" else Equation("root", top)
    equation = Equation("eq", root)
    parts = {
        "v1": v1,
        "v2": v2,
        "v3": v3,
        "v4": v4,
        "c": c,
        "plus": plus,
        "minus": minus,
        "mult": mult,
        "root": root,
    }
    return equation, parts


@pytest.mark.parametrize("input_root_kind", ROOT_KINDS)
def test_equation_args(make_args, input_root_kind):
    # C1: The equation is built.
    # Expected: Its arguments are v1 to v4, but not the constant c.
    equation, _ = _build_equation(make_args, input_root_kind)
    actual_arg_names = sorted(arg.name for arg in equation.args)
    expected_arg_names = ["v1", "v2", "v3", "v4"]
    assert actual_arg_names == expected_arg_names


@pytest.mark.parametrize("input_root_kind", ROOT_KINDS)
def test_equation_root_and_attributes(make_args, input_root_kind):
    # C1: The equation is built.
    # Expected: Its root is the given root, and each argument can be read
    # as an attribute of the same name.
    equation, parts = _build_equation(make_args, input_root_kind)
    actual_root = equation.root
    expected_root = parts["root"]
    assert actual_root is expected_root
    actual_attributes = [equation.v1, equation.v2, equation.v3, equation.v4]
    expected_attributes = [parts["v1"], parts["v2"], parts["v3"], parts["v4"]]
    assert actual_attributes == expected_attributes


# The value can be read in three equivalent ways.
@pytest.mark.parametrize("input_root_kind", ROOT_KINDS)
@pytest.mark.parametrize(
    "input_read",
    [
        # C1: The equation is called.
        # Expected: 20.
        "call",
        # C2: get_value is called.
        # Expected: 20.
        "get_value",
        # C3: The value attribute is read.
        # Expected: 20.
        "value",
    ],
)
def test_equation_value(make_args, input_root_kind, input_read):
    equation, _ = _build_equation(make_args, input_root_kind)
    reads = {
        "call": equation,
        "get_value": equation.get_value,
        "value": lambda: equation.value,
    }
    actual_value = reads[input_read]()
    expected_value = 20
    assert actual_value == expected_value


# Calling the equation with a keyword sets that argument before evaluating.
@pytest.mark.parametrize("input_root_kind", ROOT_KINDS)
@pytest.mark.parametrize(
    "input_name, input_value, expected_value",
    [
        # C1: v1 is set to 2.
        # Expected: 2.5*(2+3)*(4-2) = 25, and v1 keeps the value 2.
        ("v1", 2, 25),
        # C2: v2 is set to 0.
        # Expected: 2.5*(1+3)*(4-0) = 40, and v2 keeps the value 0.
        ("v2", 0, 40),
        # C3: v3 is set to 1.
        # Expected: 2.5*(1+1)*(4-2) = 10, and v3 keeps the value 1.
        ("v3", 1, 10),
        # C4: v4 is set to 0.
        # Expected: 2.5*(1+3)*(0-2) = -20, and v4 keeps the value 0.
        ("v4", 0, -20),
    ],
)
def test_equation_call_with_keyword(
    make_args, input_root_kind, input_name, input_value, expected_value
):
    equation, parts = _build_equation(make_args, input_root_kind)
    actual_value = equation(**{input_name: input_value})
    assert actual_value == expected_value
    actual_argument_value = parts[input_name].value
    expected_argument_value = input_value
    assert actual_argument_value == expected_argument_value


@pytest.mark.parametrize("input_root_kind", ROOT_KINDS)
def test_equation_value_cleared_when_argument_changes(
    make_args, input_root_kind
):
    # C1: The equation is evaluated, then v1 is changed.
    # Expected: The cached values of both the root and the equation are
    # cleared.
    equation, parts = _build_equation(make_args, input_root_kind)
    equation()
    parts["v1"].value = 0
    actual_cached_values = (parts["root"]._value, equation._value)
    expected_cached_values = (None, None)
    assert actual_cached_values == expected_cached_values


@pytest.mark.parametrize("input_root_kind", ROOT_KINDS)
def test_equation_swap_argument(make_args, input_root_kind):
    # C1: The equation is evaluated, then v4 is swapped for v1.
    # Expected: The cached value is cleared, v4 is no longer an argument,
    # and the equation evaluates to 2.5*(1+3)*(1-2) = -10.
    equation, parts = _build_equation(make_args, input_root_kind)
    equation()
    equation.swap(parts["v4"], parts["v1"])
    actual_cached_value = equation._value
    expected_cached_value = None
    assert actual_cached_value == expected_cached_value
    actual_arg_names = sorted(arg.name for arg in equation.args)
    expected_arg_names = ["v1", "v2", "v3"]
    assert actual_arg_names == expected_arg_names
    actual_value = equation()
    expected_value = -10
    assert actual_value == expected_value


# Swapping an argument for a literal that contains that argument would make
# the equation depend on itself.
@pytest.mark.parametrize("input_root_kind", ROOT_KINDS)
@pytest.mark.parametrize(
    "input_old_literal, input_new_literal",
    [
        # C1: v1 is swapped for the root.
        # Expected: A ValueError is raised.
        ("v1", "root"),
        # C2: v1 is swapped for the multiplication (v1+v3)*(v4-v2).
        # Expected: A ValueError is raised.
        ("v1", "mult"),
        # C3: v1 is swapped for the addition v1+v3.
        # Expected: A ValueError is raised.
        ("v1", "plus"),
        # C4: v4 is swapped for the subtraction v4-v2.
        # Expected: A ValueError is raised.
        ("v4", "minus"),
    ],
)
def test_equation_swap_cycle_bad(
    make_args, input_root_kind, input_old_literal, input_new_literal
):
    equation, parts = _build_equation(make_args, input_root_kind)
    with pytest.raises(ValueError):
        equation.swap(parts[input_old_literal], parts[input_new_literal])


@pytest.mark.parametrize("input_root_kind", ROOT_KINDS)
def test_equation_swap_root(make_args, input_root_kind):
    # C1: The equation is evaluated, then its root is swapped for v1.
    # Expected: The cached value is cleared, and the equation evaluates to
    # the value of v1.
    equation, parts = _build_equation(make_args, input_root_kind)
    equation()
    equation.swap(equation.root, parts["v1"])
    actual_cached_value = equation._value
    expected_cached_value = None
    assert actual_cached_value == expected_cached_value
    actual_value = equation()
    expected_value = parts["v1"].value
    assert actual_value == expected_value


@pytest.mark.parametrize("input_root_kind", ROOT_KINDS)
def test_equation_set_root(make_args, input_root_kind):
    # C1: The equation is evaluated, then its root is set to plus = v1+v3.
    # Expected: The cached value is cleared, the arguments are v1 and v3,
    # and the equation evaluates to 1+3 = 4.
    equation, parts = _build_equation(make_args, input_root_kind)
    equation()
    equation.set_root(parts["plus"])
    actual_cached_value = equation._value
    expected_cached_value = None
    assert actual_cached_value == expected_cached_value
    actual_arg_names = [arg.name for arg in equation.args]
    expected_arg_names = ["v1", "v3"]
    assert actual_arg_names == expected_arg_names
    actual_value = equation()
    expected_value = 4
    assert actual_value == expected_value


def test_equation_add_literal_bad(make_args):
    # C1: A literal is added to an equation.
    # Expected: A RuntimeError is raised, because an equation takes its
    # literals from its root.
    equation, parts = _build_equation(make_args, "operator")
    with pytest.raises(RuntimeError):
        equation.add_literal(parts["v1"])
