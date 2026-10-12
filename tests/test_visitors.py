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
"""Tests for the visitors module."""

import numpy
import pytest

import diffpy.srfit.equation.literals as literals
import diffpy.srfit.equation.visitors as visitors
from diffpy.srfit.equation import Equation
from diffpy.srfit.equation.visitors.visitor import Visitor


def _build_tree(make_args):
    """Build (v1+v3)*(v4-v2) with v1=1, v2=2, v3=3 and v4=4, which is 8.

    Also return v5=5, which is not part of the tree.
    """
    v1, v2, v3, v4, v5 = make_args(5)
    plus = literals.AdditionOperator()
    plus.add_literal(v1)
    plus.add_literal(v3)
    minus = literals.SubtractionOperator()
    minus.add_literal(v4)
    minus.add_literal(v2)
    mult = literals.MultiplicationOperator()
    mult.add_literal(plus)
    mult.add_literal(minus)
    parts = {
        "v1": v1,
        "v2": v2,
        "v3": v3,
        "v4": v4,
        "v5": v5,
        "plus": plus,
        "minus": minus,
    }
    return mult, parts


# The Validator counts the problems in (v1+v3)*(v4-...). Each problem is an
# operator with the wrong number of literals, or an operator missing its
# name, symbol or operation.
@pytest.mark.parametrize(
    "input_minus_literals, input_plus_missing, expected_error_count",
    [
        # C1: minus has one of its two literals, and plus has no name,
        # symbol or operation.
        # Expected: 4 errors.
        (["v4"], ["name", "symbol", "operation"], 4),
        # C2: minus is complete, and plus has no name, symbol or operation.
        # Expected: 3 errors.
        (["v4", "v3"], ["name", "symbol", "operation"], 3),
        # C3: minus is complete, and plus has no symbol or operation.
        # Expected: 2 errors.
        (["v4", "v3"], ["symbol", "operation"], 2),
        # C4: minus is complete, and plus has no operation.
        # Expected: 1 error.
        (["v4", "v3"], ["operation"], 1),
        # C5: Every operator is complete.
        # Expected: No errors.
        (["v4", "v3"], [], 0),
        # C6: minus has three literals but takes two.
        # Expected: 1 error.
        (["v4", "v3", "v1"], [], 1),
    ],
)
def test_validator(
    make_args, input_minus_literals, input_plus_missing, expected_error_count
):
    v1, v2, v3, v4 = make_args(4)
    arguments = {"v1": v1, "v2": v2, "v3": v3, "v4": v4}
    plus = literals.AdditionOperator()
    plus.add_literal(v1)
    plus.add_literal(v3)
    for attribute in input_plus_missing:
        setattr(plus, attribute, None)
    minus = literals.SubtractionOperator()
    for name in input_minus_literals:
        minus.add_literal(arguments[name])
    mult = literals.MultiplicationOperator()
    mult.add_literal(plus)
    mult.add_literal(minus)
    validator = visitors.Validator()
    mult.identify(validator)
    actual_error_count = len(validator.errors)
    assert actual_error_count == expected_error_count


def test_validator_reset(make_args):
    # C1: A Validator checks a broken tree, is reset, and checks a valid
    # tree.
    # Expected: The errors from the first check are cleared.
    mult, parts = _build_tree(make_args)
    parts["plus"].operation = None
    validator = visitors.Validator()
    mult.identify(validator)
    parts["plus"].operation = numpy.add
    validator.reset()
    mult.identify(validator)
    actual_errors = validator.errors
    expected_errors = []
    assert actual_errors == expected_errors


# get_args returns every Argument in a literal.
@pytest.mark.parametrize(
    "input_literal, expected_arg_names",
    [
        # C1: The tree (v1+v3)*(v4-v2).
        # Expected: v1, v2, v3 and v4.
        ("tree", ["v1", "v2", "v3", "v4"]),
        # C2: The single argument v1.
        # Expected: v1 only.
        ("v1", ["v1"]),
    ],
)
def test_get_args(make_args, input_literal, expected_arg_names):
    mult, parts = _build_tree(make_args)
    literal = mult if input_literal == "tree" else parts[input_literal]
    actual_arg_names = sorted(arg.name for arg in visitors.get_args(literal))
    assert actual_arg_names == expected_arg_names


def test_swap_argument(make_args):
    # C1: The tree is evaluated, then v2 is swapped for v5.
    # Expected: The cached value is cleared, minus now observes v5 instead
    # of v2, the arguments are v1, v3, v4 and v5, and the value is
    # (1+3)*(4-5) = -4.
    mult, parts = _build_tree(make_args)
    v2, v5, minus = parts["v2"], parts["v5"], parts["minus"]
    mult.get_value()
    visitors.swap(mult, v2, v5)
    actual_cached_value = mult._value
    expected_cached_value = None
    assert actual_cached_value == expected_cached_value
    actual_observed = (
        v2.has_observer(minus._flush),
        v5.has_observer(minus._flush),
    )
    expected_observed = (False, True)
    assert actual_observed == expected_observed
    actual_arg_names = sorted(arg.name for arg in visitors.get_args(mult))
    expected_arg_names = ["v1", "v3", "v4", "v5"]
    assert actual_arg_names == expected_arg_names
    actual_value = mult.value
    expected_value = -4
    assert actual_value == expected_value


def test_swap_operator(make_args):
    # C1: The tree is evaluated, then minus is swapped for an empty add
    # operator, which is then given v4 and v5.
    # Expected: The cached value is cleared, mult now observes the new
    # operator instead of minus, there is no value until the new operator
    # has its literals, and then the value is (1+3)*(4+5) = 36.
    mult, parts = _build_tree(make_args)
    minus = parts["minus"]
    mult.get_value()
    plus2 = literals.AdditionOperator()
    visitors.swap(mult, minus, plus2)
    actual_cached_value = mult._value
    expected_cached_value = None
    assert actual_cached_value == expected_cached_value
    actual_observed = (
        minus.has_observer(mult._flush),
        plus2.has_observer(mult._flush),
    )
    expected_observed = (False, True)
    assert actual_observed == expected_observed
    with pytest.raises(TypeError):
        mult.get_value()
    plus2.add_literal(parts["v4"])
    plus2.add_literal(parts["v5"])
    actual_value = mult.value
    expected_value = 36
    assert actual_value == expected_value


# get_expression writes a literal as an equation string. Embedded
# equations whose names match eqskip are written by name.
@pytest.mark.parametrize(
    "input_eqskip, expected_expression",
    [
        # C1: No eqskip is given.
        # Expected: The embedded equation "inner" is written out in full.
        (None, "(((v1 + v3) * (v4 - v2)) + v5)"),
        # C2: eqskip matches the name of the embedded equation.
        # Expected: The embedded equation is written as "inner".
        ("inner", "(inner + v5)"),
        # C3: eqskip does not match the name of the embedded equation.
        # Expected: The embedded equation is written out in full.
        ("other", "(((v1 + v3) * (v4 - v2)) + v5)"),
    ],
)
def test_get_expression(make_args, input_eqskip, expected_expression):
    mult, parts = _build_tree(make_args)
    top = literals.AdditionOperator()
    top.add_literal(Equation("inner", mult))
    top.add_literal(parts["v5"])
    actual_expression = visitors.get_expression(top, eqskip=input_eqskip)
    assert actual_expression == expected_expression


def test_pretty_print(make_args, capsys):
    # C1: The tree (v1+v3)*(v4-v2) is printed.
    # Expected: Its expression is printed on one line.
    mult, _ = _build_tree(make_args)
    visitors.pretty_print(mult)
    actual_output = capsys.readouterr().out
    expected_output = "((v1 + v3) * (v4 - v2))\n"
    assert actual_output == expected_output


# The Visitor base class leaves arguments and operators to its subclasses,
# and treats an equation as an operator.
@pytest.mark.parametrize(
    "input_method, expected_error",
    [
        # C1: on_argument is called on the base class.
        # Expected: A NotImplementedError is raised.
        ("on_argument", NotImplementedError),
        # C2: on_operator is called on the base class.
        # Expected: A NotImplementedError is raised.
        ("on_operator", NotImplementedError),
    ],
)
def test_visitor_base_bad(make_args, input_method, expected_error):
    visitor = Visitor()
    with pytest.raises(expected_error):
        getattr(visitor, input_method)(make_args(1)[0])


def test_visitor_on_equation_uses_on_operator(make_args):
    # C1: on_equation is called on a visitor that only defines
    # on_operator.
    # Expected: The equation is passed to on_operator.
    class OperatorRecorder(Visitor):
        def on_operator(self, op):
            return op

    mult, _ = _build_tree(make_args)
    equation = Equation("eq", mult)
    actual_visited = OperatorRecorder().on_equation(equation)
    assert actual_visited is equation
