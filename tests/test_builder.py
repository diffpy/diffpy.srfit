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
"""Tests for the builder module."""

import numpy
import pytest

import diffpy.srfit.equation.builder as builder
import diffpy.srfit.equation.literals as literals


def _gaussian(x, sigma):
    return numpy.sqrt(numpy.exp(-0.5 * (x / sigma) ** 2))


def _halve(v1):
    return 0.5 * v1


def _make_factory_with_function(make_args):
    """Return a factory where g(v1, v2, v3, v4) = (v1 + v2) * (v3 + v4).

    The arguments v1 to v4 have the values 1 to 4, so g is 21. The fifth
    argument from make_args, v5 = 5, is returned unregistered.
    """
    factory = builder.EquationFactory()
    v1, v2, v3, v4, v5 = make_args(5)
    for argument in (v1, v2, v3, v4):
        factory.register_argument(argument.name, argument)

    def g(v1, v2, v3, v4):
        return (v1 + v2) * (v3 + v4)

    function_builder = factory.register_function(
        "g", g, ["v1", "v2", "v3", "v4"]
    )
    return factory, function_builder.literal, (v1, v2, v3, v4, v5), g


def test_register_argument(make_args):
    # C1: An argument is registered under the name "v1".
    # Expected: The returned builder is stored under that name and wraps
    # the argument, and an equation "v1" has the argument as its only
    # argument.
    factory = builder.EquationFactory()
    v1 = make_args(1)[0]
    argument_builder = factory.register_argument("v1", v1)
    actual_builder = factory.builders["v1"]
    expected_builder = argument_builder
    assert actual_builder == expected_builder
    actual_literal = argument_builder.literal
    expected_literal = v1
    assert actual_literal == expected_literal
    actual_args = factory.make_equation("v1").args
    expected_args = [v1]
    assert actual_args == expected_args


# deregister_builder removes a name from the factory.
@pytest.mark.parametrize(
    "input_name, expected_names_removed",
    [
        # C1: The registered argument v1 is deregistered.
        # Expected: v1 is no longer a builder of the factory.
        ("v1", {"v1"}),
        # C2: A name that was never registered is deregistered.
        # Expected: Nothing happens, and no error is raised.
        ("missing", set()),
    ],
)
def test_deregister_builder(make_args, input_name, expected_names_removed):
    factory = builder.EquationFactory()
    v1 = make_args(1)[0]
    factory.register_argument("v1", v1)
    names_before = set(factory.builders)
    factory.deregister_builder(input_name)
    actual_names_removed = names_before - set(factory.builders)
    assert actual_names_removed == expected_names_removed


# With buildargs=False, the factory builds from registered arguments only.
def test_make_equation_without_building_args(make_args):
    factory = builder.EquationFactory()
    v1 = make_args(1)[0]
    factory.register_argument("v1", v1)
    actual_args = factory.make_equation("v1 + 2", False).args
    expected_args = [v1]
    assert actual_args == expected_args


# With buildargs=False, the factory may not create new arguments, so an
# equation using an unregistered name raises a ValueError.
def test_make_equation_without_building_args_bad(make_args):
    factory = builder.EquationFactory()
    v1 = make_args(1)[0]
    factory.register_argument("v1", v1)
    with pytest.raises(ValueError):
        factory.make_equation("v1 + v2", False)


# A registered operator op = v1 + v2 (= 3) can be used as a value or called
# as a function. The arguments v3 = 3 and v4 = 4 are also registered.
@pytest.mark.parametrize(
    "input_equations, expected_value",
    [
        # C1: "op" on its own.
        # Expected: The operator's value, 3.
        (["op"], 3),
        # C2: "v3*op".
        # Expected: 3 * 3 = 9.
        (["v3*op"], 9),
        # C3: "op(v3, v4)" calls the operator with new arguments.
        # Expected: 3 + 4 = 7.
        (["op(v3, v4)"], 7),
        # C4: "op(v3, v4)" is built first, then "op" on its own.
        # Expected: "op" still uses its own arguments, so 3.
        (["op(v3, v4)", "op"], 3),
    ],
)
def test_register_operator(make_args, input_equations, expected_value):
    factory = builder.EquationFactory()
    v1, v2, v3, v4 = make_args(4)
    operator = literals.AdditionOperator()
    operator.add_literal(v1)
    operator.add_literal(v2)
    factory.register_argument("v3", v3)
    factory.register_argument("v4", v4)
    factory.register_operator("op", operator)
    for equation_string in input_equations:
        equation = factory.make_equation(equation_string)
    actual_value = equation()
    assert actual_value == expected_value


def test_register_function(make_args):
    # C1: A Python function g is registered with the arguments v1 to v4.
    # Expected: The builder's operator wraps g, takes v1 to v4 and has the
    # value 21, and an equation "g" has that operator as its root.
    factory, operator, (v1, v2, v3, v4, _), g = _make_factory_with_function(
        make_args
    )
    actual_operation = operator.operation
    expected_operation = g
    assert actual_operation == expected_operation
    actual_args = operator.args
    expected_args = [v1, v2, v3, v4]
    assert actual_args == expected_args
    actual_operator_value = operator.value
    expected_value = 21
    assert actual_operator_value == expected_value
    equation = factory.make_equation("g")
    actual_root = equation.root
    expected_root = operator
    assert actual_root == expected_root
    actual_equation_value = equation()
    assert actual_equation_value == expected_value


def test_register_argument_replaces_existing(make_args):
    # C1: v5 is registered under the name "v4", which g already uses.
    # Expected: The new builder replaces the old one, g's cached value is
    # cleared, g uses v5 instead of v4, and "g" evaluates to
    # (1 + 2) * (3 + 5) = 24.
    factory, operator, (v1, v2, v3, _, v5), _ = _make_factory_with_function(
        make_args
    )
    equation = factory.make_equation("g")
    argument_builder = factory.register_argument("v4", v5)
    actual_builder = factory.builders["v4"]
    expected_builder = argument_builder
    assert actual_builder == expected_builder
    actual_literal = argument_builder.literal
    expected_literal = v5
    assert actual_literal == expected_literal
    actual_cached_value = operator._value
    expected_cached_value = None
    assert actual_cached_value == expected_cached_value
    actual_args = operator.args
    expected_args = [v1, v2, v3, v5]
    assert actual_args == expected_args
    actual_value = equation()
    expected_value = 24
    assert actual_value == expected_value


# Registering a new function under the name "g" replaces g inside equations
# that already use it.
@pytest.mark.parametrize(
    "input_replacement, expected_value",
    [
        # C1: g is replaced by the Python function 0.5 * v1.
        # Expected: "g" evaluates to 0.5.
        ("python_function", 0.5),
        # C2: g is replaced by the equation "v1 + v2".
        # Expected: "g" evaluates to 3.
        ("equation", 3),
    ],
)
def test_register_function_replaces_existing(
    make_args, input_replacement, expected_value
):
    factory, _, (v1, v2, _, _, _), _ = _make_factory_with_function(make_args)
    equation = factory.make_equation("g")
    if input_replacement == "python_function":
        function_builder = factory.register_function("g", _halve, ["v1"])
        expected_args = [v1]
    else:
        replacement = factory.make_equation("v1 + v2")
        function_builder = factory.register_function(
            "g", replacement, replacement.argdict.keys()
        )
        expected_args = [v1, v2]
    new_operator = function_builder.literal
    actual_root = equation.root
    expected_root = new_operator
    assert actual_root == expected_root
    actual_args = new_operator.args
    assert actual_args == expected_args
    actual_value = equation()
    assert actual_value == expected_value


# The factory creates the arguments of an equation string in the order they
# first appear.
@pytest.mark.parametrize(
    "input_equation, expected_arg_names",
    [
        # C1: A scalar equation with four arguments.
        # Expected: A, x, B, C.
        ("A*sin(0.5*x)+divide(B,C)", ["A", "x", "B", "C"]),
        # C2: A Gaussian that uses the constant e.
        # Expected: x, sigma, without e.
        ("sqrt(e**(-0.5*(x/sigma)**2))", ["x", "sigma"]),
    ],
)
def test_make_equation_arg_order(input_equation, expected_arg_names):
    factory = builder.EquationFactory()
    equation = factory.make_equation(input_equation)
    actual_arg_names = [arg.name for arg in equation.args]
    assert actual_arg_names == expected_arg_names


def test_make_equation_scalar():
    # C1: "A*sin(0.5*x)+divide(B,C)" with A=1, x=pi, B=4 and C=2.
    # Expected: 1*sin(pi/2) + 4/2 = 3.
    factory = builder.EquationFactory()
    equation = factory.make_equation("A*sin(0.5*x)+divide(B,C)")
    equation.A.set_value(1)
    equation.x.set_value(numpy.pi)
    equation.B.set_value(4.0)
    equation.C.set_value(2.0)
    actual_value = equation()
    expected_value = 3
    assert actual_value == expected_value


def test_make_equation_array():
    # C1: A Gaussian is evaluated on an array of x values.
    # Expected: The Gaussian at each x.
    factory = builder.EquationFactory()
    equation = factory.make_equation("sqrt(e**(-0.5*(x/sigma)**2))")
    x = numpy.arange(0, 1, 0.05)
    equation.x.set_value(x)
    equation.sigma.set_value(0.1)
    actual_values = equation()
    expected_values = _gaussian(x, 0.1)
    assert numpy.allclose(actual_values, expected_values)


def test_make_equation_with_registered_constant():
    # C1: x is registered as a constant before the Gaussian is built.
    # Expected: sigma is the only argument, and the Gaussian uses the
    # constant x.
    factory = builder.EquationFactory()
    x = numpy.arange(0, 1, 0.05)
    factory.register_constant("x", x)
    equation = factory.make_equation("sqrt(e**(-0.5*(x/sigma)**2))")
    actual_arg_names = [arg.name for arg in equation.args]
    expected_arg_names = ["sigma"]
    assert actual_arg_names == expected_arg_names
    actual_values = equation(sigma=0.1)
    expected_values = _gaussian(x, 0.1)
    assert numpy.allclose(actual_values, expected_values)


def test_make_equation_with_registered_equation():
    # C1: A Gaussian equation is registered as the function myfunc(sigma),
    # then used in "c*myfunc(sigma)".
    # Expected: The arguments are c and sigma, and the value is c times the
    # Gaussian.
    factory = builder.EquationFactory()
    x = numpy.arange(0, 1, 0.05)
    factory.register_constant("x", x)
    gaussian = factory.make_equation("sqrt(e**(-0.5*(x/sigma)**2))")
    factory.register_function("myfunc", gaussian, ["sigma"])
    equation = factory.make_equation("c*myfunc(sigma)")
    actual_arg_names = [arg.name for arg in equation.args]
    expected_arg_names = ["c", "sigma"]
    assert actual_arg_names == expected_arg_names
    actual_values = equation(c=2, sigma=0.1)
    expected_values = 2 * _gaussian(x, 0.1)
    assert numpy.allclose(actual_values, expected_values)


def test_make_equation_constant_expression():
    # C1: An equation string with only numbers, "3.12 + 2".
    # Expected: An Equation with the value 5.12, which the factory does not
    # keep, and which takes no arguments when called.
    factory = builder.EquationFactory()
    equation = factory.make_equation("3.12 + 2")
    actual_type = type(equation)
    expected_type = builder.Equation
    assert actual_type == expected_type
    actual_kept_equations = factory.equations
    expected_kept_equations = set()
    assert actual_kept_equations == expected_kept_equations
    actual_value = equation()
    expected_value = 5.12
    assert actual_value == expected_value
    with pytest.raises(ValueError):
        equation(3)


# Global builders live for the whole session. Building equations must not
# leave observers on them, or those equations would never be freed.
@pytest.mark.parametrize(
    "input_equation",
    [
        # C1: An equation that uses numpy functions.
        # Expected: No global builder has observers afterwards.
        "A*sin(0.5*x)+divide(B,C)",
        # C2: An equation that uses the constant e.
        # Expected: No global builder has observers afterwards.
        "sqrt(e**(-0.5*(x/sigma)**2))",
        # C3: An equation with only numbers.
        # Expected: No global builder has observers afterwards.
        "3.12 + 2",
    ],
)
def test_make_equation_leaves_no_global_observers(
    find_builders_with_observers, input_equation
):
    factory = builder.EquationFactory()
    factory.make_equation(input_equation)
    actual_builders = find_builders_with_observers()
    expected_builders = []
    assert actual_builders == expected_builders


def test_build_equation_from_builders():
    # C1: The equation A*sin(a*x) is written with builder objects, with
    # a=1, A=2 and x an array.
    # Expected: Its arguments are A and a, and its value is 2*sin(x).
    sin = builder.get_builder("sin")
    a = builder.ArgumentBuilder(name="a", value=1)
    A = builder.ArgumentBuilder(name="A", value=2)
    x = numpy.arange(0, numpy.pi, 0.1)
    equation = (A * sin(a * x)).get_equation()
    actual_arg_names = [arg.name for arg in equation.args]
    expected_arg_names = ["A", "a"]
    assert actual_arg_names == expected_arg_names
    actual_values = equation()
    expected_values = 2 * numpy.sin(x)
    assert numpy.array_equal(actual_values, expected_values)


def test_build_equation_with_wrapped_function():
    # C1: A Python function f(a, b) = (a - b) / (a + b) is wrapped and used
    # in sin(f(a, b)) with a=2 and b=1.
    # Expected: sin(f(2, 1)).
    def f(a, b):
        return (a - b) * 1.0 / (a + b)

    wrapped_f = builder.wrap_function("f", f, 2, 1)
    sin = builder.get_builder("sin")
    a = builder.ArgumentBuilder(name="a", value=2)
    b = builder.ArgumentBuilder(name="b", value=1)
    equation = sin(wrapped_f(a, b)).get_equation()
    actual_value = equation()
    expected_value = numpy.sin(f(2, 1))
    assert actual_value == expected_value


def test_build_equation_with_constant_argument():
    # C1: A Gaussian is written with builders, with x a constant array.
    # Expected: The Gaussian at each x.
    sqrt = builder.get_builder("sqrt")
    x_values = numpy.arange(0, 1, 0.05)
    x = builder.ArgumentBuilder(name="x", value=x_values, const=True)
    sigma = builder.ArgumentBuilder(name="sigma", value=0.1)
    equation = sqrt(numpy.e ** (-0.5 * (x / sigma) ** 2)).get_equation()
    actual_values = equation()
    expected_values = _gaussian(x_values, 0.1)
    assert numpy.allclose(actual_values, expected_values)


# A builder for a numpy function checks how many arguments it is given.
@pytest.mark.parametrize(
    "input_arg_count",
    [
        # C1: sin is called with no arguments.
        # Expected: A ValueError is raised.
        0,
        # C2: sin is called with two arguments.
        # Expected: A ValueError is raised.
        2,
    ],
)
def test_builder_wrong_arg_count_bad(input_arg_count):
    sin = builder.get_builder("sin")
    args = [builder.ArgumentBuilder(value=1)] * input_arg_count
    with pytest.raises(ValueError):
        sin(*args)


# The equation A + B, with A=2 and B=4, is wrapped as the operator E.
@pytest.mark.parametrize(
    "input_expression, expected_value",
    [
        # C1: E is used as a value, in 2*E.
        # Expected: 2 * (2 + 4) = 12.
        ("2*E", 12),
        # C2: E is called with the new arguments C=5 and D=6, in E(C, D) + 1.
        # Expected: (5 + 6) + 1 = 12.
        ("E(C, D) + 1", 12),
        # C3: E is called with its own A and the new D=6, in 3*E(A, D) - 1.
        # Expected: 3 * (2 + 6) - 1 = 23.
        ("3*E(A, D) - 1", 23),
    ],
)
def test_build_equation_with_wrapped_equation(
    input_expression, expected_value
):
    A = builder.ArgumentBuilder(name="A", value=2)
    B = builder.ArgumentBuilder(name="B", value=4)
    C = builder.ArgumentBuilder(name="C", value=5)
    D = builder.ArgumentBuilder(name="D", value=6)
    E = builder.wrap_operator("eq", (A + B).get_equation())
    names = {"A": A, "B": B, "C": C, "D": D, "E": E}
    equation = eval(input_expression, {}, names).get_equation()
    actual_value = equation()
    assert actual_value == expected_value


def test_build_equation_wrapped_equation_copies_argument_values():
    # C1: The equation A + B, with A=2 and B=4, is wrapped as E and
    # evaluated as E(C, D) with C=5 and D=6.
    # Expected: The values of C and D are copied into A and B.
    A = builder.ArgumentBuilder(name="A", value=2)
    B = builder.ArgumentBuilder(name="B", value=4)
    C = builder.ArgumentBuilder(name="C", value=5)
    D = builder.ArgumentBuilder(name="D", value=6)
    E = builder.wrap_operator("eq", (A + B).get_equation())
    E(C, D).get_equation()()
    actual_values = (A.literal.value, B.literal.value)
    expected_values = (5, 6)
    assert actual_values == expected_values


# A wrapped equation must be called with as many arguments as it has.
@pytest.mark.parametrize(
    "input_arg_names",
    [
        # C1: E(A + B) is called with one argument.
        # Expected: A ValueError is raised.
        ["A"],
        # C2: E(A + B) is called with three arguments.
        # Expected: A ValueError is raised.
        ["A", "B", "C"],
    ],
)
def test_build_equation_wrapped_equation_wrong_arg_count_bad(input_arg_names):
    arguments = {
        name: builder.ArgumentBuilder(name=name, value=1)
        for name in ("A", "B", "C")
    }
    E = builder.wrap_operator(
        "eq", (arguments["A"] + arguments["B"]).get_equation()
    )
    with pytest.raises(ValueError):
        E(*[arguments[name] for name in input_arg_names])
