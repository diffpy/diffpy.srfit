"""Tests for the deprecated camelCase names in diffpy.srfit.equation.

Each test repeats the cases of the matching test for the new name, but
calls the deprecated name and also checks that it warns.

Remove this file in version 4.0.0, together with the deprecated names.
"""

import re

import numpy
import pytest

import diffpy.srfit.equation.builder as builder
import diffpy.srfit.equation.literals as literals
import diffpy.srfit.equation.visitors as visitors
from diffpy.srfit.equation.equationmod import Equation


def expected_message(base, old_name, new_name):
    return re.escape(
        f"'{base}.{old_name}' is deprecated and will be removed in version "
        f"4.0.0. Please use '{base}.{new_name}' instead."
    )


LITERAL = "diffpy.srfit.equation.literals.Literal"
LITERALS = "diffpy.srfit.equation.literals"
OPERATOR = "diffpy.srfit.equation.literals.Operator"
EQUATION = "diffpy.srfit.equation.Equation"
EQUATION_FACTORY = "diffpy.srfit.equation.builder.EquationFactory"
BASE_BUILDER = "diffpy.srfit.equation.builder.BaseBuilder"
BUILDER = "diffpy.srfit.equation.builder"
VISITORS = "diffpy.srfit.equation.visitors"
VISITOR = "diffpy.srfit.equation.visitors.visitor.Visitor"


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


# Literal: see test_literals.py::test_get_value.
@pytest.mark.parametrize(
    "input_literal_kind",
    [
        # C1: An Argument holds the value itself.
        # Expected: Warns and returns 3.5.
        "argument",
        # C2: An Operator computes the value from its literals.
        # Expected: Warns and returns 3.5.
        "operator",
        # C3: An Equation evaluates the operator at its root.
        # Expected: Warns and returns 3.5.
        "equation",
    ],
)
def test_getValue(input_literal_kind):
    literal = _make_sum(input_literal_kind)
    with pytest.warns(
        DeprecationWarning,
        match=expected_message(LITERAL, "getValue", "get_value"),
    ):
        actual_value = literal.getValue()
    expected_value = 3.5
    assert actual_value == expected_value


# EquationFactory: see
# test_builder.py::test_register_function_replaces_existing.
@pytest.mark.parametrize(
    "input_replacement, expected_value",
    [
        # C1: g is replaced by the Python function 0.5 * v1.
        # Expected: Warns, and "g" evaluates to 0.5.
        ("python_function", 0.5),
        # C2: g is replaced by the equation "v1 + v2".
        # Expected: Warns, and "g" evaluates to 3.
        ("equation", 3),
    ],
)
def test_registerFunction(make_args, input_replacement, expected_value):
    factory = builder.EquationFactory()
    v1, v2, v3, v4 = make_args(4)
    for argument in (v1, v2, v3, v4):
        factory.register_argument(argument.name, argument)
    factory.register_function(
        "g",
        lambda v1, v2, v3, v4: (v1 + v2) * (v3 + v4),
        ["v1", "v2", "v3", "v4"],
    )
    equation = factory.make_equation("g")
    if input_replacement == "python_function":
        new_function, new_arg_names = (lambda v1: 0.5 * v1), ["v1"]
    else:
        new_function = factory.make_equation("v1 + v2")
        new_arg_names = new_function.argdict.keys()
    with pytest.warns(
        DeprecationWarning,
        match=expected_message(
            EQUATION_FACTORY, "registerFunction", "register_function"
        ),
    ):
        function_builder = factory.registerFunction(
            "g", new_function, new_arg_names
        )
    actual_root = equation.root
    assert actual_root is function_builder.literal
    actual_value = equation()
    assert actual_value == pytest.approx(expected_value)


# BaseBuilder: see test_builder.py::test_build_equation_from_builders.
def test_getEquation():
    # C1: The equation A*sin(a*x) is written with builder objects, with
    # a=1, A=2 and x an array.
    # Expected: Warns. Its arguments are A and a, and its value is
    # 2*sin(x).
    sin = builder.get_builder("sin")
    a = builder.ArgumentBuilder(name="a", value=1)
    A = builder.ArgumentBuilder(name="A", value=2)
    x = numpy.arange(0, numpy.pi, 0.1)
    with pytest.warns(
        DeprecationWarning,
        match=expected_message(BASE_BUILDER, "getEquation", "get_equation"),
    ):
        equation = (A * sin(a * x)).getEquation()
    actual_arg_names = [arg.name for arg in equation.args]
    expected_arg_names = ["A", "a"]
    assert actual_arg_names == expected_arg_names
    actual_values = equation()
    expected_values = 2 * numpy.sin(x)
    assert numpy.array_equal(actual_values, expected_values)


def _make_factory_call_args(method_name):
    """Return fresh arguments for an EquationFactory method."""
    if method_name == "make_equation":
        return ("v1 + 2",)
    if method_name == "register_constant":
        return ("c", 3)
    if method_name == "register_argument":
        return ("a", literals.Argument(name="a", value=3))
    if method_name == "register_operator":
        return ("op", _make_sum("operator"))
    if method_name == "register_builder":
        return ("b", builder.ArgumentBuilder(name="b", value=3))
    return ("v1",)


def _call_factory_method(method_name, args):
    """Call a method of a factory that has v1 = 1 registered.

    Return the names the factory has afterward and the value of what
    the method returned.
    """
    factory = builder.EquationFactory()
    factory.register_argument("v1", literals.Argument(name="v1", value=1))
    result = getattr(factory, method_name)(*args)
    if result is None:
        value = None
    elif isinstance(result, Equation):
        value = result()
    else:
        value = result.literal.value
    return set(factory.builders), value


# EquationFactory: each deprecated method does the same as its new name,
# on a factory that has v1 = 1 registered. See test_builder.py.
@pytest.mark.parametrize(
    "input_old_name, input_new_name",
    [
        # C1: makeEquation builds "v1 + 2".
        # Expected: Warns and gives the same result as make_equation.
        ("makeEquation", "make_equation"),
        # C2: registerConstant registers c = 3.
        # Expected: Warns and gives the same result as register_constant.
        ("registerConstant", "register_constant"),
        # C3: registerArgument registers the argument a = 3.
        # Expected: Warns and gives the same result as register_argument.
        ("registerArgument", "register_argument"),
        # C4: registerOperator registers the operator op = 1.5 + 2.0.
        # Expected: Warns and gives the same result as register_operator.
        ("registerOperator", "register_operator"),
        # C5: registerBuilder registers a builder for b = 3.
        # Expected: Warns and gives the same result as register_builder.
        ("registerBuilder", "register_builder"),
        # C6: deRegisterBuilder removes v1.
        # Expected: Warns and gives the same result as deregister_builder.
        ("deRegisterBuilder", "deregister_builder"),
    ],
)
def test_equation_factory_deprecations(input_old_name, input_new_name):
    with pytest.warns(
        DeprecationWarning,
        match=expected_message(
            EQUATION_FACTORY, input_old_name, input_new_name
        ),
    ):
        actual_old_result = _call_factory_method(
            input_old_name, _make_factory_call_args(input_new_name)
        )
    actual_new_result = _call_factory_method(
        input_new_name, _make_factory_call_args(input_new_name)
    )
    assert actual_old_result == actual_new_result


def _make_builder_call_args(function_name):
    """Return fresh arguments for a builder module function."""
    if function_name == "wrap_argument":
        return ("a", literals.Argument(name="a", value=3))
    if function_name == "wrap_operator":
        return ("op", _make_sum("operator"))
    if function_name == "wrap_function":
        return ("f", numpy.add, 2, 1)
    return ("sin",)


def _describe_builder(wrapped_builder):
    literal = wrapped_builder.literal
    return type(wrapped_builder), type(literal), literal.name


# builder module: each deprecated function does the same as its new name.
# See test_builder.py.
@pytest.mark.parametrize(
    "input_old_name, input_new_name",
    [
        # C1: wrapArgument wraps the argument a = 3.
        # Expected: Warns and gives the same builder as wrap_argument.
        ("wrapArgument", "wrap_argument"),
        # C2: wrapOperator wraps the operator op = 1.5 + 2.0.
        # Expected: Warns and gives the same builder as wrap_operator.
        ("wrapOperator", "wrap_operator"),
        # C3: wrapFunction wraps numpy.add with two inputs.
        # Expected: Warns and gives the same builder as wrap_function.
        ("wrapFunction", "wrap_function"),
        # C4: getBuilder gets the global builder for sin.
        # Expected: Warns and gives the same builder as get_builder.
        ("getBuilder", "get_builder"),
    ],
)
def test_builder_module_deprecations(input_old_name, input_new_name):
    with pytest.warns(
        DeprecationWarning,
        match=expected_message(BUILDER, input_old_name, input_new_name),
    ):
        actual_old_builder = getattr(builder, input_old_name)(
            *_make_builder_call_args(input_new_name)
        )
    actual_new_builder = getattr(builder, input_new_name)(
        *_make_builder_call_args(input_new_name)
    )
    actual_old_result = _describe_builder(actual_old_builder)
    actual_new_result = _describe_builder(actual_new_builder)
    assert actual_old_result == actual_new_result


def _set_equation_root(method_name):
    equation = Equation(name="eq", root=_make_sum("operator"))
    getattr(equation, method_name)(literals.Argument(name="v", value=5))
    return [arg.name for arg in equation.args], equation()


# Equation: see test_equation.py::test_equation_set_root.
def test_setRoot():
    # C1: An equation with root 1.5 + 2.0 has its root set to v = 5.
    # Expected: Warns, and the equation has the same arguments and value
    # as after set_root.
    with pytest.warns(
        DeprecationWarning,
        match=expected_message(EQUATION, "setRoot", "set_root"),
    ):
        actual_old_result = _set_equation_root("setRoot")
    actual_new_result = _set_equation_root("set_root")
    assert actual_old_result == actual_new_result


def _add_literals(method_name):
    operator = literals.AdditionOperator()
    getattr(operator, method_name)(literals.Argument(name="a", value=1))
    getattr(operator, method_name)(literals.Argument(name="b", value=2))
    return [arg.name for arg in operator.args], operator.value


# Operator: see test_literals.py::test_operator_value_after_argument_change.
def test_addLiteral():
    # C1: The arguments a = 1 and b = 2 are added to an add operator.
    # Expected: Warns, and the operator has the same arguments and value
    # as after add_literal.
    with pytest.warns(
        DeprecationWarning,
        match=expected_message(OPERATOR, "addLiteral", "add_literal"),
    ):
        actual_old_result = _add_literals("addLiteral")
    actual_new_result = _add_literals("add_literal")
    assert actual_old_result == actual_new_result


# Equation: see test_equation.py::test_equation_add_literal_bad.
def test_addLiteral_equation_bad():
    # C1: A literal is added to an equation.
    # Expected: Warns and raises a RuntimeError.
    equation = Equation(name="eq", root=_make_sum("operator"))
    with (
        pytest.warns(
            DeprecationWarning,
            match=expected_message(OPERATOR, "addLiteral", "add_literal"),
        ),
        pytest.raises(RuntimeError),
    ):
        equation.addLiteral(literals.Argument(name="c", value=3))


def _describe_operator(operator):
    return (
        type(operator),
        operator.name,
        operator.symbol,
        operator.operation,
        operator.nin,
        operator.nout,
    )


# literals module: see test_literals.py::test_operator_init.
def test_makeOperator():
    # C1: An add operator is made.
    # Expected: Warns and gives the same operator as make_operator.
    operator_kwargs = dict(
        name="add", symbol="+", operation=numpy.add, nin=2, nout=1
    )
    with pytest.warns(
        DeprecationWarning,
        match=expected_message(LITERALS, "makeOperator", "make_operator"),
    ):
        actual_old_operator = literals.makeOperator(**operator_kwargs)
    actual_new_operator = literals.make_operator(**operator_kwargs)
    actual_old_result = _describe_operator(actual_old_operator)
    actual_new_result = _describe_operator(actual_new_operator)
    assert actual_old_result == actual_new_result


def _make_embedded_tree():
    """Return the tree inner + c, where the equation inner = a + b."""
    top = literals.AdditionOperator()
    top.add_literal(Equation(name="inner", root=_make_sum("operator")))
    top.add_literal(literals.Argument(name="c", value=1.0, const=True))
    return top


def _describe_args(args):
    return [arg.name for arg in args]


# visitors module: each deprecated function does the same as its new name
# on the tree inner + c. See test_visitors.py.
@pytest.mark.parametrize(
    "input_old_name, input_new_name, input_kwargs, input_describe",
    [
        # C1: getArgs gets the arguments of the tree.
        # Expected: Warns and gives the same arguments as get_args.
        ("getArgs", "get_args", {}, _describe_args),
        # C2: getArgs gets the non-constant arguments of the tree.
        # Expected: Warns and gives the same arguments as get_args.
        ("getArgs", "get_args", {"getconsts": False}, _describe_args),
        # C3: getExpression writes the tree out in full.
        # Expected: Warns and gives the same string as get_expression.
        ("getExpression", "get_expression", {}, str),
        # C4: getExpression writes the tree with "inner" by name.
        # Expected: Warns and gives the same string as get_expression.
        ("getExpression", "get_expression", {"eqskip": "inner"}, str),
    ],
)
def test_visitors_module_deprecations(
    input_old_name, input_new_name, input_kwargs, input_describe
):
    with pytest.warns(
        DeprecationWarning,
        match=expected_message(VISITORS, input_old_name, input_new_name),
    ):
        actual_old_result = input_describe(
            getattr(visitors, input_old_name)(
                _make_embedded_tree(), **input_kwargs
            )
        )
    actual_new_result = input_describe(
        getattr(visitors, input_new_name)(
            _make_embedded_tree(), **input_kwargs
        )
    )
    assert actual_old_result == actual_new_result


# visitors module: see test_visitors.py::test_pretty_print.
def test_prettyPrint(capsys):
    # C1: The tree inner + c is printed.
    # Expected: Warns and prints the same line as pretty_print.
    with pytest.warns(
        DeprecationWarning,
        match=expected_message(VISITORS, "prettyPrint", "pretty_print"),
    ):
        visitors.prettyPrint(_make_embedded_tree())
    actual_old_output = capsys.readouterr().out
    visitors.pretty_print(_make_embedded_tree())
    actual_new_output = capsys.readouterr().out
    assert actual_old_output == actual_new_output


def _make_visited_literal(method_name):
    if method_name == "on_argument":
        return literals.Argument(name="a", value=1.5)
    if method_name == "on_operator":
        return _make_sum("operator")
    return Equation(name="eq", root=_make_sum("operator"))


# Visitor: each deprecated method on the Visitor base forwards to its new
# name in a subclass. A Printer writes the visited literal.
@pytest.mark.parametrize(
    "input_old_name, input_new_name",
    [
        # C1: onArgument writes the argument a.
        # Expected: Warns and writes the same as on_argument.
        ("onArgument", "on_argument"),
        # C2: onOperator writes the operator a + b.
        # Expected: Warns and writes the same as on_operator.
        ("onOperator", "on_operator"),
        # C3: onEquation writes the equation eq = a + b.
        # Expected: Warns and writes the same as on_equation.
        ("onEquation", "on_equation"),
    ],
)
def test_visitor_deprecations(input_old_name, input_new_name):
    with pytest.warns(
        DeprecationWarning,
        match=expected_message(VISITOR, input_old_name, input_new_name),
    ):
        actual_old_output = getattr(visitors.Printer(), input_old_name)(
            _make_visited_literal(input_new_name)
        )
    actual_new_output = getattr(visitors.Printer(), input_new_name)(
        _make_visited_literal(input_new_name)
    )
    assert actual_old_output == actual_new_output
