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
"""Tests for refinableobj module."""

import numpy
import pytest

import diffpy.srfit.equation.builder as builder
import diffpy.srfit.equation.literals as literals


def testRegisterArg(make_args, noObserversInGlobalBuilders):

    factory = builder.EquationFactory()

    v1 = make_args(1)[0]

    b1 = factory.registerArgument("v1", v1)
    assert factory.builders["v1"] is b1
    assert b1.literal is v1

    eq = factory.makeEquation("v1")

    assert v1 is eq.args[0]
    assert 1 == len(eq.args)

    # Try to parse an equation with buildargs turned off
    with pytest.raises(ValueError):
        factory.makeEquation("v1 + v2", False)

    # Make sure we can still use constants
    eq = factory.makeEquation("v1 + 2", False)
    assert v1 is eq.args[0]
    assert 1 == len(eq.args)
    assert noObserversInGlobalBuilders
    return


def testRegisterOperator(make_args, noObserversInGlobalBuilders):
    """Try to use an operator without arguments in an equation."""

    factory = builder.EquationFactory()
    v1, v2, v3, v4 = make_args(4)

    op = literals.AdditionOperator()

    op.addLiteral(v1)
    op.addLiteral(v2)

    factory.registerArgument("v3", v3)
    factory.registerArgument("v4", v4)
    factory.registerOperator("op", op)

    # Build an equation where op is treated as a terminal node
    eq = factory.makeEquation("op")
    assert 3 == pytest.approx(eq())

    eq = factory.makeEquation("v3*op")
    assert 9 == pytest.approx(eq())

    # Now use the op like a function
    eq = factory.makeEquation("op(v3, v4)")
    assert 7 == pytest.approx(eq())

    # Make sure we can still access op as itself.
    eq = factory.makeEquation("op")
    assert 3 == pytest.approx(eq())

    assert noObserversInGlobalBuilders
    return


def testSwapping(make_args, noObserversInGlobalBuilders):

    def g1(v1, v2, v3, v4):
        return (v1 + v2) * (v3 + v4)

    def g2(v1):
        return 0.5 * v1

    factory = builder.EquationFactory()
    v1, v2, v3, v4, v5 = make_args(5)

    factory.registerArgument("v1", v1)
    factory.registerArgument("v2", v2)
    factory.registerArgument("v3", v3)
    factory.registerArgument("v4", v4)
    b = factory.registerFunction("g", g1, ["v1", "v2", "v3", "v4"])

    # Now associate args with the wrapped function
    op = b.literal
    assert op.operation == g1
    assert v1 in op.args
    assert v2 in op.args
    assert v3 in op.args
    assert v4 in op.args
    assert round(abs(21 - op.value), 7) == 0

    eq1 = factory.makeEquation("g")
    assert eq1.root is op
    assert round(abs(21 - eq1()), 7) == 0

    # Swap out an argument by registering it under a taken name
    b = factory.registerArgument("v4", v5)
    assert factory.builders["v4"] is b
    assert b.literal is v5
    assert op._value is None
    assert op.args == [v1, v2, v3, v5]
    assert round(abs(24 - eq1()), 7) == 0

    # Now swap out the function
    b = factory.registerFunction("g", g2, ["v1"])
    op = b.literal
    assert op.operation == g2
    assert v1 in op.args
    assert eq1.root is op
    assert round(abs(0.5 - op.value), 7) == 0
    assert round(abs(0.5 - eq1()), 7) == 0

    # Make an equation
    eqeq = factory.makeEquation("v1 + v2")
    # Register this "g"
    b = factory.registerFunction("g", eqeq, eqeq.argdict.keys())
    op = b.literal
    assert v1 in op.args
    assert v2 in op.args
    assert eq1.root is op
    assert round(abs(3 - op.value), 7) == 0
    assert round(abs(3 - eq1()), 7) == 0
    assert noObserversInGlobalBuilders
    return


def testParseEquation(noObserversInGlobalBuilders):

    factory = builder.EquationFactory()

    # Scalar equation
    eq = factory.makeEquation("A*sin(0.5*x)+divide(B,C)")
    A = 1
    x = numpy.pi
    B = 4.0
    C = 2.0
    eq.A.setValue(A)
    eq.x.setValue(x)
    eq.B.setValue(B)
    eq.C.setValue(C)
    assert numpy.array_equal(eq(), f_equation(A, x, B, C))

    # Make sure that the arguments of eq are listed in the order in which
    # they appear in the equations.
    assert eq.args == [eq.A, eq.x, eq.B, eq.C]

    # Vector equation
    eq = factory.makeEquation("sqrt(e**(-0.5*(x/sigma)**2))")
    x = numpy.arange(0, 1, 0.05)
    sigma = 0.1
    eq.x.setValue(x)
    eq.sigma.setValue(sigma)
    assert numpy.allclose(eq(), gaussian_test(x, sigma))
    assert eq.args == [eq.x, eq.sigma]

    # Equation with constants
    factory.registerConstant("x", x)
    eq = factory.makeEquation("sqrt(e**(-0.5*(x/sigma)**2))")
    assert "sigma" in eq.argdict
    assert "x" not in eq.argdict
    assert numpy.allclose(eq(sigma=sigma), gaussian_test(x, sigma))
    assert eq.args == [eq.sigma]

    # Equation with user-defined functions
    factory.registerFunction("myfunc", eq, ["sigma"])
    eq2 = factory.makeEquation("c*myfunc(sigma)")
    assert numpy.allclose(eq2(c=2, sigma=sigma), 2 * gaussian_test(x, sigma))
    assert "sigma" in eq2.argdict
    assert "c" in eq2.argdict
    assert eq2.args == [eq2.c, eq2.sigma]
    assert noObserversInGlobalBuilders
    return


def test_parse_constant():
    """Verify parsing of constant numeric expressions."""
    factory = builder.EquationFactory()
    eq = factory.makeEquation("3.12 + 2")
    assert isinstance(eq, builder.Equation)
    assert set() == factory.equations
    assert 5.12 == eq()
    with pytest.raises(ValueError):
        eq(3)
    return


def testBuildEquation(noObserversInGlobalBuilders):

    from numpy import array_equal

    # simple equation
    sin = builder.getBuilder("sin")
    a = builder.ArgumentBuilder(name="a", value=1)
    A = builder.ArgumentBuilder(name="A", value=2)
    x = numpy.arange(0, numpy.pi, 0.1)

    beq = A * sin(a * x)
    eq = beq.getEquation()

    assert "a" in eq.argdict
    assert "A" in eq.argdict
    assert array_equal(eq(), 2 * numpy.sin(x))

    assert eq.args == [eq.A, eq.a]

    # Check the number of arguments
    with pytest.raises(ValueError):
        sin()

    # custom function
    def _f(a, b):
        return (a - b) * 1.0 / (a + b)

    f = builder.wrapFunction("f", _f, 2, 1)
    a = builder.ArgumentBuilder(name="a", value=2)
    b = builder.ArgumentBuilder(name="b", value=1)

    beq = sin(f(a, b))
    eq = beq.getEquation()
    assert eq() == numpy.sin(_f(2, 1))

    # complex function
    sqrt = builder.getBuilder("sqrt")
    e = numpy.e
    _x = numpy.arange(0, 1, 0.05)
    x = builder.ArgumentBuilder(name="x", value=_x, const=True)
    sigma = builder.ArgumentBuilder(name="sigma", value=0.1)
    beq = sqrt(e ** (-0.5 * (x / sigma) ** 2))
    eq = beq.getEquation()
    assert numpy.allclose(eq(), numpy.sqrt(numpy.exp(-0.5 * (_x / 0.1) ** 2)))

    # Equation with Equation
    A = builder.ArgumentBuilder(name="A", value=2)
    B = builder.ArgumentBuilder(name="B", value=4)
    beq = A + B
    eq = beq.getEquation()
    E = builder.wrapOperator("eq", eq)
    eq2 = (2 * E).getEquation()
    # Make sure these evaluate to the same thing
    assert eq.args == [A.literal, B.literal]
    assert 2 * eq() == eq2()
    # Pass new arguments to the equation
    C = builder.ArgumentBuilder(name="C", value=5)
    D = builder.ArgumentBuilder(name="D", value=6)
    eq3 = (E(C, D) + 1).getEquation()
    assert 12 == eq3()
    # Pass old and new arguments to the equation
    # If things work right, A has been given the value of C in the last
    # evaluation (5)
    eq4 = (3 * E(A, D) - 1).getEquation()
    assert 32 == eq4()
    # Try to pass the wrong number of arguments
    with pytest.raises(ValueError):
        E(A)
    with pytest.raises(ValueError):
        E(A, B, C)
    assert noObserversInGlobalBuilders
    return


def f_equation(a, x, b, c):
    return a * numpy.sin(0.5 * x) + numpy.divide(b, c)


def gaussian_test(x, sigma):
    return numpy.sqrt(numpy.exp(-0.5 * (x / sigma) ** 2))
