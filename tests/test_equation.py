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

import pytest

import diffpy.srfit.equation.literals as literals
from diffpy.srfit.equation import Equation


def testSimpleFunction(make_args, noObserversInGlobalBuilders):
    """Test a simple function."""

    # Make some variables
    v1, v2, v3, v4, c = make_args(5)
    c.name = "c"
    c.const = True

    # Make some operations
    mult = literals.MultiplicationOperator()
    root = mult2 = literals.MultiplicationOperator()
    plus = literals.AdditionOperator()
    minus = literals.SubtractionOperator()

    # Create the equation c*(v1+v3)*(v4-v2)
    plus.addLiteral(v1)
    plus.addLiteral(v3)
    minus.addLiteral(v4)
    minus.addLiteral(v2)
    mult.addLiteral(plus)
    mult.addLiteral(minus)
    mult2.addLiteral(mult)
    mult2.addLiteral(c)

    # Set the values of the variables.
    # The equation should evaluate to 2.5*(1+3)*(4-2) = 20
    v1.setValue(1)
    v2.setValue(2)
    v3.setValue(3)
    v4.setValue(4)
    c.setValue(2.5)

    # Make an equation and test
    eq = Equation("eq", mult2)

    assert eq._value is None
    args = eq.args
    assert v1 in args
    assert v2 in args
    assert v3 in args
    assert v4 in args
    assert c not in args
    assert root is eq.root

    assert v1 is eq.v1
    assert v2 is eq.v2
    assert v3 is eq.v3
    assert v4 is eq.v4

    assert 20 == eq()  # 20 = 2.5*(1+3)*(4-2)
    assert 20 == eq.getValue()  # same as above
    assert 20 == eq.value  # same as above
    assert 25 == eq(v1=2)  # 25 = 2.5*(2+3)*(4-2)
    assert 50 == eq(v2=0)  # 50 = 2.5*(2+3)*(4-0)
    assert 30 == eq(v3=1)  # 30 = 2.5*(2+1)*(4-0)
    assert 0 == eq(v4=0)  # 20 = 2.5*(2+1)*(0-0)

    # Try some swapping
    eq.swap(v4, v1)
    assert eq._value is None
    assert 15 == eq()  # 15 = 2.5*(2+1)*(2-0)
    args = eq.args
    assert v4 not in args

    # Try to create a dependency loop
    with pytest.raises(ValueError):
        eq.swap(v1, eq.root)

    with pytest.raises(ValueError):
        eq.swap(v1, plus)

    with pytest.raises(ValueError):
        eq.swap(v1, minus)

    with pytest.raises(ValueError):
        eq.swap(v1, mult)

    with pytest.raises(ValueError):
        eq.swap(v1, root)

    # Swap the root
    eq.swap(eq.root, v1)
    assert eq._value is None
    assert v1.value, eq()

    assert noObserversInGlobalBuilders
    return


def testEmbeddedEquation(make_args, noObserversInGlobalBuilders):
    """Test a simple function."""

    # Make some variables
    v1, v2, v3, v4, c = make_args(5)
    c.name = "c"
    c.const = True

    # Make some operations
    mult = literals.MultiplicationOperator()
    mult2 = literals.MultiplicationOperator()
    plus = literals.AdditionOperator()
    minus = literals.SubtractionOperator()

    # Create the equation c*(v1+v3)*(v4-v2)
    plus.addLiteral(v1)
    plus.addLiteral(v3)
    minus.addLiteral(v4)
    minus.addLiteral(v2)
    mult.addLiteral(plus)
    mult.addLiteral(minus)
    mult2.addLiteral(mult)
    mult2.addLiteral(c)

    # Set the values of the variables.
    # The equation should evaluate to 2.5*(1+3)*(4-2) = 20
    v1.setValue(1)
    v2.setValue(2)
    v3.setValue(3)
    v4.setValue(4)
    c.setValue(2.5)

    # Make an equation and test
    root = Equation("root", mult2)
    eq = Equation("eq", root)

    assert eq._value is None
    args = eq.args
    assert v1 in args
    assert v2 in args
    assert v3 in args
    assert v4 in args
    assert c not in args
    assert root is eq.root

    assert v1 is eq.v1
    assert v2 is eq.v2
    assert v3 is eq.v3
    assert v4 is eq.v4

    # Make sure the right messages get sent
    v1.value = 0
    assert root._value is None
    assert eq._value is None
    v1.value = 1

    assert 20 == eq()  # 20 = 2.5*(1+3)*(4-2)
    assert 20 == eq.getValue()  # same as above
    assert 20 == eq.value  # same as above
    assert 25 == eq(v1=2)  # 25 = 2.5*(2+3)*(4-2)
    assert 50 == eq(v2=0)  # 50 = 2.5*(2+3)*(4-0)
    assert 30 == eq(v3=1)  # 30 = 2.5*(2+1)*(4-0)
    assert 0 == eq(v4=0)  # 20 = 2.5*(2+1)*(0-0)

    # Try some swapping.
    eq.swap(v4, v1)
    assert eq._value is None
    assert 15 == eq()  # 15 = 2.5*(2+1)*(2-0)
    args = eq.args
    assert v4 not in args

    # Try to create a dependency loop
    with pytest.raises(ValueError):
        eq.swap(v1, eq.root)

    with pytest.raises(ValueError):
        eq.swap(v1, plus)

    with pytest.raises(ValueError):
        eq.swap(v1, minus)

    with pytest.raises(ValueError):
        eq.swap(v1, mult)

    with pytest.raises(ValueError):
        eq.swap(v1, root)

    # Swap the root
    eq.swap(eq.root, v1)
    assert eq._value is None
    assert v1.value == eq()

    assert noObserversInGlobalBuilders
    return
