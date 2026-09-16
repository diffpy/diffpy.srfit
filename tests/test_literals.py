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

import re
import unittest

import numpy
import pytest

import diffpy.srfit.equation.literals as literals
import diffpy.srfit.equation.literals.abcs as abcs
from diffpy.srfit.equation.equationmod import Equation

# ----------------------------------------------------------------------------


class TestArgument(unittest.TestCase):

    def testInit(self):
        """Test that everything initializes as expected."""
        a = literals.Argument()
        self.assertEqual(None, a._value)
        self.assertTrue(False is a.const)
        self.assertTrue(None is a.name)
        return

    def testIdentity(self):
        """Make sure an Argument is an Argument."""
        a = literals.Argument()
        self.assertTrue(issubclass(literals.Argument, abcs.ArgumentABC))
        self.assertTrue(isinstance(a, abcs.ArgumentABC))
        return

    def testValue(self):
        """Test value setting."""
        a = literals.Argument()

        self.assertEqual(None, a.get_value())

        # Test setting value
        a.set_value(3.14)
        self.assertAlmostEqual(3.14, a._value)

        a.set_value(3.14)
        self.assertAlmostEqual(3.14, a.value)
        self.assertAlmostEqual(3.14, a.get_value())
        return


# ----------------------------------------------------------------------------


class TestCustomOperator(unittest.TestCase):

    def setUp(self):
        self.op = literals.makeOperator(
            name="add", symbol="+", operation=numpy.add, nin=2, nout=1
        )
        return

    def testInit(self):
        """Test that everything initializes as expected."""
        op = self.op
        self.assertEqual("+", op.symbol)
        self.assertEqual(numpy.add, op.operation)
        self.assertEqual(2, op.nin)
        self.assertEqual(1, op.nout)
        self.assertEqual(None, op._value)
        self.assertEqual([], op.args)
        return

    def testIdentity(self):
        """Make sure an Argument is an Argument."""
        op = self.op
        self.assertTrue(issubclass(literals.Operator, abcs.OperatorABC))
        self.assertTrue(isinstance(op, abcs.OperatorABC))
        return

    def testValue(self):
        """Test value."""
        # Test addition and operations
        op = self.op
        a = literals.Argument(value=0)
        b = literals.Argument(value=0)

        op.addLiteral(a)
        op.addLiteral(b)

        self.assertAlmostEqual(0, op.value)

        # Test update from the nodes
        a.set_value(4)
        self.assertTrue(op._value is None)
        self.assertAlmostEqual(4, op.value)
        self.assertAlmostEqual(4, op.get_value())

        b.value = 2
        self.assertTrue(op._value is None)
        self.assertAlmostEqual(6, op.value)

        return

    def testAddLiteral(self):
        """Test adding a literal to an operator node."""
        op = self.op

        self.assertRaises(TypeError, op.get_value)
        op._value = 1
        self.assertEqual(op.get_value(), 1)

        # Test addition and operations
        a = literals.Argument(name="a", value=0)
        b = literals.Argument(name="b", value=0)

        op.addLiteral(a)
        self.assertRaises(TypeError, op.get_value)

        op.addLiteral(b)
        self.assertAlmostEqual(0, op.value)

        a.set_value(1)
        b.set_value(2)
        self.assertAlmostEqual(3, op.value)

        a.set_value(None)
        # Test for self-references

        # Try to add self
        op1 = literals.makeOperator(
            name="add", symbol="+", operation=numpy.add, nin=2, nout=1
        )
        op1.addLiteral(a)
        self.assertRaises(ValueError, op1.addLiteral, op1)

        # Try to add argument that contains self
        op2 = literals.makeOperator(
            name="sub", symbol="-", operation=numpy.subtract, nin=2, nout=1
        )
        op2.addLiteral(op1)
        self.assertRaises(ValueError, op1.addLiteral, op2)

        return


# ----------------------------------------------------------------------------


class TestConvolutionOperator(unittest.TestCase):

    def testValue(self):
        """Make sure the convolution operator is working properly."""
        exp = numpy.exp

        x = numpy.linspace(0, 10, 1000)

        mu1 = 4.5
        sig1 = 0.1
        mu2 = 2.5
        sig2 = 0.4

        g1 = exp(-0.5 * ((x - mu1) / sig1) ** 2)
        a1 = literals.Argument(name="g1", value=g1)
        g2 = exp(-0.5 * ((x - mu2) / sig2) ** 2)
        a2 = literals.Argument(name="g2", value=g2)

        op = literals.ConvolutionOperator()
        op.addLiteral(a1)
        op.addLiteral(a2)

        g3c = op.value

        mu3 = mu1
        sig3 = (sig1**2 + sig2**2) ** 0.5
        g3 = exp(-0.5 * ((x - mu3) / sig3) ** 2)
        g3 *= sum(g1) / sum(g3)

        self.assertAlmostEqual(sum(g3c), sum(g3))
        self.assertAlmostEqual(0, sum((g3 - g3c) ** 2))
        return


# ----------------------------------------------------------------------------


class TestArrayOperator(unittest.TestCase):

    def test_value(self):
        """Check ArrayOperator.value."""
        x = literals.Argument("x", 1.0)
        y = literals.Argument("y", 2.0)
        z = literals.Argument("z", 3.0)
        # check empty array
        op = literals.ArrayOperator()
        self.assertEqual(0, len(op.value))
        self.assertTrue(isinstance(op.value, numpy.ndarray))
        # check behavior with 2 arguments
        op.addLiteral(x)
        self.assertTrue(numpy.array_equal([1], op.value))
        op.addLiteral(y)
        op.addLiteral(z)
        self.assertTrue(numpy.array_equal([1, 2, 3], op.value))
        z.value = 7
        self.assertTrue(numpy.array_equal([1, 2, 7], op.value))
        return


# ----------------------------------------------------------------------------
# Literal.getValue is deprecated in favor of Literal.get_value. Every Literal
# in the hierarchy must keep accepting the old name, warn with a message that
# names the replacement, and dispatch to the subclass implementation of
# get_value rather than to Literal's own NotImplementedError stub.


# C1: Argument holds the value directly.
# Expected: getValue warns and returns Argument.get_value.
def test_argument_get_value_deprecated():
    expected_msg = (
        "'diffpy.srfit.equation.literals.Literal.getValue' is deprecated "
        "and will be removed in version 4.0.0. Please use "
        "'diffpy.srfit.equation.literals.Literal.get_value' instead."
    )
    expected_value = 3.5
    literal = literals.Argument(name="a", value=expected_value)

    with pytest.warns(DeprecationWarning, match=re.escape(expected_msg)):
        actual_value = literal.getValue()

    assert actual_value == expected_value


# C2: Operator computes the value from its own literals.
# Expected: getValue warns and returns Operator.get_value.
def test_operator_get_value_deprecated():
    expected_msg = (
        "'diffpy.srfit.equation.literals.Literal.getValue' is deprecated "
        "and will be removed in version 4.0.0. Please use "
        "'diffpy.srfit.equation.literals.Literal.get_value' instead."
    )
    expected_value = 3.5
    operator = literals.AdditionOperator()
    operator.addLiteral(literals.Argument(name="a", value=1.5))
    operator.addLiteral(literals.Argument(name="b", value=2.0))

    with pytest.warns(DeprecationWarning, match=re.escape(expected_msg)):
        actual_value = operator.getValue()

    assert actual_value == expected_value


# C3: Equation evaluates the operator tree at its root.
# Expected: getValue warns and returns Equation.get_value.
def test_equation_get_value_deprecated():
    expected_msg = (
        "'diffpy.srfit.equation.literals.Literal.getValue' is deprecated "
        "and will be removed in version 4.0.0. Please use "
        "'diffpy.srfit.equation.literals.Literal.get_value' instead."
    )
    expected_value = 3.5
    operator = literals.AdditionOperator()
    operator.addLiteral(literals.Argument(name="a", value=1.5))
    operator.addLiteral(literals.Argument(name="b", value=2.0))
    equation = Equation(name="eq", root=operator)

    with pytest.warns(DeprecationWarning, match=re.escape(expected_msg)):
        actual_value = equation.getValue()

    assert actual_value == expected_value


# ----------------------------------------------------------------------------

if __name__ == "__main__":
    unittest.main()
