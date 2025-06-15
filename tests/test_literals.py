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

import unittest

import numpy

import diffpy.srfit.equation.literals as literals
import diffpy.srfit.equation.literals.abcs as abcs

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

        self.assertEqual(None, a.getValue())

        # Test setting value
        a.setValue(3.14)
        self.assertAlmostEqual(3.14, a._value)

        a.setValue(3.14)
        self.assertAlmostEqual(3.14, a.value)
        self.assertAlmostEqual(3.14, a.getValue())
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
        a.setValue(4)
        self.assertTrue(op._value is None)
        self.assertAlmostEqual(4, op.value)
        self.assertAlmostEqual(4, op.getValue())

        b.value = 2
        self.assertTrue(op._value is None)
        self.assertAlmostEqual(6, op.value)

        return

    def testAddLiteral(self):
        """Test adding a literal to an operator node."""
        op = self.op

        self.assertRaises(TypeError, op.getValue)
        op._value = 1
        self.assertEqual(op.getValue(), 1)

        # Test addition and operations
        a = literals.Argument(name="a", value=0)
        b = literals.Argument(name="b", value=0)

        op.addLiteral(a)
        self.assertRaises(TypeError, op.getValue)

        op.addLiteral(b)
        self.assertAlmostEqual(0, op.value)

        a.setValue(1)
        b.setValue(2)
        self.assertAlmostEqual(3, op.value)

        a.setValue(None)
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

if __name__ == "__main__":
    unittest.main()
