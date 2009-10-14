#!/usr/bin/env python
"""Tests for refinableobj module."""

import unittest

import numpy

import diffpy.srfit.equation.literals as literals


class TestArgument(unittest.TestCase):

    def testInit(self):
        """Test that everthing initializes as expected."""
        a = literals.Argument(name = "a", value = 0)

        self.assertEqual("a", a.name)
        self.assertEqual(0, a._value)
        self.assertEqual(a, a._target)
        self.assertEqual(set(), a._proxies)
        self._proxies = set()
        return

    def testValue(self):
        """Test value setting."""

        a = literals.Argument(name = "a")

        # Test error when there is no value
        self.assertRaises(ValueError, a.getValue)

        # Test setting value
        a.setValue(3.14)
        self.assertAlmostEqual(a.getValue(), 3.14)
        return

    def testTargeting(self):
        """Test adding a target."""

        a = literals.Argument(name = "a")
        b = literals.Argument(name = "b")

        a.setValue(1)
        b.setValue(2)

        # Make sure the target gets set properly set
        a.setTarget(b)
        self.assertTrue(a in b._proxies)
        self.assertTrue(a._target is b)
        self.assertTrue(b._target is b)
        self.assertTrue(a._getDeepTarget() is b)
        self.assertTrue(b._getDeepTarget() is b)

        # Check that the values are properly set
        self.assertEqual(a.getValue(), b.getValue())

        a.setValue(3)
        self.assertEqual(1, a._value)
        self.assertEqual(3, b._value)
        self.assertEqual(a.getValue(), b.getValue())

        b.setValue(4)
        self.assertEqual(1, a._value)
        self.assertEqual(4, b._value)
        self.assertEqual(a.getValue(), b.getValue())


        # Make another target and put it deeper in the chain
        c = literals.Argument(name = "c")
        c.setValue(5)
        b.setTarget(c)
        self.assertTrue(a in b._proxies)
        self.assertTrue(b in c._proxies)
        self.assertTrue(a._target is b)
        self.assertTrue(b._target is c)
        self.assertTrue(c._target is c)
        self.assertTrue(a._getDeepTarget() is c)
        self.assertTrue(b._getDeepTarget() is c)
        self.assertTrue(c._getDeepTarget() is c)

        # check values
        self.assertEqual(1, a._value)
        self.assertEqual(4, b._value)
        self.assertEqual(5, c._value)
        self.assertEqual(5, a.getValue())
        self.assertEqual(5, b.getValue())
        self.assertEqual(5, c.getValue())

        a.setValue(6)
        self.assertEqual(1, a._value)
        self.assertEqual(4, b._value)
        self.assertEqual(6, c._value)
        self.assertEqual(6, a.getValue())
        self.assertEqual(6, b.getValue())
        self.assertEqual(6, c.getValue())

        b.setValue(7)
        self.assertEqual(1, a._value)
        self.assertEqual(4, b._value)
        self.assertEqual(7, c._value)
        self.assertEqual(7, a.getValue())
        self.assertEqual(7, b.getValue())
        self.assertEqual(7, c.getValue())

        c.setValue(8)
        self.assertEqual(1, a._value)
        self.assertEqual(4, b._value)
        self.assertEqual(8, c._value)
        self.assertEqual(8, a.getValue())
        self.assertEqual(8, b.getValue())
        self.assertEqual(8, c.getValue())

        # Detach a from b
        a.setTarget(a)
        self.assertTrue(a._target is a)
        self.assertTrue(b._target is c)
        self.assertTrue(c._target is c)
        self.assertTrue(a._getDeepTarget() is a)
        self.assertTrue(b._getDeepTarget() is c)
        self.assertTrue(c._getDeepTarget() is c)

        # Check the values again
        self.assertRaises(ValueError, a.getValue)
        self.assertEqual(None, a._value)
        self.assertEqual(4, b._value)
        self.assertEqual(8, c._value)
        self.assertEqual(8, b.getValue())
        self.assertEqual(8, c.getValue())

        # Detect loops
        self.assertRaises(ValueError, c.setTarget, b)

        a.setTarget(b)
        self.assertTrue(a in b._proxies)
        self.assertTrue(b in c._proxies)
        self.assertRaises(ValueError, c.setTarget, b)
        self.assertRaises(ValueError, c.setTarget, a)
        self.assertRaises(ValueError, b.setTarget, a)

        a.setTarget(c)
        self.assertTrue(a not in b._proxies)
        self.assertTrue(a in c._proxies)
        self.assertTrue(b in c._proxies)
        self.assertRaises(ValueError, c.setTarget, b)
        self.assertRaises(ValueError, c.setTarget, a)


        return

class TestOperator(unittest.TestCase):

    def testInit(self):
        """Test that everthing initializes as expected."""
        op = literals.Operator(name = "add", symbol = "+", operation =
                numpy.add, nin = 2)

        self.assertEqual("add", op.name)
        self.assertEqual("+", op.symbol)
        self.assertEqual(numpy.add, op.operation)
        self.assertEqual(2, op.nin)
        self.assertEqual(1, op.nout)
        self.assertEqual([], op.args)
        self.assertEqual(None, op._value)
        self.assertEqual(op, op._target)
        self.assertEqual(set(), op._proxies)
        return

    def testAddLiteral(self):
        """Test adding a literal to an operator node."""
        op = literals.Operator(name = "add", symbol = "+", operation =
                numpy.add, nin = 2)

        self.assertRaises(ValueError, op.getValue)
        op._value = 1
        self.assertEqual(op.getValue(), 1)

        # Test addition and operations
        a = literals.Argument(name = "a", value = 0)
        b = literals.Argument(name = "b", value = 0)

        op.addLiteral(a)
        self.assertRaises(ValueError, op.getValue)

        op.addLiteral(b)
        self.assertAlmostEqual(op.getValue(), 0)

        a.setValue(1)
        b.setValue(2)
        self.assertTrue(op._value is None)
        self.assertAlmostEqual(op.getValue(), 3)

        a.setValue(None)
        self.assertRaises(ValueError, op.getValue)

        # Test for self-references

        # Try to add self
        op = literals.Operator(name = "add", symbol = "+", operation =
                numpy.add, nin = 2)
        op.addLiteral(a)
        self.assertRaises(ValueError, op.addLiteral, op)

        # Try to add proxy to self
        b.setTarget(op)
        self.assertRaises(ValueError, op.addLiteral, b)

        # Try to add more distant proxy
        c = literals.Argument(name = "c")
        c.setTarget(b)
        self.assertRaises(ValueError, op.addLiteral, c)

        # Try to add argument that contains self
        op2 = literals.Operator(name = "sub", symbol = "-", operation =
                numpy.subtract, nin = 2)
        op2.addLiteral(op)
        self.assertRaises(ValueError, op.addLiteral, op2)

        # Try to add argument that proxies self
        c.setTarget(op2)
        self.assertRaises(ValueError, op.addLiteral, c)

        # Try to add argument that contains proxy to self
        op3 = literals.Operator(name = "mul", symbol = "*", operation =
                numpy.multiply, nin = 2)
        op3.addLiteral(c)
        self.assertRaises(ValueError, op2.addLiteral, op3)

        # Continue making the equation
        a = literals.Argument(name = "a", value = 1)
        b = literals.Argument(name = "b", value = 2)
        c = literals.Argument(name = "c", value = 3)
        op = literals.Operator(name = "add", symbol = "+", operation =
                numpy.add, nin = 2)
        op2 = literals.Operator(name = "sub", symbol = "-", operation =
                numpy.subtract, nin = 2)

        op.addLiteral(op2)
        op.addLiteral(a)
        op2.addLiteral(b)
        op2.addLiteral(c)
        self.assertAlmostEqual(-1, op2.getValue())
        self.assertAlmostEqual(0, op.getValue())

        c.setValue(0)
        self.assertTrue(op2._value is None)
        self.assertTrue(op._value is None)
        self.assertAlmostEqual(2, op2.getValue())
        self.assertAlmostEqual(3, op.getValue())

        return



if __name__ == "__main__":

    unittest.main()

