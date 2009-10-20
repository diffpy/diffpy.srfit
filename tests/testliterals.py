#!/usr/bin/env python
"""Tests for refinableobj module."""

import unittest

import numpy

import diffpy.srfit.equation.literals as literals


class TestArgument(unittest.TestCase):

    def testInit(self):
        """Test that everthing initializes as expected."""
        a = literals.Argument(value = 0)
        self.assertEqual(0, a._value)
        return

    def testValue(self):
        """Test value setting."""

        a = literals.Argument()

        # Test error when there is no value
        self.assertRaises(ValueError, a.evaluate, None)

        # Test setting value
        a.setValue(3.14)
        self.assertAlmostEqual(a._value, 3.14)
        return

class TestOperator(unittest.TestCase):

    def testInit(self):
        """Test that everthing initializes as expected."""
        op = literals.Operator(symbol = "+", operation = numpy.add, nin = 2)

        self.assertEqual("+", op.symbol)
        self.assertEqual(numpy.add, op.operation)
        self.assertEqual(2, op.nin)
        self.assertEqual(1, op.nout)
        self.assertEqual(None, op._value)
        return

    def testValue(self):
        """Test value."""
        # Test addition and operations
        ns = {}
        op = literals.Operator(symbol = "+", operation = numpy.add, nin = 2)
        ns["a"] = literals.Argument(value = 0)
        ns["b"] = literals.Argument(value = 0)
        ns["op"] = op

        a = literals.Node("a", ns)
        b = literals.Node("b", ns)
        node = literals.Node("op", ns)
        node.addNode(a)
        node.addNode(b)

        op.evaluate(node)
        self.assertAlmostEquals(0, op._value)

        # Test update from the nodes
        a.setValue(4)
        self.assertTrue(op._value is None)
        self.assertAlmostEqual(4, op.evaluate(node))

        b.value = 2
        self.assertTrue(op._value is None)
        self.assertAlmostEqual(6, op.evaluate(node))

        return


class TestNode(unittest.TestCase):

    def testInit(self):
        """Test that everthing initializes as expected."""

        ns = {}
        ns["a"] = literals.Argument(0)

        node = literals.Node("a", ns)
        self.assertEqual("a", node.name)
        self.assertEqual(ns, node.namespace)
        self.assertAlmostEqual(0, node.value)
        self.assertTrue(ns["a"] is node.target)
        self.assertEqual(node.args, [])
        return

    def testAddNode(self):
        """Test adding a node."""
        ns = {}
        ns["a"] = literals.Argument(0)
        ns["b"] = literals.Argument(0)
        anode = literals.Node("a", ns)
        bnode = literals.Node("b", ns)
        ns["add"] = literals.Operator(symbol = "+", operation = numpy.add, nin = 2)
        addnode = literals.Node("add", ns)
        addnode.addNode(anode)
        addnode.addNode(bnode)
        self.assertEquals([anode, bnode], addnode.args)
        return

    def testValue(self):
        """Test setting and getting the value."""
        ns = {}
        ns["a"] = literals.Argument(0)
        ns["b"] = literals.Argument(0)
        anode = literals.Node("a", ns)
        bnode = literals.Node("b", ns)
        self.assertAlmostEqual(0, anode.value)
        self.assertAlmostEqual(0, bnode.value)
        anode.setValue(1)
        self.assertAlmostEqual(1, anode.value)
        self.assertAlmostEqual(1, ns["a"]._value)

        # Add two arguments
        ns["add"] = literals.Operator(symbol = "+", operation = numpy.add, nin = 2)
        addnode = literals.Node("add", ns)
        addnode.addNode(anode)
        addnode.addNode(bnode)

        self.assertAlmostEqual(1, addnode.value)

        bnode.value = 2
        self.assertAlmostEqual(3, addnode.value)
        
        # Add an operator and another argument.
        ns["c"] = literals.Argument(0)
        cnode = literals.Node("c", ns)
        ns["sub"] = literals.Operator(symbol = "-", operation = numpy.subtract,
                nin = 2)
        subnode = literals.Node("sub", ns)
        subnode.addNode(addnode)
        subnode.addNode(cnode)

        anode.value = 2
        self.assertTrue(ns["add"]._value is None)
        self.assertTrue(ns["sub"]._value is None)
        cnode.value = 0
        self.assertAlmostEqual(4, subnode.value)

        # Make sure that not changing the value does not modify the operators
        anode.value = 2
        self.assertFalse(ns["add"]._value is None)
        return

    def testModifyNamespace(self):
        """Test modifying the namespace."""
        ns = {}
        ns["a"] = a = literals.Argument(1)
        b = literals.Argument(2)
        anode = literals.Node("a", ns)
        bnode = literals.Node("b", ns)

        self.assertAlmostEqual(1, anode.value)
        ns["a"] = b
        self.assertAlmostEqual(2, anode.value)

        ns["a"] = a
        ns["b"] = b

        ns["add"] = literals.Operator(symbol = "+", operation = numpy.add, nin = 2)
        addnode = literals.Node("add", ns)
        addnode.addNode(anode)
        addnode.addNode(bnode)
        self.assertAlmostEqual(3, addnode.value)

        ns["add"] = literals.Operator(symbol = "-", operation = numpy.subtract, nin = 2)
        self.assertAlmostEqual(-1, addnode.value)

        return


if __name__ == "__main__":

    unittest.main()

