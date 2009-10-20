#!/usr/bin/env python
"""Tests for refinableobj module."""

import diffpy.srfit.equation.visitors as visitors
import diffpy.srfit.equation.literals as literals
import unittest

from utils import _makeArgs, _makeNodes

class TestValidator(unittest.TestCase):

    def testSimpleFunction(self):
        """Test a simple function."""

        # Make some variables
        v1, v2, v3, v4 = _makeNodes(4)

        # Make some operations
        mult = literals.MultiplicationOperator()
        plus = literals.AdditionOperator()
        minus = literals.SubtractionOperator()

        # Let's hobble the plus operator
        plus.name = None
        plus.symbol = None
        plus.operation = None

        # Partially define the equation (v1+v3)*(v4-v2). Let's only give one
        # variable to the '-' operation.
        plus.addLiteral(v1)
        plus.addLiteral(v3)
        minus.addLiteral(v4)
        mult.addLiteral(plus)
        mult.addLiteral(minus)

        # Now validate
        validator = visitors.Validator()
        mult.identify(validator)
        print validator.errors
        self.assertEqual(4, len(validator.errors))

        # Fix the equation
        minus.addLiteral(v3)
        validator.reset()
        mult.identify(validator)
        self.assertEqual(3, len(validator.errors))

        # Fix the name of plus
        plus.name = "add"
        validator.reset()
        mult.identify(validator)
        self.assertEqual(2, len(validator.errors))

        # Fix the symbol of plus
        plus.symbol = "+"
        validator.reset()
        mult.identify(validator)
        self.assertEqual(1, len(validator.errors))

        # Fix the operation of plus
        import numpy
        plus.operation = numpy.add
        validator.reset()
        mult.identify(validator)
        self.assertEqual(0, len(validator.errors))

        # Add another literal to minus
        minus.addLiteral(v1)
        validator.reset()
        mult.identify(validator)
        self.assertEqual(1, len(validator.errors))

        return

class TestArgFinder(unittest.TestCase):

    def testSimpleFunction(self):
        """Test a simple function."""

        # Make some variables
        v1, v2, v3, v4 = _makeNodes(4)

        # Make some operations
        ns = v1.namespace
        ns["mult"] = mult = literals.MultiplicationOperator()
        ns["plus"] = plus = literals.AdditionOperator()
        ns["minus"] = minus = literals.SubtractionOperator()
        nmult = literals.Node(mult, ns)
        nplus = literals.Node(plus, ns)
        nminus = literals.Node(minus, ns)

        # Create the equation (v1+v3)*(v4-v2)
        plus.addLiteral(v1.target)
        plus.addLiteral(v3.target)
        minus.addLiteral(v4.target)
        minus.addLiteral(v2.target)
        mult.addLiteral(plus)
        mult.addLiteral(minus)

        # Set the values of the variables.
        # The equation should evaluate to (1+3)*(4-2) = 8
        v1.setValue(1)
        v2.setValue(2)
        v3.setValue(3)
        v4.setValue(4)

        # now get the args
        argfinder = visitors.ArgFinder()
        mult.identify(argfinder)
        args = argfinder.args
        self.assertEqual(4, len(args))
        self.assertTrue(v1.target in args)
        self.assertTrue(v2.target in args)
        self.assertTrue(v3.target in args)
        self.assertTrue(v4.target in args)
        return

    def testArg(self):
        """Test just an Argument equation."""
        # Make some variables
        v1 = _makeArgs(1)[0]

        argfinder = visitors.ArgFinder()
        v1.identify(argfinder)
        args = argfinder.args

        self.assertEquals(1, len(args))
        self.assertTrue(args[0] is v1)
        return

if __name__ == "__main__":

    unittest.main()

