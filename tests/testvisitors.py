#!/usr/bin/env python
"""Tests for refinableobj module."""

import diffpy.srfit.equation.visitors as visitors
import diffpy.srfit.equation.literals as literals
import unittest

def _makeArgs(num):
    args = []
    for i in xrange(num):
        args.append(literals.Argument())
    return args


class TestEvaluator(unittest.TestCase):

    def testSimpleFunction(self):
        """Test a simple function."""

        # Make some variables
        v1, v2, v3, v4 = _makeArgs(4)

        # Make some operations
        mult = literals.MultiplicationOperator()
        plus = literals.AdditionOperator()
        minus = literals.SubtractionOperator()

        # Create the equation (v1+v3)*(v4-v2)
        plus.addLiteral(v1)
        plus.addLiteral(v3)
        minus.addLiteral(v4)
        minus.addLiteral(v2)
        mult.addLiteral(plus)
        mult.addLiteral(minus)

        # Set the values of the variables.
        # The equation should evaluate to (1+3)*(4-2) = 8
        v1.setValue(1)
        v2.setValue(2)
        v3.setValue(3)
        v4.setValue(4)

        # Evaluate this
        evaluator = visitors.Evaluator()
        mult.identify(evaluator)
        evaluator.clicker.click()
        self.assertEqual(8, evaluator.value)
        self.assertTrue(evaluator.clicker > plus.clicker)
        self.assertTrue(evaluator.clicker > minus.clicker)
        self.assertTrue(evaluator.clicker > mult.clicker)
        self.assertTrue(evaluator.clicker > v1.clicker)
        self.assertTrue(evaluator.clicker > v2.clicker)
        self.assertTrue(evaluator.clicker > v3.clicker)
        self.assertTrue(evaluator.clicker > v4.clicker)

        # Change one of the variables
        v1.setValue(7)
        self.assertTrue(v1.clicker > evaluator.clicker)
        self.assertTrue(plus.clicker > evaluator.clicker)
        self.assertTrue(mult.clicker > evaluator.clicker)
        self.assertTrue(v1.clicker == plus.clicker)
        self.assertTrue(v1.clicker == mult.clicker)
        mult.identify(evaluator)
        evaluator.clicker.click()
        self.assertEqual(20, evaluator.value)
        return
        
    def testCustomFunction(self):
        """Test a custom function."""
        evaluator = visitors.Evaluator()

        # Make some variables
        v1, v2, v3 = _makeArgs(3)

        # Make some operations
        mult = literals.MultiplicationOperator()
        plus = literals.AdditionOperator()
        minus = literals.SubtractionOperator()
        import numpy
        sin = literals.UfuncOperator(numpy.sin)

        # Create the equation v1*sin(v2) + v3
        sin.addLiteral(v2)
        mult.addLiteral(sin)
        mult.addLiteral(v1)
        plus.addLiteral(v3)
        plus.addLiteral(mult)

        # Give the variables values. The equation should evaluate to
        # 2*sin(pi/6)+3
        v1.setValue(2)
        v2.setValue(numpy.pi/6)
        v3.setValue(3)

        # Evaluate this
        plus.identify(evaluator)
        evaluator.clicker.click()
        self.assertAlmostEqual(4, evaluator.value)
        return

    def testArray(self):
        """Test functions operating on a numpy array."""
        evaluator = visitors.Evaluator()

        # Make some variables
        v1, v2 = _makeArgs(2)

        # Make some operations
        mult = literals.MultiplicationOperator()
        import numpy
        sum = literals.Operator()
        sum.name = sum.symbol = "sum"
        sum.nin = 1
        sum.operation = numpy.sum

        # Create the equation sum(v1*v2)
        sum.addLiteral(mult)
        mult.addLiteral(v1)
        mult.addLiteral(v2)

        # Give the variables values. 
        # v1 = 5
        # v2 = array([0, 1, 2, 3, 4])
        v1.setValue(5)
        v2.setValue(numpy.arange(5))

        # Evaluate this. It should give 50.
        sum.identify(evaluator)
        evaluator.clicker.click()
        self.assertAlmostEqual(50, evaluator.value)
        return


class TestValidator(unittest.TestCase):

    def testSimpleFunction(self):
        """Test a simple function."""

        # Make some variables
        v1, v2, v3, v4 = _makeArgs(4)

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

if __name__ == "__main__":

    unittest.main()

