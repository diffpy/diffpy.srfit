#!/usr/bin/env python
"""Tests for refinableobj module."""

import diffpy.srfit.equation.visitors as visitors
import diffpy.srfit.equation.literals as literals
import unittest


class TestEvaluator(unittest.TestCase):

    def _makeArgs(self, num):
        args = []
        for i in xrange(num):
            args.append(literals.Argument())
        return args

    def testSimpleFunction(self):
        """Test a simple function."""
        evaluator = visitors.Evaluator()

        # Make some variables
        v1, v2, v3, v4 = self._makeArgs(4)

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
        v1, v2, v3 = self._makeArgs(3)

        # Make some operations
        mult = literals.MultiplicationOperator()
        plus = literals.AdditionOperator()
        minus = literals.SubtractionOperator()
        import numpy
        sin = literals.Operator()
        sin.name = sin.symbol = "sin"
        sin.numargs = 1
        sin.operation = numpy.sin

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
        v1, v2 = self._makeArgs(2)

        # Make some operations
        mult = literals.MultiplicationOperator()
        import numpy
        sum = literals.Operator()
        sum.name = sum.symbol = "sum"
        sum.numargs = 1
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



if __name__ == "__main__":

    unittest.main()

