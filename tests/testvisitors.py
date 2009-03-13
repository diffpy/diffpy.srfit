#!/usr/bin/env python
"""Tests for refinableobj module."""

import diffpy.srfit.equation.visitors as visitors
import diffpy.srfit.equation.literals as literals
import unittest

from utils import _makeArgs


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
        evaluator.click()
        self.assertEqual(8, evaluator.value)
        self.assertTrue(evaluator._clicker > plus.clicker)
        self.assertTrue(evaluator._clicker > minus.clicker)
        self.assertTrue(evaluator._clicker > mult.clicker)
        self.assertTrue(evaluator._clicker > v1.clicker)
        self.assertTrue(evaluator._clicker > v2.clicker)
        self.assertTrue(evaluator._clicker > v3.clicker)
        self.assertTrue(evaluator._clicker > v4.clicker)

        # Change one of the variables
        v1.setValue(7)
        self.assertTrue(v1.clicker > evaluator._clicker)
        self.assertTrue(plus.clicker > evaluator._clicker)
        self.assertTrue(mult.clicker > evaluator._clicker)
        self.assertTrue(v1.clicker == plus.clicker)
        self.assertTrue(v1.clicker == mult.clicker)
        mult.identify(evaluator)
        evaluator.click()
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
        sin = literals.UFuncOperator(numpy.sin)

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
        evaluator.click()
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
        evaluator.click()
        self.assertAlmostEqual(50, evaluator.value)
        return

    def testSimplePartition(self):
        """Test a function with a partition in it."""
        p1 = literals.Partition("p1")

        # Make the partition 1|2
        v1, v2, v3 = _makeArgs(3)
        p1.addArgument(v1)
        p1.addArgument(v2)

        # make the equation sin(3*partition).
        mult = literals.MultiplicationOperator()
        import numpy
        sin = literals.UFuncOperator(numpy.sin)

        mult.addLiteral(p1)
        mult.addLiteral(v3)
        sin.addLiteral(mult)

        # With no operations able to combine, the equation will be combined
        # after the root operation.
        # sin(3*1) + sin(3*2)
        evaluator = visitors.Evaluator()
        sin.identify(evaluator)
        evaluator.click()
        self.assertEqual(numpy.sin(3*1)+numpy.sin(3*2), evaluator.value)

        # Change the value of an argument.
        v1.setValue(2)
        sin.identify(evaluator)
        evaluator.click()
        self.assertEqual(numpy.sin(3*2)+numpy.sin(3*2), evaluator.value)
        v1.setValue(1)

        # Now let the '*' operator combine. This should give
        # sin(3*1 + 3*2)
        mult.setCombine()
        sin.identify(evaluator)
        evaluator.click()
        self.assertEqual(numpy.sin(3*1 + 3*2), evaluator.value)

        # Do that again to prove that it works.
        sin.identify(evaluator)
        evaluator.click()
        self.assertEqual(numpy.sin(3*1 + 3*2), evaluator.value)
        return

    def testCombinePartition(self):
        """Test a function with a partition in it and the CombineOperator."""

        p1 = literals.Partition("p1")
        # Make the partition 1|2
        v1, v2, v3 = _makeArgs(3)
        p1.addArgument(v1)
        p1.addArgument(v2)

        # make the equation sin(3*partition).
        mult = literals.MultiplicationOperator()
        import numpy
        sin = literals.UFuncOperator(numpy.sin)
        comb = literals.CombineOperator()

        # This combines after the '*' operator. This should give
        # sin(3*1 + 3*2)
        mult.addLiteral(p1)
        mult.addLiteral(v3)
        comb.addLiteral(mult)
        sin.addLiteral(comb)

        evaluator = visitors.Evaluator()
        sin.identify(evaluator)
        evaluator.click()
        self.assertEqual(numpy.sin(3*1+3*2), evaluator.value)

        evaluator = visitors.Evaluator()
        sin.identify(evaluator)
        self.assertEqual(numpy.sin(3*1 + 3*2), evaluator.value)
        return

    def testTwoPartitions(self):
        """Test a function with two partitions."""

        p1 = literals.Partition("p1")
        p2 = literals.Partition("p2")
        v1, v2, v3, v4, v5 = _makeArgs(5)

        # make the equation (4|5) + ( (1|2) + 3 )
        p1.addArgument(v1)
        p1.addArgument(v2)

        p2.addArgument(v4)
        p2.addArgument(v5)

        addroot = literals.AdditionOperator()
        add = literals.AdditionOperator()

        add.addLiteral(p1)
        add.addLiteral(v3)
        addroot.addLiteral(p2)
        addroot.addLiteral(add)

        # The partitions colide in the top '+' and are therefore collected
        # before then. This should give
        # (4|5) + (4|5) = 9 + 9 = 18
        evaluator = visitors.Evaluator()
        addroot.identify(evaluator)
        evaluator.click()
        self.assertEqual(18, evaluator.value)
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

class TestArgFinder(unittest.TestCase):

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

        # now get the args
        argfinder = visitors.ArgFinder()
        mult.identify(argfinder)
        args = argfinder.args
        self.assertEqual(4, len(args))
        self.assertTrue(v1 in args)
        self.assertTrue(v2 in args)
        self.assertTrue(v3 in args)
        self.assertTrue(v4 in args)

        return

    def testSimplePartition(self):
        """Test an equation with a partition."""
        p1 = literals.Partition("p1")

        # Make the partition 1|2
        v1, v2, v3 = _makeArgs(3)
        p1.addArgument(v1)
        p1.addArgument(v2)

        # make the equation sin(3*partition).
        mult = literals.MultiplicationOperator()
        import numpy
        sin = literals.UFuncOperator(numpy.sin)

        mult.addLiteral(p1)
        mult.addLiteral(v3)
        sin.addLiteral(mult)

        argfinder = visitors.ArgFinder()
        sin.identify(argfinder)
        args = argfinder.args
        self.assertEqual(3, len(args))
        self.assertTrue(v1 in args)
        self.assertTrue(v2 in args)
        self.assertTrue(v3 in args)
        return

if __name__ == "__main__":

    unittest.main()

