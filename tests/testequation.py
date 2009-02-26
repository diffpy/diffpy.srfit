#!/usr/bin/env python
"""Tests for refinableobj module."""

import diffpy.srfit.equation.visitors as visitors
import diffpy.srfit.equation.literals as literals
from diffpy.srfit.equation import Equation
import unittest

from utils import _makeArgs

class TestEquation(unittest.TestCase):

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

        # Make the equation
        eq = Equation(mult)
        args = eq.args.values()
        self.assertTrue(v1 in args)
        self.assertTrue(v2 in args)
        self.assertTrue(v3 in args)
        self.assertTrue(v4 in args)
        self.assertTrue(mult is eq.root)

        self.assertEqual(v1, eq.v1)
        self.assertEqual(v2, eq.v2)
        self.assertEqual(v3, eq.v3)
        self.assertEqual(v4, eq.v4)

        self.assertEqual(8, eq())
        self.assertEqual(10, eq(v1=2))
        self.assertEqual(20, eq(v2=0))
        self.assertEqual(12, eq(v3=1))
        self.assertEqual(0, eq(v4=0))
        return

if __name__ == "__main__":

    unittest.main()

