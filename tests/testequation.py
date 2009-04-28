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
        v1, v2, v3, v4, c = _makeArgs(5)
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

        # Make the equation
        eq = Equation(mult2)
        args = eq.args
        self.assertTrue(v1 in args)
        self.assertTrue(v2 in args)
        self.assertTrue(v3 in args)
        self.assertTrue(v4 in args)
        self.assertTrue(c not in args)
        self.assertTrue(mult2 is eq.root)

        self.assertEqual(v1, eq.v1)
        self.assertEqual(v2, eq.v2)
        self.assertEqual(v3, eq.v3)
        self.assertEqual(v4, eq.v4)

        self.assertEqual(20, eq()) # 20 = 2.5*(1+3)*(4-2)
        self.assertEqual(25, eq(v1=2)) # 25 = 2.5*(2+3)*(4-2)
        self.assertEqual(50, eq(v2=0)) # 50 = 2.5*(2+3)*(4-0)
        self.assertEqual(30, eq(v3=1)) # 30 = 2.5*(2+1)*(4-0)
        self.assertEqual(0, eq(v4=0)) # 20 = 2.5*(2+1)*(0-0)
        return

if __name__ == "__main__":

    unittest.main()

