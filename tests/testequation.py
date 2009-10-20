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
        v1.value = 1
        v2.value = 2
        v3.value = 3
        v4.value = 4
        c.value = 2.5

        # Make the equation
        eq = Equation(mult2)
        args = eq.args
        # The arguments are found depth-first
        self.assertEquals(args, [v1, v3, v4, v2])
        self.assertTrue(v1 in args)
        self.assertTrue(v2 in args)
        self.assertTrue(v3 in args)
        self.assertTrue(v4 in args)
        self.assertTrue(c not in args)
        self.assertTrue(mult2 is eq.root)

        self.assertAlmostEqual(20, eq()) # 20 = 2.5*(1+3)*(4-2)
        self.assertAlmostEqual(25, eq(2)) # 25 = 2.5*(2+3)*(4-2)
        self.assertAlmostEqual(10, eq(2, 0)) # 10 = 2.5*(2+0)*(4-2)
        self.assertAlmostEqual(-5, eq(2, 0, 1)) # -5 = 2.5*(2+0)*(1-2)
        self.assertAlmostEqual(-5, eq(2, 0, 0, 1)) # -5 = 2.5*(2+0)*(0-1)

        return

if __name__ == "__main__":

    unittest.main()

