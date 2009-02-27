#!/usr/bin/env python
"""Tests for refinableobj module."""

import diffpy.srfit.equation.utils as utils

import unittest


class TestEquationParser(unittest.TestCase):

    def testParseEquation(self):

        eq = utils.makeEquation("A*sin(0.5*x)+divide(B,C)")
        eq.A.setValue(1)
        eq.x.setValue(3.14159)
        eq.B.setValue(4.0)
        eq.C.setValue(2.0)
        self.assertAlmostEqual(3, eq())

        eq = utils.makeEquation("A*sin(0.5*x)+B%C")
        eq.A.setValue(1)
        eq.x.setValue(3.14159)
        eq.B.setValue(5.0)
        eq.C.setValue(3.0)
        self.assertAlmostEqual(3, eq())

if __name__ == "__main__":

    unittest.main()

