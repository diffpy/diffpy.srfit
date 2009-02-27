#!/usr/bin/env python
"""Tests for refinableobj module."""

import diffpy.srfit.equation.utils as utils

import unittest

import numpy


class TestEquationParser(unittest.TestCase):

    def testParseEquation(self):

        from numpy import exp, sin, divide, sqrt, array_equal, e

        eq = utils.makeEquation("A*sin(0.5*x)+divide(B,C)")
        A = 1
        x = numpy.pi
        B = 4.0
        C = 2.0
        eq.A.setValue(A)
        eq.x.setValue(x)
        eq.B.setValue(B)
        eq.C.setValue(C)
        f = lambda A, x, B, C: A*sin(0.5*x)+divide(B,C)
        self.assertAlmostEqual(eq(), f(A,x,B,C))

        eq = utils.makeEquation("sqrt(e**(-0.5*(x/sigma)**2))")
        x = numpy.arange(0, 20, 0.05)
        sigma = 0.01
        eq.x.setValue(x)
        eq.sigma.setValue(sigma)
        f = lambda x, sigma : sqrt(e**(-0.5*(x/sigma)**2))
        self.assertTrue(array_equal(eq(), f(x,sigma)))

if __name__ == "__main__":

    unittest.main()

