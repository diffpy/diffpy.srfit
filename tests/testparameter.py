#!/usr/bin/env python
"""Tests for refinableobj module."""

import unittest

from diffpy.srfit.fitbase.parameter import Parameter, ParameterProxy

class TestParameter(unittest.TestCase):

    def testSetValue(self):
        """Test initialization."""
        l = Parameter("l")

        l.setValue(3.14)
        self.assertAlmostEqual(l.getValue(), 3.14)

        # Try array
        import numpy
        x = numpy.arange(0, 10, 0.1)
        l.setValue(x)
        self.assertTrue( l.getValue() is x )

        # Change the array
        y = numpy.arange(0, 10, 0.5)
        l.setValue(y)
        self.assertTrue( l.getValue() is y )

        # Back to scalar
        l.setValue(1.01)
        self.assertAlmostEqual(l.getValue(), 1.01)
        return

class TestParameterProxy(unittest.TestCase):

    def testProxy(self):
        """Test the ParameterProxy class."""
        l = Parameter("l", 3.14)

        # Try Accessor adaptation
        la = ParameterProxy("l2", l)

        self.assertEqual("l2", la.name)
        self.assertEqual(l.getValue(), la.getValue())

        # Change the parameter
        l.setValue(2.3)
        self.assertEqual(l.getValue(), la.getValue())

        # Change the proxy
        la.setValue(3.2)
        self.assertEqual(l.getValue(), la.getValue())

        return


if __name__ == "__main__":

    unittest.main()

