#!/usr/bin/env python
"""Tests for refinableobj module."""

import unittest

from diffpy.srfit.fitbase.parameter import Parameter


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

if __name__ == "__main__":

    unittest.main()

