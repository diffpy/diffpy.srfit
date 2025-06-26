#!/usr/bin/env python
##############################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2010 The Trustees of Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
##############################################################################
"""Tests for refinableobj module."""

import unittest

from diffpy.srfit.fitbase.parameter import (
    Parameter,
    ParameterAdapter,
    ParameterProxy,
)


class TestParameter(unittest.TestCase):

    def testSetValue(self):
        """Test initialization."""
        par_l = Parameter("l")

        par_l.setValue(3.14)
        self.assertAlmostEqual(3.14, par_l.getValue())

        # Try array
        import numpy

        x = numpy.arange(0, 10, 0.1)
        par_l.setValue(x)
        self.assertTrue(par_l.getValue() is x)
        self.assertTrue(par_l.value is x)

        # Change the array
        y = numpy.arange(0, 10, 0.5)
        par_l.value = y
        self.assertTrue(par_l.getValue() is y)
        self.assertTrue(par_l.value is y)

        # Back to scalar
        par_l.setValue(1.01)
        self.assertAlmostEqual(1.01, par_l.getValue())
        self.assertAlmostEqual(1.01, par_l.value)
        return


class TestParameterProxy(unittest.TestCase):

    def testProxy(self):
        """Test the ParameterProxy class."""
        par_l = Parameter("l", 3.14)

        # Try Accessor adaptation
        la = ParameterProxy("l2", par_l)

        self.assertEqual("l2", la.name)
        self.assertEqual(par_l.getValue(), la.getValue())

        # Change the parameter
        par_l.value = 2.3
        self.assertEqual(par_l.getValue(), la.getValue())
        self.assertEqual(par_l.value, la.value)

        # Change the proxy
        la.value = 3.2
        self.assertEqual(par_l.getValue(), la.getValue())
        self.assertEqual(par_l.value, la.value)

        return


class TestParameterAdapter(unittest.TestCase):

    def testWrapper(self):
        """Test the adapter.

        This adapts a Parameter to the Parameter interface. :)
        """
        par_l = Parameter("l", 3.14)

        # Try Accessor adaptation
        la = ParameterAdapter(
            "l", par_l, getter=Parameter.getValue, setter=Parameter.setValue
        )

        self.assertEqual(par_l.name, la.name)
        self.assertEqual(par_l.getValue(), la.getValue())

        # Change the parameter
        par_l.setValue(2.3)
        self.assertEqual(par_l.getValue(), la.getValue())

        # Change the adapter
        la.setValue(3.2)
        self.assertEqual(par_l.getValue(), la.getValue())

        # Try Attribute adaptation
        la = ParameterAdapter("l", par_l, attr="value")

        self.assertEqual(par_l.name, la.name)
        self.assertEqual("value", la.attr)
        self.assertEqual(par_l.getValue(), la.getValue())

        # Change the parameter
        par_l.setValue(2.3)
        self.assertEqual(par_l.getValue(), la.getValue())

        # Change the adapter
        la.setValue(3.2)
        self.assertEqual(par_l.getValue(), la.getValue())

        return


if __name__ == "__main__":
    unittest.main()
