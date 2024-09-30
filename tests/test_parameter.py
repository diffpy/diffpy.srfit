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

from diffpy.srfit.fitbase.parameter import Parameter, ParameterAdapter, ParameterProxy


class TestParameter(unittest.TestCase):
    def testSetValue(self):
        """Test initialization."""
        l_parameter = Parameter("l")

        l_parameter.setValue(3.14)
        self.assertAlmostEqual(3.14, l_parameter.getValue())

        # Try array
        import numpy

        x = numpy.arange(0, 10, 0.1)
        l_parameter.setValue(x)
        self.assertTrue(l_parameter.getValue() is x)
        self.assertTrue(l_parameter.value is x)

        # Change the array
        y = numpy.arange(0, 10, 0.5)
        l_parameter.value = y
        self.assertTrue(l_parameter.getValue() is y)
        self.assertTrue(l_parameter.value is y)

        # Back to scalar
        l_parameter.setValue(1.01)
        self.assertAlmostEqual(1.01, l_parameter.getValue())
        self.assertAlmostEqual(1.01, l_parameter.value)
        return


class TestParameterProxy(unittest.TestCase):
    def testProxy(self):
        """Test the ParameterProxy class."""
        l_parameter = Parameter("l", 3.14)

        # Try Accessor adaptation
        la = ParameterProxy("l2", l_parameter)

        self.assertEqual("l2", la.name)
        self.assertEqual(l_parameter.getValue(), la.getValue())

        # Change the parameter
        l_parameter.value = 2.3
        self.assertEqual(l_parameter.getValue(), la.getValue())
        self.assertEqual(l_parameter.value, la.value)

        # Change the proxy
        la.value = 3.2
        self.assertEqual(l_parameter.getValue(), la.getValue())
        self.assertEqual(l_parameter.value, la.value)

        return


class TestParameterAdapter(unittest.TestCase):
    def testWrapper(self):
        """Test the adapter.

        This adapts a Parameter to the Parameter interface. :)
        """
        l_parameter = Parameter("l", 3.14)

        # Try Accessor adaptation
        la = ParameterAdapter("l", l_parameter, getter=Parameter.getValue, setter=Parameter.setValue)

        self.assertEqual(l_parameter.name, la.name)
        self.assertEqual(l_parameter.getValue(), la.getValue())

        # Change the parameter
        l_parameter.setValue(2.3)
        self.assertEqual(l_parameter.getValue(), la.getValue())

        # Change the adapter
        la.setValue(3.2)
        self.assertEqual(l_parameter.getValue(), la.getValue())

        # Try Attribute adaptation
        la = ParameterAdapter("l", l_parameter, attr="value")

        self.assertEqual(l_parameter.name, la.name)
        self.assertEqual("value", la.attr)
        self.assertEqual(l_parameter.getValue(), la.getValue())

        # Change the parameter
        l_parameter.setValue(2.3)
        self.assertEqual(l_parameter.getValue(), la.getValue())

        # Change the adapter
        la.setValue(3.2)
        self.assertEqual(l_parameter.getValue(), la.getValue())

        return


if __name__ == "__main__":
    unittest.main()
