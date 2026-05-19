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

import numpy as np
import pytest

from diffpy.srfit.fitbase.parameter import (
    Parameter,
    ParameterAdapter,
    ParameterProxy,
)


class TestParameter(unittest.TestCase):

    def testSetValue(self):
        """Test initialization."""
        par_l = Parameter("l")

        par_l.set_value(3.14)
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
        par_l.set_value(1.01)
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
            "l", par_l, getter=Parameter.getValue, setter=Parameter.set_value
        )

        self.assertEqual(par_l.name, la.name)
        self.assertEqual(par_l.getValue(), la.getValue())

        # Change the parameter
        par_l.set_value(2.3)
        self.assertEqual(par_l.getValue(), la.getValue())

        # Change the adapter
        la.set_value(3.2)
        self.assertEqual(par_l.getValue(), la.getValue())

        # Try Attribute adaptation
        la = ParameterAdapter("l", par_l, attr="value")

        self.assertEqual(par_l.name, la.name)
        self.assertEqual("value", la.attr)
        self.assertEqual(par_l.getValue(), la.getValue())

        # Change the parameter
        par_l.set_value(2.3)
        self.assertEqual(par_l.getValue(), la.getValue())

        # Change the adapter
        la.set_value(3.2)
        self.assertEqual(par_l.getValue(), la.getValue())

        return


@pytest.mark.parametrize(
    "lower, upper, expected",
    [
        # User sets both lower and upper bounds explicitly.
        (1, 10, [1, 10]),
        # User sets only a lower bound.
        (2, None, [2, np.inf]),
        # User sets only an upper bound.
        (None, 8, [-np.inf, 8]),
        # User overwrites existing bounds.
        (2, 6, [2, 6]),
    ],
)
def test_bound_range(lower, upper, expected):
    p = Parameter("a", value=5)
    # If testing overwrite, pre-set bounds to see overwrite effect
    if expected == [2, 6]:
        p.bound_range(0, 10)
    p.bound_range(lower_bound=lower, upper_bound=upper)
    actual = p.bounds
    assert actual == expected


@pytest.mark.parametrize(
    "lower, upper, expected",
    [
        # User sets both lower and upper bounds explicitly.
        (1, 10, [1, 10]),
        # User sets only a lower bound.
        (2, None, [2, np.inf]),
        # User sets only an upper bound.
        (None, 8, [-np.inf, 8]),
        # User overwrites existing bounds.
        (2, 6, [2, 6]),
    ],
)
def test_boundRange(lower, upper, expected):
    p = Parameter("a", value=5)
    # If testing overwrite, pre-set bounds to see overwrite effect
    if expected == [2, 6]:
        p.boundRange(0, 10)
    p.boundRange(lower_bound=lower, upper_bound=upper)
    actual = p.bounds
    assert actual == expected


@pytest.mark.parametrize(
    "value, lower_radius, upper_radius, expected",
    [
        # Symmetric radius (upper_radius None, uses lower_radius)
        (10, 2, None, [8, 12]),
        # Asymmetric radius
        (10, 3, 5, [7, 15]),
        # Zero radius
        (4, 0, None, [4, 4]),
        # Current value updated before bounding
        (20, 2, None, [18, 22]),
    ],
)
def test_bound_window(value, lower_radius, upper_radius, expected):
    p = Parameter("a", value=5)
    if value != 5:
        p.set_value(value)
    p.bound_window(lower_radius=lower_radius, upper_radius=upper_radius)
    actual = p.bounds
    assert actual == expected


@pytest.mark.parametrize(
    "value, lower_radius, upper_radius, expected",
    [
        # Symmetric radius (upper_radius None, uses lower_radius)
        (10, 2, None, [8, 12]),
        # Asymmetric radius
        (10, 3, 5, [7, 15]),
        # Zero radius
        (4, 0, None, [4, 4]),
        # Current value updated before bounding
        (20, 2, None, [18, 22]),
    ],
)
def test_boundWindow(value, lower_radius, upper_radius, expected):
    p = Parameter("a", value=5)
    if value != 5:
        p.set_value(value)
    p.boundWindow(lr=lower_radius, ur=upper_radius)
    actual = p.bounds
    assert actual == expected


if __name__ == "__main__":
    unittest.main()
