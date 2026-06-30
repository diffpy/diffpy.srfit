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

import io
import re
import unittest

import pytest
from numpy import allclose, arange, array, array_equal, ones_like

from diffpy.srfit.exceptions import SrFitError
from diffpy.srfit.fitbase import ProfileParser
from diffpy.srfit.fitbase.profile import Profile


class TestProfile(unittest.TestCase):

    def setUp(self):
        self.profile = Profile()
        return

    def testInit(self):
        profile = self.profile
        self.assertTrue(profile.xobs is None)
        self.assertTrue(profile.yobs is None)
        self.assertTrue(profile.dyobs is None)
        self.assertTrue(profile.x is None)
        self.assertTrue(profile.y is None)
        self.assertTrue(profile.dy is None)
        self.assertTrue(profile.ycalc is None)
        self.assertEqual(profile.meta, {})
        return

    def test_set_observed_profile(self):
        """Test the set_observed_profile method."""
        # Make a profile with defined dy

        x = arange(0, 10, 0.1)
        y = x
        dy = x

        prof = self.profile
        prof.set_observed_profile(x, y, dy)

        self.assertTrue(array_equal(x, prof.xobs))
        self.assertTrue(array_equal(y, prof.yobs))
        self.assertTrue(array_equal(dy, prof.dyobs))

        # Make a profile with undefined dy
        x = arange(0, 10, 0.1)
        y = x
        dy = None

        self.profile.set_observed_profile(x, y, dy)

        self.assertTrue(array_equal(x, prof.xobs))
        self.assertTrue(array_equal(y, prof.yobs))
        self.assertTrue(array_equal(ones_like(prof.xobs), prof.dyobs))

        # Get the ranged profile to make sure its the same
        self.assertTrue(array_equal(x, prof.x))
        self.assertTrue(array_equal(y, prof.y))
        self.assertTrue(array_equal(ones_like(prof.xobs), prof.dy))

        return

    def testSetObservedProfile(self):
        """Test the deprecated setObservedProfile method.

        Remove this test when setObservedProfile is removed in 4.0.0.
        """
        # Make a profile with defined dy

        x = arange(0, 10, 0.1)
        y = x
        dy = x

        prof = self.profile
        prof.set_observed_profile(x, y, dy)

        self.assertTrue(array_equal(x, prof.xobs))
        self.assertTrue(array_equal(y, prof.yobs))
        self.assertTrue(array_equal(dy, prof.dyobs))

        # Make a profile with undefined dy
        x = arange(0, 10, 0.1)
        y = x
        dy = None

        self.profile.setObservedProfile(x, y, dy)

        self.assertTrue(array_equal(x, prof.xobs))
        self.assertTrue(array_equal(y, prof.yobs))
        self.assertTrue(array_equal(ones_like(prof.xobs), prof.dyobs))

        # Get the ranged profile to make sure its the same
        self.assertTrue(array_equal(x, prof.x))
        self.assertTrue(array_equal(y, prof.y))
        self.assertTrue(array_equal(ones_like(prof.xobs), prof.dy))

        return

    def test_set_calculation_range(self):
        """Test the set_calculation_range method."""
        x = arange(2, 9.6, 0.5)
        y = array(x)
        dy = array(x)
        prof = self.profile
        # Check call before data arrays are present
        self.assertRaises(AttributeError, prof.set_calculation_range)
        self.assertRaises(AttributeError, prof.set_calculation_range, 0)
        self.assertRaises(AttributeError, prof.set_calculation_range, 0, 5)
        self.assertRaises(
            AttributeError, prof.set_calculation_range, 0, 5, 0.2
        )
        # assign data
        prof.set_observed_profile(x, y, dy)
        # Test normal execution w/o arguments
        self.assertTrue(array_equal(x, prof.x))
        self.assertTrue(array_equal(y, prof.y))
        self.assertTrue(array_equal(dy, prof.dy))
        prof.set_calculation_range()
        self.assertTrue(array_equal(x, prof.x))
        self.assertTrue(array_equal(y, prof.y))
        self.assertTrue(array_equal(dy, prof.dy))
        # Test a lower bound < xmin
        prof.set_calculation_range(xmin=0)
        self.assertTrue(array_equal(x, prof.x))
        self.assertTrue(array_equal(y, prof.y))
        self.assertTrue(array_equal(dy, prof.dy))
        # Test an upper bound > xmax
        prof.set_calculation_range(xmax=100)
        self.assertTrue(array_equal(x, prof.x))
        self.assertTrue(array_equal(y, prof.y))
        self.assertTrue(array_equal(dy, prof.dy))
        # Test xmin > xmax
        self.assertRaises(
            ValueError, prof.set_calculation_range, xmin=10, xmax=3
        )
        # Test xmax - xmin < dx
        self.assertRaises(
            ValueError, prof.set_calculation_range, xmin=3, xmax=3.9, dx=1.0
        )
        # Test dx <= 0
        self.assertRaises(ValueError, prof.set_calculation_range, dx=0)
        self.assertRaises(ValueError, prof.set_calculation_range, dx=-0.000001)
        # using string other than 'obs'
        self.assertRaises(ValueError, prof.set_calculation_range, xmin="oobs")
        self.assertRaises(ValueError, prof.set_calculation_range, xmax="oobs")
        self.assertRaises(ValueError, prof.set_calculation_range, dx="oobs")
        # This should be alright
        prof.set_calculation_range(3, 5)
        prof.set_calculation_range(xmin="obs", xmax=7, dx=0.001)
        self.assertEqual(5001, len(prof.x))
        self.assertEqual(len(prof.x), len(prof.y))
        self.assertEqual(len(prof.x), len(prof.dy))
        # Test an internal bound
        prof.set_calculation_range(4, 7, dx="obs")
        self.assertTrue(array_equal(prof.x, arange(4, 7.1, 0.5)))
        self.assertTrue(array_equal(prof.y, arange(4, 7.1, 0.5)))
        self.assertTrue(array_equal(prof.y, arange(4, 7.1, 0.5)))
        # test setting only one of the bounds
        prof.set_calculation_range(xmin=3)
        self.assertTrue(array_equal(prof.x, arange(3, 7.1, 0.5)))
        self.assertTrue(array_equal(prof.x, prof.y))
        self.assertTrue(array_equal(prof.x, prof.dy))
        prof.set_calculation_range(xmax=5.1)
        self.assertTrue(array_equal(prof.x, arange(3, 5.1, 0.5)))
        self.assertTrue(array_equal(prof.x, prof.y))
        self.assertTrue(array_equal(prof.x, prof.dy))
        prof.set_calculation_range(dx=1)
        self.assertTrue(array_equal(prof.x, arange(3, 5.1)))
        self.assertTrue(array_equal(prof.x, prof.y))
        self.assertTrue(array_equal(prof.x, prof.dy))
        # Test a new grid
        prof.set_calculation_range(4.2, 7, 0.3)
        self.assertTrue(array_equal(prof.x, arange(4.2, 6.901, 0.3)))
        self.assertTrue(allclose(prof.x, prof.y))
        self.assertTrue(allclose(prof.x, prof.dy))
        prof.set_calculation_range(xmin=4.2, xmax=6.001)
        self.assertTrue(array_equal(prof.x, arange(4.2, 6.001, 0.3)))
        # resample on a clipped grid
        prof.set_calculation_range(dx=0.5)
        self.assertTrue(array_equal(prof.x, arange(4.5, 6.1, 0.5)))
        return

    def testSetCalculationRange(self):
        """Test the deprecated setCalculationRange method.

        Remove this test when setCalculationRange is removed in 4.0.0.
        """
        x = arange(2, 9.6, 0.5)
        y = array(x)
        dy = array(x)
        prof = self.profile
        prof.set_observed_profile(x, y, dy)
        # Test normal execution w/o arguments
        self.assertTrue(array_equal(x, prof.x))
        self.assertTrue(array_equal(y, prof.y))
        self.assertTrue(array_equal(dy, prof.dy))
        prof.setCalculationRange()
        self.assertTrue(array_equal(x, prof.x))
        self.assertTrue(array_equal(y, prof.y))
        self.assertTrue(array_equal(dy, prof.dy))
        return

    def test_set_calculation_points(self):
        """Test the set_calculation_points method."""
        prof = self.profile

        x = arange(2, 10.5, 0.5)
        y = array(x)
        dy = array(x)

        # Test without data
        xcalc = arange(3, 12.2, 0.2)
        prof.set_calculation_points(xcalc)
        self.assertTrue(array_equal(xcalc, prof.x))

        # Add the data. This should change the bounds of the calculation array.
        prof.set_observed_profile(x, y, dy)
        self.assertTrue(array_equal(arange(3, 10.1, 0.2), prof.x))

        return

    def testSetCalculationPoints(self):
        """Test the deprecated setCalculationPoints method.

        Remove this test when setCalculationPoints is removed in 4.0.0.
        """
        prof = self.profile

        x = arange(2, 10.5, 0.5)
        y = array(x)
        dy = array(x)

        # Test without data
        xcalc = arange(3, 12.2, 0.2)
        prof.setCalculationPoints(xcalc)
        self.assertTrue(array_equal(xcalc, prof.x))

        # Add the data. This should change the bounds of the calculation array.
        prof.set_observed_profile(x, y, dy)
        self.assertTrue(array_equal(arange(3, 10.1, 0.2), prof.x))

        return

    def test_savetxt(self):
        "Check the savetxt method."
        prof = self.profile
        self.assertRaises(SrFitError, prof.savetxt, "foo")
        xobs = arange(-2, 3.01, 0.25)
        yobs = xobs**2
        prof.set_observed_profile(xobs, yobs)
        prof.ycalc = yobs.copy()
        fp = io.BytesIO()
        prof.savetxt(fp)
        txt = fp.getvalue().decode()
        self.assertTrue(re.match(r"^# x +ycalc +y +dy\b", txt))
        nlines = len(txt.strip().split("\n"))
        self.assertEqual(22, nlines)
        return


def testLoadtxt(datafile):
    """Test the loadtxt method."""
    prof = Profile()
    data = datafile("testdata.txt")

    def _test(p):
        assert 1e-2 == pytest.approx(p.x[0])
        assert 1.105784e-1 == pytest.approx(p.y[0])
        assert 1.802192e-3 == pytest.approx(p.dy[0])

    # Test normal load
    prof.loadtxt(data, usecols=(0, 1, 3))
    _test(prof)

    # Test trying to not set unpack
    prof.loadtxt(data, usecols=(0, 1, 3), unpack=False)
    _test(prof)
    prof.loadtxt(data, float, "#", None, None, 0, (0, 1, 3), False)
    _test(prof)

    # Try not including dy
    prof.loadtxt(data, usecols=(0, 1))
    assert 1e-2 == pytest.approx(prof.x[0])
    assert 1.105784e-1 == pytest.approx(prof.y[0])
    assert 1 == pytest.approx(prof.dy[0])

    # Try to include too little
    with pytest.raises(ValueError):
        prof.loadtxt(data, usecols=(0,))
    return


def test_load_parsed_data(parser_datafiles):
    """Test the load_parsed_data method."""
    prof = Profile()
    parser = ProfileParser()
    datafile = parser_datafiles / "four_col.gr"
    parser.parse_file(datafile)
    prof.load_parsed_data(parser)
    expected_xobs = [1.0, 1.1, 1.2]
    expected_yobs = [2.0, 2.1, 2.2]
    expected_dyobs = [0.2, 0.4, 0.6]
    assert prof.xobs.tolist() == expected_xobs
    assert prof.yobs.tolist() == expected_yobs
    assert prof.dyobs.tolist() == expected_dyobs


if __name__ == "__main__":
    unittest.main()
