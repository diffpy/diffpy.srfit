#!/usr/bin/env python
##############################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2010 Trustees of the Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
##############################################################################
"""Tests for refinableobj module."""

import unittest

from numpy import array, arange, array_equal, ones_like

from diffpy.srfit.fitbase.profile import Profile
from diffpy.srfit.tests.utils import datafile


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
        self.assertEquals(profile.meta, {})
        return

    def testSetObservedProfile(self):
        """Test the setObservedProfile method."""
        # Make a profile with defined dy

        x = arange(0, 10, 0.1)
        y = x
        dy = x

        prof = self.profile
        prof.setObservedProfile(x, y, dy)

        self.assertTrue( array_equal(x, prof.xobs) )
        self.assertTrue( array_equal(y, prof.yobs) )
        self.assertTrue( array_equal(dy, prof.dyobs) )

        # Make a profile with undefined dy
        x = arange(0, 10, 0.1)
        y = x
        dy = None

        self.profile.setObservedProfile(x, y, dy)

        self.assertTrue( array_equal(x, prof.xobs) )
        self.assertTrue( array_equal(y, prof.yobs) )
        self.assertTrue( array_equal(ones_like(prof.xobs), prof.dyobs))

        # Get the ranged profile to make sure its the same
        self.assertTrue( array_equal(x, prof.x) )
        self.assertTrue( array_equal(y, prof.y) )
        self.assertTrue( array_equal(ones_like(prof.xobs), prof.dy))

        return

    def testSetCalculationRange(self):
        """Test the setCalculationRange method."""
        x = arange(2, 10, 0.5)
        y = array(x)
        dy = array(x)

        prof = self.profile

        # Check call before data arrays are present
        self.assertRaises(AttributeError, prof.setCalculationRange)
        self.assertRaises(AttributeError, prof.setCalculationRange, 0)
        self.assertRaises(AttributeError, prof.setCalculationRange, 0,
                10)
        self.assertRaises(AttributeError, prof.setCalculationRange, 0,
                10, 0.2)

        prof.setObservedProfile(x, y, dy)

        # Test normal execution w/o arguments
        prof.setCalculationRange()
        self.assertTrue( array_equal(x, prof.x) )
        self.assertTrue( array_equal(y, prof.y) )
        self.assertTrue( array_equal(dy, prof.dy) )

        # Test a lower bound < xmin
        prof.setCalculationRange(xmin = 0)
        self.assertTrue( array_equal(x, prof.x) )
        self.assertTrue( array_equal(y, prof.y) )
        self.assertTrue( array_equal(dy, dprof.y) )

        # Test an upper bound > xmax
        prof.setCalculationRange(xmax = 100)
        prof.x, prof.y, dprof.y = prof.getRangedProfile()
        self.assertTrue( array_equal(x, prof.x) )
        self.assertTrue( array_equal(y, prof.y) )
        self.assertTrue( array_equal(dy, dprof.y) )

        # Test xmin > xmax
        self.assertRaises(ValueError, profile.setCalculationRange, xmin = 10,
                xmax = 3)

        # Test xmax - xmin < dx
        self.assertRaises(ValueError, profile.setCalculationRange, xmin = 3,
                xmax = 3 + 0.4, dx = 0.5)

        # Test dx <= 0
        self.assertRaises(ValueError, profile.setCalculationRange, dx = 0)
        self.assertRaises(ValueError, profile.setCalculationRange, dx =
                -0.000001)
        # This should be alright
        profile.setCalculationRange(dx = 0.000001)

        # Test an internal bound
        prof.setCalculationRange(4, 7)
        prof.x, prof.y, dprof.y = prof.getRangedProfile()
        self.assertTrue( array_equal(prof.x, arange(4, 7.5, 0.5) ) )
        self.assertTrue( array_equal(prof.y, arange(4, 7.5, 0.5) ) )
        self.assertTrue( array_equal(dprof.y, arange(4, 7.5, 0.5) ) )

        # Test a new grid
        prof.setCalculationRange(4, 7, 0.1)
        prof.x, prof.y, dprof.y = prof.getRangedProfile()
        self.assertTrue( array_equal(prof.x, arange(4, 7.1, 0.1) ) )
        self.assertAlmostEqual( 0, sum(prof.y- arange(4, 7.1, 0.1))**2 )
        self.assertAlmostEqual( 0, sum(dprof.y- arange(4, 7.1, 0.1))**2 )

        return

    def testSetCalculationRange(self):
        """Test the setCalculationRange method."""
        prof = self.profile

        x = arange(2, 10.5, 0.5)
        y = array(x)
        dy = array(x)

        # Test without data
        xcalc = arange(3, 12.2, 0.2)
        prof.setCalculationPoints(xcalc)
        self.assertTrue( array_equal(xcalc, prof.x) )

        # Add the data. This should change the bounds of the calculation array.
        prof.setObservedProfile(x, y, dy)
        self.assertTrue( array_equal(arange(3, 10.1, 0.2), prof.x ) )

        return

    def testLoadtxt(self):
        """Test the loadtxt method"""

        prof = self.profile
        data = datafile("testdata.txt")

        def _test(p):
            self.assertAlmostEqual(1e-2, p.x[0])
            self.assertAlmostEqual(1.105784e-1, p.y[0])
            self.assertAlmostEqual(1.802192e-3, p.dy[0])

        # Test normal load
        prof.loadtxt(data, usecols=(0,1,3))
        _test(prof)

        # Test trying to not set unpack
        prof.loadtxt(data, usecols=(0,1,3), unpack = False)
        _test(prof)
        prof.loadtxt(data, float, '#', None, None, 0, (0,1,3), False)
        _test(prof)

        # Try not including dy
        prof.loadtxt(data, usecols=(0,1))
        self.assertAlmostEqual(1e-2, prof.x[0])
        self.assertAlmostEqual(1.105784e-1, prof.y[0])
        self.assertAlmostEqual(1, prof.dy[0])

        # Try to include too little
        self.assertRaises(ValueError, prof.loadtxt, data, usecols=(0,))
        return


if __name__ == "__main__":
    unittest.main()
