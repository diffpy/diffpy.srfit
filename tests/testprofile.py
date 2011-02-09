#!/usr/bin/env python
"""Tests for refinableobj module."""

import unittest
import os.path

from numpy import array, arange, array_equal, ones_like

from diffpy.srfit.fit.profile import Profile

thisfile = locals().get('__file__', 'testprofile.py')
tests_dir = os.path.dirname(os.path.abspath(thisfile))
testdata_dir = os.path.join(tests_dir, 'testdata')

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
        self.assertEquals(profile.meta, {})
        return

    def testSetObservedProfile(self):
        """Test the setObserved method."""
        # Make a profile with defined dy

        x = arange(0, 10, 0.1)
        y = x
        dy = x

        prof = self.profile
        prof.setObserved(x, y, dy)

        self.assertTrue( array_equal(x, prof.xobs) )
        self.assertTrue( array_equal(y, prof.yobs) )
        self.assertTrue( array_equal(dy, prof.dyobs) )

        # Make a profile with undefined dy
        x = arange(0, 10, 0.1)
        y = x
        dy = None

        self.profile.setObserved(x, y, dy)

        self.assertTrue( array_equal(x, prof.xobs) )
        self.assertTrue( array_equal(y, prof.yobs) )
        self.assertTrue( array_equal(ones_like(prof.xobs), prof.dyobs))

        # Get the ranged profile to make sure its the same
        self.assertTrue( array_equal(x, prof.x) )
        self.assertTrue( array_equal(y, prof.y) )
        self.assertTrue( array_equal(ones_like(prof.xobs), prof.dy))

        return

    def testSetRange(self):
        """Test the setRange method."""
        x = arange(2, 10, 0.5)
        y = array(x)
        dy = array(x)

        prof = self.profile

        # Check call before data arrays are present
        self.assertRaises(AttributeError, prof.setRange)
        self.assertRaises(AttributeError, prof.setRange, 0)
        self.assertRaises(AttributeError, prof.setRange, 0,
                10)
        self.assertRaises(AttributeError, prof.setRange, 0,
                10, 0.2)

        prof.setObserved(x, y, dy)

        # Test normal execution w/o arguments
        prof.setRange()
        self.assertTrue( array_equal(x, prof.x) )
        self.assertTrue( array_equal(y, prof.y) )
        self.assertTrue( array_equal(dy, prof.dy) )

        # Test a lower bound < xmin
        prof.setRange(xmin = 0)
        self.assertTrue( array_equal(x, prof.x) )
        self.assertTrue( array_equal(y, prof.y) )
        self.assertTrue( array_equal(dy, prof.dy) )

        # Test an upper bound > xmax
        prof.setRange(xmax = 100)
        self.assertTrue( array_equal(x, prof.x) )
        self.assertTrue( array_equal(y, prof.y) )
        self.assertTrue( array_equal(dy, prof.dy) )

        # Test xmin > xmax
        self.assertRaises(ValueError, prof.setRange, xmin = 10,
                xmax = 3)

        # Test xmax - xmin < dx
        self.assertRaises(ValueError, prof.setRange, xmin = 3,
                xmax = 3 + 0.4, dx = 0.5)

        # Test dx <= 0
        self.assertRaises(ValueError, prof.setRange, dx = 0)
        self.assertRaises(ValueError, prof.setRange, dx = -0.000001)
        # This should be alright
        prof.setRange(dx = 0.000001)

        # Test an internal bound
        prof.setRange(4, 7)
        self.assertTrue( array_equal(prof.x, arange(4, 7.5, 0.5) ) )
        self.assertTrue( array_equal(prof.y, arange(4, 7.5, 0.5) ) )
        self.assertTrue( array_equal(prof.dy, arange(4, 7.5, 0.5) ) )

        # Test a new grid
        prof.setRange(4, 7, 0.1)
        self.assertTrue( array_equal(prof.x, arange(4, 7.1, 0.1) ) )
        self.assertAlmostEqual( 0, sum(prof.y - arange(4, 7.1, 0.1))**2 )
        self.assertAlmostEqual( 0, sum(prof.dy - arange(4, 7.1, 0.1))**2 )

        return

    def testSetObserved(self):
        """Test the setObserved method."""
        prof = self.profile

        x = arange(2, 10.5, 0.5)
        y = array(x)
        dy = array(x)

        # Test without data
        xcalc = arange(3, 12.2, 0.2)
        prof.setPoints(xcalc)
        self.assertTrue( array_equal(xcalc, prof.x) )

        # Add the data. This should change the bounds of the calculation array.
        prof.setObserved(x, y, dy)
        self.assertTrue( array_equal(arange(3, 10.1, 0.2), prof.x ) )

        return

if __name__ == "__main__":

    unittest.main()

