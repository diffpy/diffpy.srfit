#!/usr/bin/env python
"""Tests for refinableobj module."""

import unittest

from numpy import arange, array_equal

from diffpy.srfit.fitbase.profilegenerator import ProfileGenerator
from diffpy.srfit.fitbase.profile import Profile

class TestProfileGenerator(unittest.TestCase):

    def setUp(self):
        self.gen = ProfileGenerator("test")
        self.profile = Profile()
        x = arange(0, 10, 0.1)
        self.profile.setCalculationPoints(x)
        self.gen.setProfile(self.profile)
        return

    def testOperation(self):
        """Test the operation method."""
        gen = self.gen
        prof = self.profile

        # Try the direct evaluation
        gen.operation()
        self.assertTrue(array_equal(prof.x, prof.ycalc))

        # Try evaluation through __call__
        gen(prof.x)
        self.assertTrue(array_equal(prof.x, prof.ycalc))
        return

    def testUpdate(self):
        """Update and change the profile to make sure generator is flushed."""
        gen = self.gen
        prof = self.profile

        # Make sure attributes get updated with a change in the calculation
        # points.
        x = arange(0, 9, 0.1)
        prof.setCalculationPoints(x)
        self.assertTrue(gen._value is None)
        val = gen.value
        self.assertTrue(array_equal(x, prof.ycalc))
        self.assertTrue(array_equal(prof.x, prof.ycalc))
        self.assertTrue(array_equal(val, prof.ycalc))
        self.assertTrue(array_equal(gen._value, prof.ycalc))

        # Make sure attributes get updated with a new profile.
        x = arange(0, 8, 0.1)
        prof = Profile()
        prof.setCalculationPoints(x)
        gen.setProfile(prof)
        self.assertTrue(gen._value is None)
        val = gen.value
        self.assertTrue(array_equal(x, prof.ycalc))
        self.assertTrue(array_equal(prof.x, prof.ycalc))
        self.assertTrue(array_equal(val, prof.ycalc))
        self.assertTrue(array_equal(gen._value, prof.ycalc))
        return


if __name__ == "__main__":
    unittest.main()
