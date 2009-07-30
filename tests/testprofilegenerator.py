#!/usr/bin/env python
"""Tests for refinableobj module."""

import unittest

from numpy import arange, array_equal

from diffpy.srfit.fitbase.profilegenerator import ProfileGenerator
from diffpy.srfit.fitbase.profile import Profile
from diffpy.srfit.util.clicker import Clicker


class TestProfileGenerator(unittest.TestCase):

    def setUp(self):
        self.calc = ProfileGenerator("test")
        self.profile = Profile()
        x = arange(0, 10, 0.1)
        self.profile.setCalculationPoints(x)
        self.calc.setProfile(self.profile)
        return

    def testEval(self):
        """Test the eval method."""
        calc = self.calc
        prof = self.profile

        # Try the direct evaluation
        calc.eval()
        self.assertTrue(array_equal(prof.x, prof.ycalc))

        # Try evaluation through the generate method
        calc.generate(Clicker())
        self.assertTrue(array_equal(prof.x, prof.ycalc))
        ylit = calc.literal.getValue()
        self.assertTrue(array_equal(prof.ycalc, ylit))

        return

if __name__ == "__main__":

    unittest.main()

