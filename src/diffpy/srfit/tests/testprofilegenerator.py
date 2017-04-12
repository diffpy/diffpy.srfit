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

"""Tests for profilegenerator module."""

import pickle
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
        val = gen.operation()
        self.assertTrue(array_equal(prof.x, val))

        # Try evaluation through __call__
        val = gen(2 * prof.x)
        self.assertTrue(array_equal(2 * prof.x, val))
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
        self.assertTrue(array_equal(x, val))

        # Verify generated value listens to changes in profile.x.
        x3 = x + 3
        prof.x = x3
        self.assertTrue(array_equal(x3, gen.value))

        # Make sure attributes get updated with a new profile.
        x = arange(0, 8, 0.1)
        prof = Profile()
        prof.setCalculationPoints(x)
        gen.setProfile(prof)
        self.assertTrue(gen._value is None)
        self.assertTrue(array_equal(x, gen.value))
        return


    def test_pickling(self):
        """Test pickling of ProfileGenerator.
        """
        data = pickle.dumps(self.gen)
        gen2 = pickle.loads(data)
        self.assertEqual('test', gen2.name)
        x = self.profile.x
        self.assertTrue(array_equal(x, gen2.operation()))
        self.assertTrue(array_equal(3 * x, gen2(3 * x)))
        return

# End of class TestProfileGenerator

if __name__ == "__main__":
    unittest.main()
