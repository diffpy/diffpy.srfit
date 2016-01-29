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

from numpy import arange, dot, array_equal

from diffpy.srfit.fitbase.fitcontribution import FitContribution
from diffpy.srfit.fitbase.profilegenerator import ProfileGenerator
from diffpy.srfit.fitbase.profile import Profile
from diffpy.srfit.fitbase.parameter import Parameter


class TestContribution(unittest.TestCase):

    def setUp(self):
        self.gen = ProfileGenerator("test")
        self.profile = Profile()
        self.fitcontribution = FitContribution("test")
        return

    def testSetProfile(self):
        fc = self.fitcontribution
        profile = self.profile
        fc.setProfile(self.profile)

        self.assertTrue(fc.profile is profile)
        self.assertTrue(fc.x.par is profile.xpar)
        self.assertTrue(fc.y.par is profile.ypar)
        self.assertTrue(fc.dy.par is profile.dypar)

        self.assertTrue(fc._eq is None)
        self.assertTrue(fc._reseq is None)
        return

    def testAddProfileGenerator(self):
        fc = self.fitcontribution
        gen = self.gen
        fc.addProfileGenerator(gen, "gen")

        xobs = arange(0, 10, 0.5)
        self.assertTrue(array_equal(xobs, gen(xobs)))

        self.assertTrue(gen.profile is None)
        self.assertTrue(fc._eq is not None)
        return

    def testInteraction(self):
        """Test the interaction between the profile and profile generator."""
        fc = self.fitcontribution
        profile = self.profile
        gen = self.gen

        # Add the calculator and profile
        fc.setProfile(profile)
        fc.addProfileGenerator(gen, "I")

        # Check attributes are created
        self.assertTrue(fc.profile is profile)
        self.assertTrue(fc.x.par is profile.xpar)
        self.assertTrue(fc.y.par is profile.ypar)
        self.assertTrue(fc.dy.par is profile.dypar)
        self.assertTrue(fc._eq is not None)
        self.assertTrue(fc._reseq is not None)
        self.assertTrue(gen.profile is profile)

        # create some data
        xobs = arange(0, 10, 0.5)
        yobs = xobs
        profile.setObservedProfile(xobs, yobs)

        # Make sure this is where it's supposed to be
        self.assertTrue(gen.profile.xobs is xobs)
        self.assertTrue(gen.profile.yobs is yobs)
        self.assertTrue(array_equal(fc.x.value, xobs))
        self.assertTrue(array_equal(gen.profile.x, xobs))
        self.assertTrue(array_equal(xobs, gen(xobs)))

        # Now evaluate the profile equation
        self.assertTrue(array_equal(fc._eq(), yobs))
        # And the residual equation
        self.assertAlmostEquals(0, sum(fc._reseq()))

        return

    def testReplacements(self):
        """Test attribute integrity when objects get replaced."""
        fc = self.fitcontribution
        xobs = arange(0, 10, 0.5)
        yobs = xobs
        profile = self.profile
        profile.setObservedProfile(xobs, yobs)
        xobs2 = arange(0, 10, 0.8)
        yobs2 = 0.5*xobs2
        profile2 = Profile()
        profile2.setObservedProfile(xobs2, yobs2)
        gen = self.gen

        # Validate equations
        fc.setProfile(profile)
        fc.addProfileGenerator(gen, "I")
        self.assertTrue(array_equal(gen.value, xobs))
        self.assertTrue(array_equal(fc._eq(), xobs))
        self.assertAlmostEquals(0, sum(fc._reseq()))
        eq = fc._eq
        reseq = fc._reseq

        # Now set a different profile
        fc.setProfile(profile2)
        self.assertTrue(fc.profile is profile2)
        self.assertTrue(gen.profile is profile2)
        self.assertTrue(fc._eq is eq)
        self.assertTrue(fc._reseq is reseq)
        self.assertTrue(fc._eq._value is None)
        self.assertTrue(fc._reseq._value is None)

        # Validate equations
        self.assertTrue(array_equal(xobs2, gen.value))
        self.assertTrue(array_equal(fc._eq(), gen.value))
        return

    def testResidual(self):
        """Test the residual, which requires all other methods."""
        fc = self.fitcontribution
        profile = self.profile
        gen = self.gen

        # Add the calculator and profile
        fc.setProfile(profile)
        self.assertTrue(fc.profile is profile)
        fc.addProfileGenerator(gen, "I")
        self.assertTrue(fc._eq._value is None)
        self.assertTrue(fc._reseq._value is None)
        self.assertEquals(1, len(fc._generators))
        self.assertTrue(gen.name in fc._generators)

        # Let's create some data
        xobs = arange(0, 10, 0.5)
        yobs = xobs
        profile.setObservedProfile(xobs, yobs)

        # Check our fitting equation.
        self.assertTrue(array_equal(fc._eq(), gen(xobs)))

        # Now calculate the residual
        chiv = fc.residual()
        self.assertAlmostEqual(0, dot(chiv, chiv))

        # Now change the equation
        fc.setEquation("2*I")
        self.assertTrue(fc._eq._value is None)
        self.assertTrue(fc._reseq._value is None)
        chiv = fc.residual()
        self.assertAlmostEqual(dot(yobs, yobs), dot(chiv, chiv))

        # Try to add a parameter
        c = Parameter("c", 2)
        fc._addParameter(c)
        fc.setEquation("c*I")
        self.assertTrue(fc._eq._value is None)
        self.assertTrue(fc._reseq._value is None)
        chiv = fc.residual()
        self.assertAlmostEqual(dot(yobs, yobs), dot(chiv, chiv))

        # Try something more complex
        c.setValue(3)
        fc.setEquation("c**2*sin(I)")
        self.assertTrue(fc._eq._value is None)
        self.assertTrue(fc._reseq._value is None)
        from numpy import sin
        xobs = arange(0, 10, 0.5)
        yobs = 9*sin(xobs)
        profile.setObservedProfile(xobs, yobs)
        self.assertTrue(fc._eq._value is None)
        self.assertTrue(fc._reseq._value is None)

        chiv = fc.residual()
        self.assertAlmostEqual(0, dot(chiv, chiv))

        # Choose a new residual.
        fc.setEquation("2*I")
        fc.setResidualEquation("resv")
        chiv = fc.residual()
        self.assertAlmostEqual(sum((2*xobs-yobs)**2)/sum(yobs**2),
                dot(chiv, chiv))

        # Make a custom residual.
        fc.setResidualEquation("abs(eq-y)**0.5")
        chiv = fc.residual()
        self.assertEqual(sum(abs(2*xobs-yobs)), dot(chiv, chiv))

        return


if __name__ == "__main__":
    unittest.main()
