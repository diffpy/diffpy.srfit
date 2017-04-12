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

from numpy import arange, dot, array_equal, sin

from diffpy.srfit.fitbase.fitcontribution import FitContribution
from diffpy.srfit.fitbase.profilegenerator import ProfileGenerator
from diffpy.srfit.fitbase.profile import Profile
from diffpy.srfit.fitbase.parameter import Parameter
from diffpy.srfit.exceptions import SrFitError
from diffpy.srfit.tests.utils import noObserversInGlobalBuilders


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
        # verify standard profile setup
        self.assertTrue(fc.profile is profile)
        self.assertTrue(fc.x.par is profile.xpar)
        self.assertTrue(fc.y.par is profile.ypar)
        self.assertTrue(fc.dy.par is profile.dypar)
        self.assertTrue(fc._eq is None)
        self.assertTrue(fc._reseq is None)
        # check type checking
        fc1 = FitContribution('test1')
        self.assertRaises(TypeError, fc1.setProfile, 'invalid')
        # check if residual equation is set up when possible
        fc2 = FitContribution('test2')
        fc2.setEquation('A * x')
        fc2.setProfile(profile)
        self.assertFalse(fc2._reseq is None)
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
        self.assertAlmostEqual(0, sum(fc._reseq()))

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
        self.assertAlmostEqual(0, sum(fc._reseq()))
        self.assertEqual(len(xobs), len(fc.residual()))
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
        self.assertEqual(len(xobs2), len(fc.residual()))
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
        self.assertEqual(1, len(fc._generators))
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

        # Test configuration checks
        fc1 = FitContribution('test1')
        self.assertRaises(SrFitError, fc1.setResidualEquation, 'chiv')
        fc1.setProfile(self.profile)
        self.assertRaises(SrFitError, fc1.setResidualEquation, 'chiv')
        fc1.setEquation('A * x')
        fc1.setResidualEquation('chiv')
        self.assertTrue(noObserversInGlobalBuilders())
        return


    def test_setEquation(self):
        """Check replacement of removed parameters."""
        fc = self.fitcontribution
        fc.setEquation("x + 5")
        fc.x.setValue(2)
        self.assertEqual(7, fc.evaluate())
        fc.removeParameter(fc.x)
        x = arange(0, 10, 0.5)
        fc.newParameter('x', x)
        self.assertTrue(array_equal(5 + x, fc.evaluate()))
        self.assertTrue(noObserversInGlobalBuilders())
        return


    def test_getEquation(self):
        """Check getting the current profile simulation formula."""
        fc = self.fitcontribution
        self.assertEqual('', fc.getEquation())
        fc.setEquation("A * sin(x + 5)")
        self.assertEqual('(A * sin((x + 5)))', fc.getEquation())
        self.assertTrue(noObserversInGlobalBuilders())
        return


    def test_getResidualEquation(self):
        """Check getting the current formula for residual equation."""
        fc = self.fitcontribution
        self.assertEqual('', fc.getResidualEquation())
        fc.setProfile(self.profile)
        fc.setEquation('A * x + B')
        self.assertEqual('((eq - y) / dy)', fc.getResidualEquation())
        fc.setResidualEquation('2 * (eq - y)')
        self.assertEqual('(2 * (eq - y))', fc.getResidualEquation())
        return


    def test_releaseOldEquations(self):
        """Ensure EquationFactory does not hold to obsolete Equations.
        """
        fc = self.fitcontribution
        self.assertEqual(0, len(fc._eqfactory.equations))
        for i in range(5):
            fc.setEquation('A * x + B')
        self.assertEqual(1, len(fc._eqfactory.equations))
        fc.setProfile(self.profile)
        for i in range(5):
            fc.setResidualEquation('chiv')
        self.assertEqual(2, len(fc._eqfactory.equations))
        return


    def test_registerFunction(self):
        """Ensure registered function works after second setEquation call.
        """
        fc = self.fitcontribution
        fsquare = lambda x : x**2
        fc.registerFunction(fsquare, name='fsquare')
        fc.setEquation('fsquare')
        fc.x.setValue(5)
        self.assertEqual(25, fc.evaluate())
        fc.x << 6
        self.assertEqual(36, fc.evaluate())
        fc.setEquation('fsquare + 5')
        self.assertEqual(41, fc.evaluate())
        fc.x << -1
        self.assertEqual(6, fc.evaluate())
        return


if __name__ == "__main__":
    unittest.main()
