#!/usr/bin/env python
"""Tests for refinableobj module."""

import unittest

from numpy import arange, dot

from diffpy.srfit.fitbase.fitcontribution import FitContribution
from diffpy.srfit.fitbase.profilegenerator import ProfileGenerator
from diffpy.srfit.fitbase.profile import Profile
from diffpy.srfit.fitbase.parameter import Parameter


class TestContribution(unittest.TestCase):

    def setUp(self):
        self.calc = ProfileGenerator("test")
        self.profile = Profile()
        self.fitcontribution = FitContribution("test")
        return

    def testResidual(self):
        """Test the residual, which requires all other methods."""

        # Add the calculator and profile
        self.fitcontribution.addProfileGenerator(self.calc, "I")
        self.fitcontribution.setProfile(self.profile)

        # Let's create some data
        xobs = arange(0, 10, 0.5)
        yobs = xobs
        self.profile.setObservedProfile(xobs, yobs)

        # Now calculate the residual
        chiv = self.fitcontribution.residual()
        self.assertAlmostEqual(0, dot(chiv, chiv))

        # Now change the equation
        self.fitcontribution.setEquation("2*I")
        chiv = self.fitcontribution.residual()
        self.assertAlmostEqual(dot(yobs, yobs), dot(chiv, chiv))

        # Try to add a parameter
        c = Parameter("c", 2)
        self.fitcontribution._addParameter(c)
        self.fitcontribution.setEquation("c*I")
        chiv = self.fitcontribution.residual()
        self.assertAlmostEqual(dot(yobs, yobs), dot(chiv, chiv))

        # Try something more complex
        c.setValue(3)
        self.fitcontribution.setEquation("c**2*sin(I)")
        from numpy import sin
        xobs = arange(0, 10, 0.5)
        from numpy import sin
        yobs = 9*sin(xobs)
        self.profile.setObservedProfile(xobs, yobs)
        
        chiv = self.fitcontribution.residual()
        self.assertAlmostEqual(0, dot(chiv, chiv))

        # Choose a new residual. Note that the calculator can be treated like a
        # function, hence the I() below
        self.fitcontribution.setEquation("2*I")
        self.fitcontribution.setResidualEquation("resv")
        chiv = self.fitcontribution.residual()
        self.assertAlmostEqual(sum((2*xobs-yobs)**2)/sum(yobs**2), dot(chiv, chiv))

        # Make a custom residual.
        self.fitcontribution.setResidualEquation("abs(eq-y)**0.5")
        chiv = self.fitcontribution.residual()
        self.assertEqual(sum(abs(2*xobs-yobs)), dot(chiv, chiv))

        return

if __name__ == "__main__":

    unittest.main()

