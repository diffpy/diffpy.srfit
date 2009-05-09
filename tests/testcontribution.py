#!/usr/bin/env python
"""Tests for refinableobj module."""

import unittest

from numpy import arange, dot

from diffpy.srfit.fitbase.contribution import Contribution
from diffpy.srfit.fitbase.calculator import Calculator
from diffpy.srfit.fitbase.profile import Profile
from diffpy.srfit.fitbase.parameter import Parameter


class TestContribution(unittest.TestCase):

    def setUp(self):
        self.calc = Calculator("test")
        self.profile = Profile()
        self.contribution = Contribution("test")
        return

    def testResidual(self):
        """Test the residual, which requires all other methods."""

        # Add the calculator and profile
        self.contribution.setCalculator(self.calc, "I")
        self.contribution.setProfile(self.profile)

        # Let's create some data
        xobs = arange(0, 10, 0.5)
        yobs = xobs
        self.profile.setObservedProfile(xobs, yobs)

        # Now calculate the residual
        chiv = self.contribution.residual()
        self.assertAlmostEqual(0, dot(chiv, chiv))

        # Now change the equation
        self.contribution.setEquation("2*I")
        chiv = self.contribution.residual()
        self.assertAlmostEqual(dot(yobs, yobs), dot(chiv, chiv))

        # Try to add a parameter
        c = Parameter("c", 2)
        self.contribution._addParameter(c)
        self.contribution.setEquation("c*I")
        chiv = self.contribution.residual()
        self.assertAlmostEqual(dot(yobs, yobs), dot(chiv, chiv))

        # Try something more complex
        c.setValue(3)
        self.contribution.setEquation("c**2*sin(I)")
        from numpy import sin
        xobs = arange(0, 10, 0.5)
        from numpy import sin
        yobs = 9*sin(xobs)
        self.profile.setObservedProfile(xobs, yobs)
        
        chiv = self.contribution.residual()
        self.assertAlmostEqual(0, dot(chiv, chiv))

        # Choose a new residual. Note that the calculator can be treated like a
        # function, hence the I() below
        self.contribution.setEquation("2*I()")
        self.contribution.setResidualEquation("resv")
        chiv = self.contribution.residual()
        self.assertEqual(sum((2*xobs-yobs)**2)/sum(yobs**2), dot(chiv, chiv))

        # Make a custom residual.
        self.contribution.setResidualEquation("abs(eq()-y)**0.5")
        chiv = self.contribution.residual()
        self.assertEqual(sum(abs(2*xobs-yobs)), dot(chiv, chiv))

        return

if __name__ == "__main__":

    unittest.main()

