#!/usr/bin/env python
"""Tests for refinableobj module."""

import unittest

from numpy import linspace, array_equal, pi, sin, dot

from diffpy.srfit.fitbase.fitmodel import FitModel
from diffpy.srfit.fitbase.contribution import Contribution
from diffpy.srfit.fitbase.calculator import Calculator
from diffpy.srfit.fitbase.profile import Profile
from diffpy.srfit.fitbase.parameter import Parameter

class TestFitModel(unittest.TestCase):

    def setUp(self):
        self.model = FitModel("model")
        # Set up the calculator
        self.calc = Calculator("calc")

        # Set up the Profile
        self.profile = Profile()
        x = linspace(0, pi, 10)
        y = sin(x)
        self.profile.setObservedProfile(x, y)

        # Set up the Contribution
        self.contribution = Contribution("cont")
        self.contribution._newParameter("A", 1)
        self.contribution._newParameter("k", 1)
        self.contribution._newParameter("c", 0)
        self.contribution.setCalculator(self.calc, "x")
        self.contribution.setEquation("A*sin(k*x + c)")
        self.contribution.setProfile(self.profile, yname = "y", dyname = "dy")

        self.model.addContribution(self.contribution)
        return

    def testResidual(self):
        """Test the residual and everything that can change it."""

        # With thing set up as they are, the residual should be 0
        res = self.model.residual()
        self.assertAlmostEquals(0, dot(res, res))

        # Change the c value to 1 so that the equation evaluates as sin(x+1)
        x = self.profile.x
        y = sin(x+1)
        self.model.cont.c.setValue(1)
        res = self.model.residual()
        self.assertTrue( array_equal(y-self.profile.y, res) )

        # Try some constraints
        # Make c = 2*A
        self.model.addVar(self.model.cont.A)
        self.model.constrain(self.contribution.c, "2*A")
        # This should evaluate to sin(x+2)
        x = self.profile.x
        y = sin(x+2)
        res = self.model.residual([self.model.cont.A.getValue()])
        self.assertTrue( array_equal(y-self.profile.y, res) )

        # Now try some restraints. We want c to be exactly zero. It should give
        # a penalty of (c-0)**2, which is 4 in this case
        r1 = self.model.restrain(self.contribution.c, 0, 0, 1, 2)
        res = self.model.residual([self.model.cont.A.getValue()])
        chi2 = 4 + dot(y - self.profile.y, y - self.profile.y)
        self.assertAlmostEqual(chi2, dot(res, res) )

        # Clear the constraint and restore the value of c to 0. This should
        # give us chi2 = 0 again.
        self.model.unconstrain(self.contribution.c)
        self.contribution.c.setValue(0)
        res = self.model.residual([self.model.cont.A.getValue()])
        chi2 = 0
        self.assertAlmostEqual(chi2, dot(res, res) )

        # Remove the restraint and variable
        self.model.unrestrain(r1)
        self.model.delVar(self.contribution.A)
        res = self.model.residual([])
        chi2 = 0
        self.assertAlmostEqual(chi2, dot(res, res) )

        # Add constraints at the contribution level. This would normally be
        # done before handing the contribution to a FitModel, so we must call
        # _prepare() manually.
        self.contribution.constrain(self.contribution.c, "2*A")
        self.model._prepare()
        # This should evaluate to sin(x+2)
        x = self.profile.x
        y = sin(x+2)
        res = self.model.residual([])
        self.assertTrue( array_equal(y-self.profile.y, res) )

        # Add a restraint at the contribution level. This would normally be
        # done before handing the contribution to a FitModel, so we must call
        # _prepare() manually.
        r1 = self.contribution.restrain(self.contribution.c, 0, 0, 1, 2)
        self.model._prepare()
        # The chi2 is the same as above, plus 4
        res = self.model.residual([])
        x = self.profile.x
        y = sin(x+2)
        chi2 = 4 + dot(y - self.profile.y, y - self.profile.y)
        self.assertAlmostEqual(chi2, dot(res, res) )

        # Remove those
        self.contribution.unrestrain(r1)
        self.contribution.unconstrain(self.contribution.c)
        self.contribution.c.setValue(0)
        self.model._prepare()
        res = self.model.residual([])
        chi2 = 0
        self.assertAlmostEqual(chi2, dot(res, res) )

        # Now try to use the observed profile inside of the equation
        # Set the equation equal to the data
        self.contribution.setEquation("y")
        res = self.model.residual([])
        self.assertAlmostEquals(0, dot(res, res))

        # Now add the uncertainty. This should give dy/dy = 1 for the residual
        self.contribution.setEquation("y+dy")
        res = self.model.residual([])
        self.assertAlmostEquals(len(res), dot(res, res))

        return


if __name__ == "__main__":

    unittest.main()

