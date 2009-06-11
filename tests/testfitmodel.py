#!/usr/bin/env python
"""Tests for refinableobj module."""

import unittest

from numpy import linspace, array_equal, pi, sin, dot

from diffpy.srfit.fitbase.fitmodel import FitModel
from diffpy.srfit.fitbase.contribution import Contribution
from diffpy.srfit.fitbase.profile import Profile
from diffpy.srfit.fitbase.parameter import Parameter

class TestFitModel(unittest.TestCase):

    def setUp(self):
        self.model = FitModel("model")
        self.model.fithook.verbose = 0

        # Set up the Profile
        self.profile = Profile()
        x = linspace(0, pi, 10)
        y = sin(x)
        self.profile.setObservedProfile(x, y)

        # Set up the Contribution
        self.contribution = Contribution("cont")
        self.contribution.setProfile(self.profile)
        self.contribution.setEquation("A*sin(k*x + c)")
        self.contribution.A.setValue(1)
        self.contribution.k.setValue(1)
        self.contribution.c.setValue(0)

        self.model.addContribution(self.contribution)
        return

    def testVars(self):
        """Test to see if variables are added and removed properly."""
        model = self.model
        con = self.contribution

        model.addVar(con.A, 2)
        model.addVar(con.k, 1)
        model.addVar(con.c, 0)

        names = model.getNames()
        self.assertEquals(names, ["A", "k", "c"])
        values = model.getValues()
        self.assertEquals(values, [2, 1, 0])

        model.fixVar(model.k)

        # Try to fix a variable that is not there
        self.assertRaises(ValueError, model.fixVar, "")
        self.assertRaises(ValueError, model.fixVar, "k")
        self.assertRaises(ValueError, model.fixVar, None)
        self.assertRaises(ValueError, model.fixVar, con.A)
        names = model.getNames()
        self.assertEquals(names, ["A", "c"])
        values = model.getValues()
        self.assertEquals(values, [2, 0])

        model.fixAll()
        names = model.getNames()
        self.assertEquals(names, [])
        values = model.getValues()
        self.assertEquals(values, [])

        # The order is no longer valid
        model.freeAll()
        names = model.getNames()
        self.assertEquals(3, len(names))
        self.assertTrue("A" in names)
        self.assertTrue("k" in names)
        self.assertTrue("c" in names)
        values = model.getValues()
        self.assertEquals(3, len(values))
        self.assertTrue(0 in values)
        self.assertTrue(1 in values)
        self.assertTrue(2 in values)
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
        self.model.delVar(self.model.A)
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

