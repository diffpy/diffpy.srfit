#!/usr/bin/env python
"""Tests for refinableobj module."""

import unittest

from numpy import linspace, array_equal, pi, sin, dot

from diffpy.srfit.fitbase.fitrecipe import FitRecipe
from diffpy.srfit.fitbase.fitcontribution import FitContribution
from diffpy.srfit.fitbase.profile import Profile
from diffpy.srfit.fitbase.parameter import Parameter

class TestFitRecipe(unittest.TestCase):

    def setUp(self):
        self.recipe = FitRecipe("recipe")
        self.recipe.fithook.verbose = 0

        # Set up the Profile
        self.profile = Profile()
        x = linspace(0, pi, 10)
        y = sin(x)
        self.profile.setObservedProfile(x, y)

        # Set up the FitContribution
        self.fitcontribution = FitContribution("cont")
        self.fitcontribution.setProfile(self.profile)
        self.fitcontribution.setEquation("A*sin(k*x + c)")
        self.fitcontribution.A.setValue(1)
        self.fitcontribution.k.setValue(1)
        self.fitcontribution.c.setValue(0)

        self.recipe.addContribution(self.fitcontribution)
        return

    def testVars(self):
        """Test to see if variables are added and removed properly."""
        recipe = self.recipe
        con = self.fitcontribution

        recipe.addVar(con.A, 2)
        recipe.addVar(con.k, 1)
        recipe.addVar(con.c, 0)
        recipe.newVar("B", 0)

        names = recipe.getNames()
        self.assertEquals(names, ["A", "k", "c", "B"])
        values = recipe.getValues()
        self.assertEquals(values, [2, 1, 0, 0])

        # Constrain a parameter to the B-variable to give it a value
        p = Parameter("Bpar", -1)
        recipe.constrain(recipe.B, p)
        values = recipe.getValues()
        self.assertEquals(values, [2, 1, 0])
        recipe.delVar(recipe.B)

        recipe.fixVar(recipe.k)

        names = recipe.getNames()
        self.assertEquals(names, ["A", "c"])
        values = recipe.getValues()
        self.assertEquals(values, [2, 0])

        recipe.fixAll()
        names = recipe.getNames()
        self.assertEquals(names, [])
        values = recipe.getValues()
        self.assertEquals(values, [])

        # The order is no longer valid
        recipe.freeAll()
        names = recipe.getNames()
        self.assertEquals(3, len(names))
        self.assertTrue("A" in names)
        self.assertTrue("k" in names)
        self.assertTrue("c" in names)
        values = recipe.getValues()
        self.assertEquals(3, len(values))
        self.assertTrue(0 in values)
        self.assertTrue(1 in values)
        self.assertTrue(2 in values)
        return

    def testResidual(self):
        """Test the residual and everything that can change it."""

        # With thing set up as they are, the residual should be 0
        res = self.recipe.residual()
        self.assertAlmostEquals(0, dot(res, res))

        # Change the c value to 1 so that the equation evaluates as sin(x+1)
        x = self.profile.x
        y = sin(x+1)
        self.recipe.cont.c.setValue(1)
        res = self.recipe.residual()
        self.assertTrue( array_equal(y-self.profile.y, res) )

        # Try some constraints
        # Make c = 2*A, A = Avar
        var = self.recipe.newVar("Avar")
        self.recipe.constrain(self.fitcontribution.c, "2*A",
                {"A" : self.fitcontribution.A})
        self.assertEquals(2, self.fitcontribution.c.value)
        self.recipe.constrain(self.fitcontribution.A, var)
        self.assertEquals(1, var.getValue())
        self.assertEquals(self.recipe.cont.A.getValue(), var.getValue())
        # c is constrained to a constrained parameter.
        self.assertEquals(2, self.fitcontribution.c.value)
        # The equation should evaluate to sin(x+2)
        x = self.profile.x
        y = sin(x+2)
        res = self.recipe.residual()
        self.assertTrue( array_equal(y-self.profile.y, res) )

        # Now try some restraints. We want c to be exactly zero. It should give
        # a penalty of (c-0)**2, which is 4 in this case
        r1 = self.recipe.restrain(self.fitcontribution.c, 0, 0, 1, 2)
        self.recipe._ready = False
        res = self.recipe.residual()
        chi2 = 4 + dot(y - self.profile.y, y - self.profile.y)
        self.assertAlmostEqual(chi2, dot(res, res) )

        # Clear the constraint and restore the value of c to 0. This should
        # give us chi2 = 0 again.
        self.recipe.unconstrain(self.fitcontribution.c)
        self.fitcontribution.c.setValue(0)
        res = self.recipe.residual([self.recipe.cont.A.getValue()])
        chi2 = 0
        self.assertAlmostEqual(chi2, dot(res, res) )

        # Remove the restraint and variable
        self.recipe.unrestrain(r1)
        self.recipe.delVar(self.recipe.Avar)
        self.recipe._ready = False
        res = self.recipe.residual()
        chi2 = 0
        self.assertAlmostEqual(chi2, dot(res, res) )

        # Add constraints at the fitcontribution level. 
        self.fitcontribution.constrain(self.fitcontribution.c, "2*A")
        # This should evaluate to sin(x+2)
        x = self.profile.x
        y = sin(x+2)
        res = self.recipe.residual()
        self.assertTrue( array_equal(y-self.profile.y, res) )

        # Add a restraint at the fitcontribution level.
        r1 = self.fitcontribution.restrain(self.fitcontribution.c, 0, 0, 1, 2)
        self.recipe._ready = False
        # The chi2 is the same as above, plus 4
        res = self.recipe.residual()
        x = self.profile.x
        y = sin(x+2)
        chi2 = 4 + dot(y - self.profile.y, y - self.profile.y)
        self.assertAlmostEqual(chi2, dot(res, res) )

        # Remove those
        self.fitcontribution.unrestrain(r1)
        self.recipe._ready = False
        self.fitcontribution.unconstrain(self.fitcontribution.c)
        self.fitcontribution.c.setValue(0)
        res = self.recipe.residual()
        chi2 = 0
        self.assertAlmostEqual(chi2, dot(res, res) )

        # Now try to use the observed profile inside of the equation
        # Set the equation equal to the data
        self.fitcontribution.setEquation("y")
        res = self.recipe.residual()
        self.assertAlmostEquals(0, dot(res, res))

        # Now add the uncertainty. This should give dy/dy = 1 for the residual
        self.fitcontribution.setEquation("y+dy")
        res = self.recipe.residual()
        self.assertAlmostEquals(len(res), dot(res, res))

        return

if __name__ == "__main__":

    unittest.main()

