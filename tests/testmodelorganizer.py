#!/usr/bin/env python
"""Tests for refinableobj module."""

import unittest

from diffpy.srfit.fitbase.modelorganizer import ModelOrganizer
from diffpy.srfit.fitbase.modelorganizer import equationFromString
from diffpy.srfit.fitbase.parameter import Parameter
from diffpy.srfit.equation.builder import EquationFactory



class TestEquationFromString(unittest.TestCase):

    def testEquationFromString(self):
        """Test the equationFromString method."""

        p1 = Parameter("p1", 1)
        p2 = Parameter("p2", 2)
        p3 = Parameter("p3", 3)
        p4 = Parameter("p4", 4)

        factory = EquationFactory()

        factory.registerArgument("p1", p1)
        factory.registerArgument("p2", p2)

        # Check usage where all parameters are registered with the factory
        eq = equationFromString("p1+p2", factory)

        self.assertEqual(2, len(eq.args))
        self.assertTrue(p1 in eq.arglist)
        self.assertTrue(p2 in eq.arglist)
        self.assertEqual(3, eq())

        # Try to use a parameter that is not registered
        self.assertRaises(ValueError, equationFromString, "p1+p2+p3", factory)

        # Pass that argument in the ns dictionary
        eq = equationFromString("p1+p2+p3", factory, {"p3":p3})
        self.assertEqual(3, len(eq.args))
        self.assertTrue(p1 in eq.arglist)
        self.assertTrue(p2 in eq.arglist)
        self.assertTrue(p3 in eq.arglist)
        self.assertEqual(6, eq())

        # Make sure that there are no remnants of p3 in the factory
        self.assertTrue("p3" not in factory.builders)

        # Pass and use an unregistered parameter
        self.assertRaises(ValueError, equationFromString, "p1+p2+p3+p4", 
                factory, {"p3":p3})

        # Try to overload a registered parameter
        self.assertRaises(ValueError, equationFromString, "p1+p2",
                factory, {"p2":p3})

        return


class TestModelOrganizer(unittest.TestCase):

    def setUp(self):
        self.m = ModelOrganizer()
        return

    def testConstrain(self):
        """Test the constrain method."""

        p1 = Parameter("p1", 1)
        p2 = Parameter("p2", 2)
        p3 = Parameter("p3", 3)
        self.m._eqfactory.registerArgument("p1", p1)
        self.m._eqfactory.registerArgument("p2", p2)

        self.assertEquals(0, len(self.m.constraints))
        self.m.constrain(p1, "2*p2")

        self.assertTrue(p1 in self.m.constraints)
        self.assertEquals(1, len(self.m.constraints))

        p2.setValue(10)
        self.m.constraints[p1].update()
        self.assertEquals(20, p1.getValue())

        # Check errors on unregistered parameters
        self.assertRaises(ValueError, self.m.constrain, p1, "2*p3")
        self.assertRaises(ValueError, self.m.constrain, p1, "2*p2", {"p2":p3})

        # Remove the constraint
        self.m.unconstrain(p1)
        self.assertEquals(0, len(self.m.constraints))
        return

    def testRestrain(self):
        """Test the restrain method."""

        p1 = Parameter("p1", 1)
        p2 = Parameter("p2", 2)
        p3 = Parameter("p3", 3)
        self.m._eqfactory.registerArgument("p1", p1)
        self.m._eqfactory.registerArgument("p2", p2)

        self.assertEquals(0, len(self.m.restraints))
        r = self.m.restrain("p1+p2", ub = 10)
        self.assertEquals(1, len(self.m.restraints))
        p2.setValue(10)
        self.assertEquals(1, r.penalty())
        self.m.unrestrain(r)
        self.assertEquals(0, len(self.m.restraints))

        r = self.m.restrain(p1, ub = 10)
        self.assertEquals(1, len(self.m.restraints))
        p1.setValue(11)
        self.assertEquals(1, r.penalty())

        # Check errors on unregistered parameters
        self.assertRaises(ValueError, self.m.restrain, "2*p3")
        self.assertRaises(ValueError, self.m.restrain, "2*p2", ns = {"p2":p3})
        return

    def testGetConstraints(self):
        """Test the getConstraints method."""
        m2 = ModelOrganizer()
        self.m.suborganizers.append(m2)

        p1 = Parameter("p1", 1)
        p2 = Parameter("p2", 2)
        p3 = Parameter("p3", 3)
        p4 = Parameter("p4", 4)

        self.m._eqfactory.registerArgument("p1", p1)
        self.m._eqfactory.registerArgument("p2", p2)

        m2._eqfactory.registerArgument("p3", p3)
        m2._eqfactory.registerArgument("p4", p4)

        self.m.constrain(p1, "p2")
        m2.constrain(p3, "p4")

        cons = self.m.getConstraints()
        self.assertTrue(p1 in cons)
        self.assertTrue(p3 in cons)
        self.assertEquals(2, len(cons))
        return

    def testGetRestraints(self):
        """Test the getConstraints method."""
        m2 = ModelOrganizer()
        self.m.suborganizers.append(m2)

        p1 = Parameter("p1", 1)
        p2 = Parameter("p2", 2)
        p3 = Parameter("p3", 3)
        p4 = Parameter("p4", 4)

        self.m._eqfactory.registerArgument("p1", p1)
        self.m._eqfactory.registerArgument("p2", p2)

        m2._eqfactory.registerArgument("p3", p3)
        m2._eqfactory.registerArgument("p4", p4)

        r1 = self.m.restrain("p1 + p2")
        r2 = m2.restrain("2*p3 + p4")

        res = self.m.getRestraints()
        self.assertTrue(r1 in res)
        self.assertTrue(r2 in res)
        self.assertEquals(2, len(res))
        return

if __name__ == "__main__":

    unittest.main()

