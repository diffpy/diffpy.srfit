#!/usr/bin/env python
"""Tests for refinableobj module."""

import unittest

from diffpy.srfit.fitbase.recipeorganizer import RecipeOrganizer
from diffpy.srfit.fitbase.recipeorganizer import equationFromString
from diffpy.srfit.fitbase.calculator import Calculator
from diffpy.srfit.fitbase.parameter import Parameter
from diffpy.srfit.equation.builder import EquationFactory

import numpy

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
        self.assertTrue(p1 in eq.args)
        self.assertTrue(p2 in eq.args)
        self.assertEqual(3, eq())

        # Try to use a parameter that is not registered
        self.assertRaises(ValueError, equationFromString, "p1+p2+p3", factory)

        # Pass that argument in the ns dictionary
        eq = equationFromString("p1+p2+p3", factory, {"p3":p3})
        self.assertEqual(3, len(eq.args))
        self.assertTrue(p1 in eq.args)
        self.assertTrue(p2 in eq.args)
        self.assertTrue(p3 in eq.args)
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


class TestRecipeOrganizer(unittest.TestCase):

    def setUp(self):
        self.m = RecipeOrganizer("test")
        return

    def testAddParameter(self):
        """Test the AddParameter method."""

        p1 = Parameter("p1", 1)
        p2 = Parameter("p1", 2)

        # Check normal insert
        self.m._addParameter(p1)
        self.assertTrue(self.m.p1 is p1)

        # Try to insert another parameter with the same name
        self.assertRaises(ValueError, self.m._addParameter, p2)

        return

    def testAddOrganizer(self):
        """Test the AddOrganizer method."""
        m2 = RecipeOrganizer("m2")
        p1 = Parameter("m2", 1)

        self.m._addOrganizer(m2)
        self.assertTrue(self.m.m2 is m2)

        self.assertRaises(ValueError, self.m._addOrganizer, p1)

        p1.name = "p1"
        m2._addParameter(p1)

        self.assertTrue(self.m.m2.p1 is p1)

        return

    def testClickers(self):
        """Test make sure that objects are observed by the organizer."""
        m2 = RecipeOrganizer("m2")
        p1 = Parameter("p1", 1)
        p2 = Parameter("p2", 1)
        m2._addParameter(p1)
        self.m._addOrganizer(m2)
        self.m._addParameter(p2)

        self.assertTrue(self.m.clicker >= p1.clicker)
        self.assertTrue(self.m.clicker >= m2.clicker)

        p1.setValue(1.234)
        self.assertTrue(self.m.clicker >= p1.clicker)
        self.assertTrue(self.m.clicker >= m2.clicker)

        p2.setValue(5.678)
        self.assertTrue(self.m.clicker >= p2.clicker)
        self.assertTrue(self.m.clicker > m2.clicker)

        return

    def testConstrain(self):
        """Test the constrain method."""

        p1 = Parameter("p1", 1)
        p2 = Parameter("p2", 2)
        p3 = Parameter("p3", 3)
        self.m._eqfactory.registerArgument("p1", p1)
        self.m._eqfactory.registerArgument("p2", p2)

        self.assertEquals(0, len(self.m._constraints))
        self.m.constrain(p1, "2*p2")

        self.assertTrue(p1 in self.m._constraints)
        self.assertEquals(1, len(self.m._constraints))

        p2.setValue(10)
        self.m._constraints[p1].update()
        self.assertEquals(20, p1.getValue())

        # Check errors on unregistered parameters
        self.assertRaises(ValueError, self.m.constrain, p1, "2*p3")
        self.assertRaises(ValueError, self.m.constrain, p1, "2*p2", {"p2":p3})

        # Remove the constraint
        self.m.unconstrain(p1)
        self.assertEquals(0, len(self.m._constraints))

        # Try an straight constraint
        self.m.constrain(p1, p2)
        p2.setValue(7)
        self.m._constraints[p1].update()
        self.assertEquals(7, p1.getValue())
        return

    def testRestrain(self):
        """Test the restrain method."""

        p1 = Parameter("p1", 1)
        p2 = Parameter("p2", 2)
        p3 = Parameter("p3", 3)
        self.m._eqfactory.registerArgument("p1", p1)
        self.m._eqfactory.registerArgument("p2", p2)

        self.assertEquals(0, len(self.m._restraints))
        r = self.m.restrain("p1+p2", ub = 10)
        self.assertEquals(1, len(self.m._restraints))
        p2.setValue(10)
        self.assertEquals(1, r.penalty())
        self.m.unrestrain(r)
        self.assertEquals(0, len(self.m._restraints))

        r = self.m.restrain(p1, ub = 10)
        self.assertEquals(1, len(self.m._restraints))
        p1.setValue(11)
        self.assertEquals(1, r.penalty())

        # Check errors on unregistered parameters
        self.assertRaises(ValueError, self.m.restrain, "2*p3")
        self.assertRaises(ValueError, self.m.restrain, "2*p2", ns = {"p2":p3})
        return

    def testGetConstraints(self):
        """Test the _getConstraints method."""
        m2 = RecipeOrganizer("m2")
        self.m._addOrganizer(m2)

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

        cons = self.m._getConstraints()
        self.assertTrue(p1 in cons)
        self.assertTrue(p3 in cons)
        self.assertEquals(2, len(cons))
        return

    def testGetRestraints(self):
        """Test the _getRestraints method."""
        m2 = RecipeOrganizer("m2")
        self.m._addOrganizer(m2)

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

        res = self.m._getRestraints()
        self.assertTrue(r1 in res)
        self.assertTrue(r2 in res)
        self.assertEquals(2, len(res))
        return

    def testLocateParameter(self):
        """Test the AddOrganizer method."""
        m1 = self.m
        p1 = Parameter("p1", 1)
        m1._addParameter(p1)

        m2 = RecipeOrganizer("m2")
        p2 = Parameter("p2", 2)
        m2._addParameter(p2)

        m1._addOrganizer(m2)

        p3 = Parameter("p3", 3)

        # Locate m2 in m1 (m1.m2)
        loc = m1._locateChild(m2)
        self.assertEquals(loc, [m1, m2])

        # Locate p1 (m1.p1)
        loc = m1._locateChild(p1)
        self.assertEquals(loc, [m1, p1])

        # Locate p2 in m2 (m2.p2)
        loc = m2._locateChild(p2)
        self.assertEquals(loc, [m2, p2])

        # Locate p2 in m1 (m1.m2.p2)
        loc = m1._locateChild(p2)
        self.assertEquals(loc, [m1, m2, p2])

        # Locate p3 in m1 (not there)
        loc = m1._locateChild(p3)
        self.assertEquals(loc, [])

        # Locate p3 in m2 (not there)
        loc = m2._locateChild(p3)
        self.assertEquals(loc, [])

        return

    def testRegisterCalculator(self):

        class GCalc(Calculator):

            def __init__(self, name):
                Calculator.__init__(self, name)
                self.newParameter("A", 1.0)
                self.newParameter("center", 0.0)
                self.newParameter("width", 0.1)
                return

            def __call__(self, x):
                A = self.A.getValue()
                c = self.center.getValue()
                w = self.width.getValue()
                return A * numpy.exp(-0.5*((x-c)/w)**2)

        # End class GCalc

        g = GCalc("g")

        self.m.registerCalculator(g)

        x = numpy.arange(0.5, 10, 0.5)
        self.m.x.setValue(x)

        self.m.g.center.setValue(3.0)

        self.assertTrue(numpy.array_equal(numpy.exp(-0.5*((x-3.0)/0.1)**2),
            g(x)))

        # Use this in another equation

        eq = self.m.registerStringFunction("g/x - 1", "pdf")
        self.assertTrue(numpy.array_equal(g(x)/x - 1, eq()))

        # Now swap out what g is
        g2 = self.m.registerStringFunction("x", "g")
        self.assertTrue(numpy.array_equal(numpy.zeros_like(x), eq()))

        return

    def testRegisterFunction(self):
        """Test registering various functions."""
        def g(A, c, w, x):
            return A * numpy.exp(-0.5*((x-c)/w)**2)

        eq = self.m.registerFunction(g)

        x = numpy.arange(0.5, 10, 0.5)
        self.m.x.setValue(x)
        self.m.A.setValue(1.0)
        self.m.c.setValue(3.0)
        self.m.w.setValue(0.1)

        self.assertTrue(numpy.array_equal(numpy.exp(-0.5*((x-3.0)/0.1)**2),
            eq()))

        # Use this in another equation

        eq2 = self.m.registerStringFunction("g/x - 1", "pdf")
        self.assertTrue(numpy.array_equal(eq(x)/x - 1, eq2()))

        # Now swap out what g is
        g2 = self.m.registerStringFunction("x", "g")
        self.assertTrue(numpy.array_equal(numpy.zeros_like(x), eq2()))
        return

    def testRegisterStringFunction(self):
        """Test registering string functions in various ways."""

        # Make an equation.
        eq1 = self.m.registerStringFunction("x**2 + 3", "eq1")
        eq1.x.setValue(0)

        # Add a parameter
        self.m._newParameter("y", 3.0)

        # Make sure that x and y are in the organizer
        self.assertEquals(0, self.m.x.getValue())
        self.assertEquals(3.0, self.m.y.getValue())

        # Use eq1 in some equations

        # x**2 (x**2 + 3 - 3)
        eq2 = self.m.registerStringFunction("eq1 - 3", "eq2")
        # y**2
        eq3 = self.m.registerStringFunction("eq1(y) - 3", "eq3")

        # Test these equations.
        self.assertEquals(0, eq2())
        self.assertEquals(9.0, eq3())
        # Note that eq3 injects the value of y into the argument of eq1, which
        # is x. Thus, calling eq3 sets x to 3.
        self.assertEquals(9.0, eq2())

        # One more level of embedding
        # 2*y**2
        eq4 = self.m.registerStringFunction("2*eq3", "eq4")
        self.assertEquals(18.0, eq4())

        # Replace eq1 with some other equation
        eq1b = self.m.registerStringFunction("x + 3", "eq1")

        # Make sure that instances of eq1 are replaced in eq2, eq3 and eq4.
        # x
        self.assertEquals(3, eq2())
        # y
        self.assertEquals(3, eq3())
        # 2*y
        self.assertEquals(6, eq4())

        return
        

#

if __name__ == "__main__":

    unittest.main()

