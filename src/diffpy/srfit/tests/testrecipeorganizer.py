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

import sys
import unittest
import cStringIO

from diffpy.srfit.equation.builder import EquationFactory
from diffpy.srfit.fitbase.calculator import Calculator
from diffpy.srfit.fitbase.parameter import Parameter
from diffpy.srfit.fitbase.recipeorganizer import equationFromString
from diffpy.srfit.fitbase.recipeorganizer import RecipeContainer
from diffpy.srfit.fitbase.recipeorganizer import RecipeOrganizer

import numpy

# ----------------------------------------------------------------------------

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
                factory, {"p2":p4})

        return

# ----------------------------------------------------------------------------

class TestRecipeContainer(unittest.TestCase):

    def setUp(self):
        self.m = RecipeContainer("test")

        # Add another managed dictionary
        self.m._containers = {}
        self.m._manage(self.m._containers)

        return

    def testAccessors(self):
        """Test accessors."""
        m1 = self.m
        p1 = Parameter("p1", 1)
        m1._addObject(p1, m1._parameters)

        m2 = RecipeContainer("m2")
        p2 = Parameter("p2", 2)
        m2._addObject(p2, m2._parameters)

        m1._addObject(m2, m1._containers)

        self.assertTrue(m1.m2 is m2)
        self.assertTrue(m1.p1 is p1)
        self.assertTrue(m2.p2 is p2)
        self.assertTrue(m1.m2.p2 is p2)

        self.assertTrue(m1[0] is p1)
        self.assertTrue(m1[0:] == [p1,])
        self.assertTrue(m2[0] is p2)

        self.assertEqual(1, len(m1))
        self.assertEqual(1, len(m2))
        return


    def testLocateManagedObject(self):
        """Test the locateManagedObject method."""
        m1 = self.m
        p1 = Parameter("p1", 1)
        m1._addObject(p1, m1._parameters)

        m2 = RecipeContainer("m2")
        p2 = Parameter("p2", 2)
        m2._addObject(p2, m2._parameters)

        m1._addObject(m2, m1._containers)

        p3 = Parameter("p3", 3)

        # Locate m2 in m1 (m1.m2)
        loc = m1._locateManagedObject(m2)
        self.assertEqual(loc, [m1, m2])

        # Locate p1 (m1.p1)
        loc = m1._locateManagedObject(p1)
        self.assertEqual(loc, [m1, p1])

        # Locate p2 in m2 (m2.p2)
        loc = m2._locateManagedObject(p2)
        self.assertEqual(loc, [m2, p2])

        # Locate p2 in m1 (m1.m2.p2)
        loc = m1._locateManagedObject(p2)
        self.assertEqual(loc, [m1, m2, p2])

        # Locate p3 in m1 (not there)
        loc = m1._locateManagedObject(p3)
        self.assertEqual(loc, [])

        # Locate p3 in m2 (not there)
        loc = m2._locateManagedObject(p3)
        self.assertEqual(loc, [])

        return

# ----------------------------------------------------------------------------

class TestRecipeOrganizer(unittest.TestCase):

    def setUp(self):
        self.m = RecipeOrganizer("test")
        # Add a managed container so we can do more in-depth tests.
        self.m._containers = {}
        self.m._manage(self.m._containers)
        return

    def tearDown(self):
        sys.stdout = sys.__stdout__
        return


    def testNewParameter(self):
        """Test the addParameter method."""

        m = self.m

        p1 = Parameter("p1", 1)
        m._addParameter(p1)

        # Test duplication of Parameters
        self.assertRaises(ValueError, m._newParameter, "p1", 0)

        # Add a new Parameter
        p2 = m._newParameter("p2", 0)
        self.assertTrue(p2 is m.p2)
        return


    def testAddParameter(self):
        """Test the addParameter method."""

        m = self.m

        p1 = Parameter("p1", 1)
        p2 = Parameter("p1", 2)

        # Check normal insert
        m._addParameter(p1)
        self.assertTrue(m.p1 is p1)
        self.assertTrue(p1.name in m._eqfactory.builders)

        # Try to insert another parameter with the same name
        self.assertRaises(ValueError, m._addParameter, p2)

        # Now allow this
        m._addParameter(p2, check=False)
        self.assertTrue(m.p1 is p2)
        self.assertTrue(p1.name in m._eqfactory.builders)

        # Try to insert a Parameter when a RecipeContainer with the same name
        # is already inside.
        c = RecipeContainer("test")
        m._addObject(c, m._containers)

        p3 = Parameter("test", 0)
        self.assertRaises(ValueError, m._addParameter, p3)

        p4 = Parameter("xyz", 0)
        m._addParameter(p4)

        # Check order
        self.assertEqual(m._parameters.keys(), ["p1", "xyz"])
        self.assertEqual(m._parameters.values(), [p2, p4])

        return

    def testRemoveParameter(self):
        """Test removeParameter method."""

        m = self.m

        p1 = Parameter("p1", 1)
        p2 = Parameter("p1", 2)

        m._addParameter(p1)

        # Check for bad remove
        self.assertRaises(ValueError, m._removeParameter, p2)

        # Remove p1
        m._removeParameter(p1)
        self.assertTrue(p1.name not in m._eqfactory.builders)

        # Try to remove it again
        self.assertRaises(ValueError, m._removeParameter, p1)

        # Try to remove a RecipeContainer
        c = RecipeContainer("test")
        self.assertRaises(ValueError, m._removeParameter, c)
        return

    def testConstrain(self):
        """Test the constrain method."""

        p1 = self.m._newParameter("p1", 1)
        p2 = self.m._newParameter("p2", 2)
        p3 = Parameter("p3", 3)

        self.assertFalse(p1.constrained)
        self.assertEqual(0, len(self.m._constraints))
        self.m.constrain(p1, "2*p2")


        self.assertTrue(p1.constrained)
        self.assertTrue(p1 in self.m._constraints)
        self.assertEqual(1, len(self.m._constraints))
        self.assertTrue(self.m.isConstrained(p1))

        p2.setValue(10)
        self.m._constraints[p1].update()
        self.assertEqual(20, p1.getValue())

        # Check errors on unregistered parameters
        self.assertRaises(ValueError, self.m.constrain, p1, "2*p3")
        self.assertRaises(ValueError, self.m.constrain, p1, "2*p2", {"p2":p3})

        # Remove the constraint
        self.m.unconstrain(p1)
        self.assertFalse(p1.constrained)
        self.assertEqual(0, len(self.m._constraints))
        self.assertFalse(self.m.isConstrained(p1))

        # Try an straight constraint
        self.m.constrain(p1, p2)
        p2.setValue(7)
        self.m._constraints[p1].update()
        self.assertEqual(7, p1.getValue())
        return

    def testRestrain(self):
        """Test the restrain method."""

        p1 = Parameter("p1", 1)
        p2 = Parameter("p2", 2)
        p3 = Parameter("p3", 3)
        self.m._eqfactory.registerArgument("p1", p1)
        self.m._eqfactory.registerArgument("p2", p2)

        self.assertEqual(0, len(self.m._restraints))
        r = self.m.restrain("p1+p2", ub = 10)
        self.assertEqual(1, len(self.m._restraints))
        p2.setValue(10)
        self.assertEqual(1, r.penalty())
        self.m.unrestrain(r)
        self.assertEqual(0, len(self.m._restraints))

        r = self.m.restrain(p1, ub = 10)
        self.assertEqual(1, len(self.m._restraints))
        p1.setValue(11)
        self.assertEqual(1, r.penalty())

        # Check errors on unregistered parameters
        self.assertRaises(ValueError, self.m.restrain, "2*p3")
        self.assertRaises(ValueError, self.m.restrain, "2*p2", ns = {"p2":p3})
        return

    def testGetConstraints(self):
        """Test the _getConstraints method."""
        m2 = RecipeOrganizer("m2")
        self.m._organizers = {}
        self.m._manage(self.m._organizers)
        self.m._addObject(m2, self.m._organizers)

        p1 = Parameter("p1", 1)
        p2 = Parameter("p2", 2)
        p3 = Parameter("p3", 3)
        p4 = Parameter("p4", 4)

        self.m._addParameter(p1)
        self.m._addParameter(p2)

        m2._addParameter(p3)
        m2._addParameter(p4)

        self.m.constrain(p1, "p2")
        m2.constrain(p3, "p4")

        cons = self.m._getConstraints()
        self.assertTrue(p1 in cons)
        self.assertTrue(p3 in cons)
        self.assertEqual(2, len(cons))
        return

    def testGetRestraints(self):
        """Test the _getRestraints method."""
        m2 = RecipeOrganizer("m2")
        self.m._organizers = {}
        self.m._manage(self.m._organizers)
        self.m._addObject(m2, self.m._organizers)

        p1 = Parameter("p1", 1)
        p2 = Parameter("p2", 2)
        p3 = Parameter("p3", 3)
        p4 = Parameter("p4", 4)

        self.m._addParameter(p1)
        self.m._addParameter(p2)

        m2._addParameter(p3)
        m2._addParameter(p4)

        r1 = self.m.restrain("p1 + p2")
        r2 = m2.restrain("2*p3 + p4")

        res = self.m._getRestraints()
        self.assertTrue(r1 in res)
        self.assertTrue(r2 in res)
        self.assertEqual(2, len(res))
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

        self.m.g.center.setValue(5.0)

        self.assertTrue(numpy.array_equal(numpy.exp(-0.5*((x-5.0)/0.1)**2),
            g(x)))

        # Use this in another equation

        eq = self.m.registerStringFunction("g/x - 1", "pdf")
        self.assertTrue(numpy.array_equal(g(x)/x - 1, eq()))

        return

    def testRegisterFunction(self):
        """Test registering various functions."""
        def g1(A, c, w, x):
            return A * numpy.exp(-0.5*((x-c)/w)**2)
        def g2(A):
            return A+1

        eq = self.m.registerFunction(g1, "g")

        for p in eq.args:
            self.assertTrue(p in self.m._parameters.values())

        x = numpy.arange(0.5, 10, 0.5)
        self.m.x.setValue(x)
        self.m.A.setValue(1.0)
        self.m.c.setValue(3.0)
        self.m.w.setValue(0.1)

        self.assertTrue(numpy.array_equal(numpy.exp(-0.5*((x-3.0)/0.1)**2),
            eq()))

        # Use this in another equation
        eq2 = self.m.registerStringFunction("g/x - 1", "pdf")
        self.assertTrue(numpy.array_equal(eq()/x - 1, eq2()))

        # Make sure we can swap out "g".
        self.m.registerFunction(g2, "g")
        self.assertAlmostEqual(2.0, eq())

        # Try a bound method
        class temp(object):
            def eval(self): return 1.23
            def __call__(self): return 4.56

        t = temp()
        eq = self.m.registerFunction(t.eval, "eval")
        self.assertAlmostEqual(1.23, eq())

        # Now the callable
        eq2 = self.m.registerFunction(t, "eval2")
        self.assertAlmostEqual(4.56, eq2())

        return

    def testRegisterStringFunction(self):
        """Test registering string functions in various ways."""

        # Make an equation.
        eq1 = self.m.registerStringFunction("x**2 + 3", "eq1")
        eq1.x.setValue(0)

        for p in eq1.args:
            self.assertTrue(p in self.m._parameters.values())

        # Add a parameter
        self.m._newParameter("y", 3.0)

        # Make sure that x and y are in the organizer
        self.assertEqual(0, self.m.x.getValue())
        self.assertEqual(3.0, self.m.y.getValue())

        # Use eq1 in some equations

        # x**2 (x**2 + 3 - 3)
        eq2 = self.m.registerStringFunction("eq1 - 3", "eq2")
        # y**2
        eq3 = self.m.registerStringFunction("eq1(y) - 3", "eq3")

        # Test these equations.
        self.assertEqual(0, eq2())
        self.assertEqual(9.0, eq3())
        # Note that eq3 injects the value of y into the argument of eq1, which
        # is x. Thus, calling eq3 sets x to 3.
        self.assertEqual(9.0, eq2())

        # One more level of embedding
        # 2*y**2
        eq4 = self.m.registerStringFunction("2*eq3", "eq4")
        self.assertEqual(18.0, eq4())

        return


    def test_releaseOldEquations(self):
        """Verify EquationFactory does not hold temporary equations.
        """
        self.m._newParameter('x', 12)
        self.assertEqual(36, self.m.evaluateEquation('3 * x'))
        self.assertEqual(0, len(self.m._eqfactory.equations))
        return


    def test_show(self):
        """Verify output from the show function.
        """
        def capture_show(*args, **kwargs):
            sys.stdout = cStringIO.StringIO()
            self.m.show(*args, **kwargs)
            rv = sys.stdout.getvalue()
            sys.stdout = sys.__stdout__
            return rv
        self.assertEqual('', capture_show())
        self.m._newParameter('x', 1)
        self.m._newParameter('y', 2)
        out1 = capture_show()
        lines1 = out1.strip().split('\n')
        self.assertEqual(4, len(lines1))
        self.assertTrue('Parameters' in lines1)
        self.assertFalse('Constraints' in lines1)
        self.assertFalse('Restraints' in lines1)
        self.m._newParameter('z', 7)
        self.m.constrain('y', '3 * z')
        out2 = capture_show()
        lines2 = out2.strip().split('\n')
        self.assertEqual(9, len(lines2))
        self.assertTrue('Parameters' in lines2)
        self.assertTrue('Constraints' in lines2)
        self.assertFalse('Restraints' in lines2)
        self.m.restrain('z', lb=2, ub=3, sig=0.001)
        out3 = capture_show()
        lines3 = out3.strip().split('\n')
        self.assertEqual(13, len(lines3))
        self.assertTrue('Parameters' in lines3)
        self.assertTrue('Constraints' in lines3)
        self.assertTrue('Restraints' in lines3)
        out4 = capture_show(pattern='x')
        lines4 = out4.strip().split('\n')
        self.assertEqual(9, len(lines4))
        out5 = capture_show(pattern='^')
        self.assertEqual(out3, out5)
        # check output with another level of hierarchy
        self.m._addObject(RecipeOrganizer("foo"), self.m._containers)
        self.m.foo._newParameter("bar", 13)
        out6 = capture_show()
        self.assertTrue("foo.bar" in out6)
        # filter out foo.bar
        out7 = capture_show('^(?!foo).')
        self.assertFalse("foo.bar" in out7)
        self.assertEqual(out3, out7)
        return

# ----------------------------------------------------------------------------

if __name__ == "__main__":
    unittest.main()
