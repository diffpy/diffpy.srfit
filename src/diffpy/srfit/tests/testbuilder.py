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

import diffpy.srfit.equation.builder as builder
import diffpy.srfit.equation.literals as literals

import unittest

import numpy

from diffpy.srfit.tests.utils import _makeArgs
from diffpy.srfit.tests.utils import noObserversInGlobalBuilders

class TestBuilder(unittest.TestCase):

    def testRegisterArg(self):

        factory = builder.EquationFactory()

        v1 = _makeArgs(1)[0]

        b1 = factory.registerArgument("v1", v1)
        self.assertTrue(factory.builders["v1"] is b1)
        self.assertTrue(b1.literal is v1)

        eq = factory.makeEquation("v1")

        self.assertTrue(v1 is eq.args[0])
        self.assertEqual(1, len(eq.args))

        # Try to parse an equation with buildargs turned off
        self.assertRaises(ValueError, factory.makeEquation, "v1 + v2", False)

        # Make sure we can still use constants
        eq = factory.makeEquation("v1 + 2", False)
        self.assertTrue(v1 is eq.args[0])
        self.assertEqual(1, len(eq.args))
        self.assertTrue(noObserversInGlobalBuilders())
        return


    def testRegisterOperator(self):
        """Try to use an operator without arguments in an equation."""

        factory = builder.EquationFactory()
        v1, v2, v3, v4 = _makeArgs(4)

        op = literals.AdditionOperator()

        op.addLiteral(v1)
        op.addLiteral(v2)

        factory.registerArgument("v3", v3)
        factory.registerArgument("v4", v4)
        factory.registerOperator("op", op)

        # Build an equation where op is treated as a terminal node
        eq = factory.makeEquation("op")
        self.assertAlmostEqual(3, eq())

        eq = factory.makeEquation("v3*op")
        self.assertAlmostEqual(9, eq())

        # Now use the op like a function
        eq = factory.makeEquation("op(v3, v4)")
        self.assertAlmostEqual(7, eq())

        # Make sure we can still access op as itself.
        eq = factory.makeEquation("op")
        self.assertAlmostEqual(3, eq())

        self.assertTrue(noObserversInGlobalBuilders())
        return


    def testSwapping(self):

        def g1(v1, v2, v3, v4):
            return (v1 + v2) * (v3 + v4)
        def g2(v1):
            return 0.5*v1

        factory = builder.EquationFactory()
        v1, v2, v3, v4, v5 = _makeArgs(5)

        factory.registerArgument("v1", v1)
        factory.registerArgument("v2", v2)
        factory.registerArgument("v3", v3)
        factory.registerArgument("v4", v4)
        b = factory.registerFunction("g", g1, ["v1", "v2", "v3", "v4"])

        # Now associate args with the wrapped function
        op = b.literal
        self.assertTrue(op.operation == g1)
        self.assertTrue(v1 in op.args)
        self.assertTrue(v2 in op.args)
        self.assertTrue(v3 in op.args)
        self.assertTrue(v4 in op.args)
        self.assertAlmostEqual(21, op.value)

        eq1 = factory.makeEquation("g")
        self.assertTrue(eq1.root is op)
        self.assertAlmostEqual(21, eq1())

        # Swap out an argument by registering it under a taken name
        b = factory.registerArgument("v4", v5)
        self.assertTrue(factory.builders["v4"] is b)
        self.assertTrue(b.literal is v5)
        self.assertTrue(op._value is None)
        self.assertTrue(op.args == [v1, v2, v3, v5])
        self.assertAlmostEqual(24, eq1())

        # Now swap out the function
        b = factory.registerFunction("g", g2, ["v1"])
        op = b.literal
        self.assertTrue(op.operation == g2)
        self.assertTrue(v1 in op.args)
        self.assertTrue(eq1.root is op)
        self.assertAlmostEqual(0.5, op.value)
        self.assertAlmostEqual(0.5, eq1())

        # Make an equation
        eqeq = factory.makeEquation("v1 + v2")
        # Register this "g"
        b = factory.registerFunction("g", eqeq, eqeq.argdict.keys())
        op = b.literal
        self.assertTrue(v1 in op.args)
        self.assertTrue(v2 in op.args)
        self.assertTrue(eq1.root is op)
        self.assertAlmostEqual(3, op.value)
        self.assertAlmostEqual(3, eq1())

        self.assertTrue(noObserversInGlobalBuilders())
        return

    def testParseEquation(self):

        from numpy import sin, divide, sqrt, array_equal, e

        factory = builder.EquationFactory()

        # Scalar equation
        eq = factory.makeEquation("A*sin(0.5*x)+divide(B,C)")
        A = 1
        x = numpy.pi
        B = 4.0
        C = 2.0
        eq.A.setValue(A)
        eq.x.setValue(x)
        eq.B.setValue(B)
        eq.C.setValue(C)
        f = lambda A, x, B, C: A*sin(0.5*x)+divide(B,C)
        self.assertTrue(array_equal(eq(), f(A,x,B,C)))

        # Make sure that the arguments of eq are listed in the order in which
        # they appear in the equations.
        self.assertEqual(eq.args, [eq.A, eq.x, eq.B, eq.C])

        # Vector equation
        eq = factory.makeEquation("sqrt(e**(-0.5*(x/sigma)**2))")
        x = numpy.arange(0, 1, 0.05)
        sigma = 0.1
        eq.x.setValue(x)
        eq.sigma.setValue(sigma)
        f = lambda x, sigma : sqrt(e**(-0.5*(x/sigma)**2))
        self.assertTrue(numpy.allclose(eq(), f(x,sigma)))

        self.assertEqual(eq.args, [eq.x, eq.sigma])

        # Equation with constants
        factory.registerConstant("x", x)
        eq = factory.makeEquation("sqrt(e**(-0.5*(x/sigma)**2))")
        self.assertTrue("sigma" in eq.argdict)
        self.assertTrue("x" not in eq.argdict)
        self.assertTrue(numpy.allclose(eq(sigma=sigma), f(x,sigma)))

        self.assertEqual(eq.args, [eq.sigma])

        # Equation with user-defined functions
        factory.registerFunction("myfunc", eq, ["sigma"])
        eq2 = factory.makeEquation("c*myfunc(sigma)")
        self.assertTrue(numpy.allclose(eq2(c=2, sigma=sigma), 2*f(x,sigma)))
        self.assertTrue("sigma" in eq2.argdict)
        self.assertTrue("c" in eq2.argdict)
        self.assertEqual(eq2.args, [eq2.c, eq2.sigma])

        self.assertTrue(noObserversInGlobalBuilders())
        return


    def test_parse_constant(self):
        """Verify parsing of constant numeric expressions.
        """
        factory = builder.EquationFactory()
        eq = factory.makeEquation('3.12 + 2')
        self.assertTrue(isinstance(eq, builder.Equation))
        self.assertEqual(set(), factory.equations)
        self.assertEqual(5.12, eq())
        self.assertRaises(ValueError, eq, 3)
        return


    def testBuildEquation(self):

        from numpy import array_equal

        # simple equation
        sin = builder.getBuilder("sin")
        a = builder.ArgumentBuilder(name="a", value = 1)
        A = builder.ArgumentBuilder(name="A", value = 2)
        x = numpy.arange(0, numpy.pi, 0.1)

        beq = A*sin(a*x)
        eq = beq.getEquation()

        self.assertTrue("a" in eq.argdict)
        self.assertTrue("A" in eq.argdict)
        self.assertTrue(array_equal(eq(), 2*numpy.sin(x)))

        self.assertEqual(eq.args, [eq.A, eq.a])

        # Check the number of arguments
        self.assertRaises(ValueError, sin)

        # custom function
        def _f(a, b):
            return (a-b)*1.0/(a+b)

        f = builder.wrapFunction("f", _f, 2, 1)
        a = builder.ArgumentBuilder(name="a", value = 2)
        b = builder.ArgumentBuilder(name="b", value = 1)

        beq = sin(f(a,b))
        eq = beq.getEquation()
        self.assertEqual(eq(), numpy.sin(_f(2, 1)))

        # complex function
        sqrt = builder.getBuilder("sqrt")
        e = numpy.e
        _x = numpy.arange(0, 1, 0.05)
        x = builder.ArgumentBuilder(name="x", value = _x, const = True)
        sigma = builder.ArgumentBuilder(name="sigma", value = 0.1)
        beq = sqrt(e**(-0.5*(x/sigma)**2))
        eq = beq.getEquation()
        f = lambda x, sigma : sqrt(e**(-0.5*(x/sigma)**2))
        self.assertTrue(numpy.allclose(eq(), numpy.sqrt(e**(-0.5*(_x/0.1)**2))))

        # Equation with Equation
        A = builder.ArgumentBuilder(name="A", value = 2)
        B = builder.ArgumentBuilder(name="B", value = 4)
        beq = A + B
        eq = beq.getEquation()
        E = builder.wrapOperator("eq", eq)
        eq2 = (2*E).getEquation()
        # Make sure these evaulate to the same thing
        self.assertEqual(eq.args, [A.literal, B.literal])
        self.assertEqual(2*eq(), eq2())
        # Pass new arguments to the equation
        C = builder.ArgumentBuilder(name="C", value = 5)
        D = builder.ArgumentBuilder(name="D", value = 6)
        eq3 = (E(C, D)+1).getEquation()
        self.assertEqual(12, eq3())
        # Pass old and new arguments to the equation
        # If things work right, A has been given the value of C in the last
        # evaluation (5)
        eq4 = (3*E(A, D)-1).getEquation()
        self.assertEqual(32, eq4())
        # Try to pass the wrong number of arguments
        self.assertRaises(ValueError, E, A)
        self.assertRaises(ValueError, E, A, B, C)

        self.assertTrue(noObserversInGlobalBuilders())
        return


if __name__ == "__main__":
    unittest.main()
