#!/usr/bin/env python
"""Tests for refinableobj module."""

import diffpy.srfit.equation.builder as builder
import diffpy.srfit.equation.literals as literals

import unittest

import numpy

from utils import _makeArgs

class TestRegistration(unittest.TestCase):

    def testRegisterArg(self):

        factory = builder.EquationFactory()

        v1 = _makeArgs(1)[0]

        factory.registerArgument("v1", v1)

        eq = factory.makeEquation("v1")

        self.assertTrue(v1 is eq.args[0])
        self.assertEquals(1, len(eq.args))

        # Try to parse an equation with buildargs turned off
        self.assertRaises(ValueError, factory.makeEquation, "v1 + v2", False)

        # Make sure we can still use constants
        eq = factory.makeEquation("v1 + 2", False)
        self.assertTrue(v1 is eq.args[0])
        self.assertEquals(1, len(eq.args))
        return



class TestEquationParser(unittest.TestCase):

    def testParseEquation(self):

        from numpy import exp, sin, divide, sqrt, array_equal, e

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


        # Vector equation
        eq = factory.makeEquation("sqrt(e**(-0.5*(x/sigma)**2))")
        x = numpy.arange(0, 1, 0.05)
        sigma = 0.1
        eq.x.setValue(x)
        eq.sigma.setValue(sigma)
        f = lambda x, sigma : sqrt(e**(-0.5*(x/sigma)**2))
        self.assertTrue(array_equal(eq(), f(x,sigma)))


        # Equation with constants
        factory.registerConstant("x", x)
        eq = factory.makeEquation("sqrt(e**(-0.5*(x/sigma)**2))")
        self.assertTrue("sigma" in eq.argdict)
        self.assertTrue("x" not in eq.argdict)
        self.assertTrue(array_equal(eq(sigma=sigma), f(x,sigma)))


        # Equation with user-defined functions
        factory.registerFunction("myfunc", eq, 1)
        eq2 = factory.makeEquation("c*myfunc(sigma)")
        self.assertTrue(array_equal(eq2(c=2, sigma=sigma), 2*f(x,sigma)))
        self.assertTrue("sigma" in eq2.argdict)
        self.assertTrue("c" in eq2.argdict)

        # Equation with partition
        p1 = literals.Partition("p1")
        v1, v2 = _makeArgs(2)
        p1.addArgument(v1)
        p1.addArgument(v2)
        factory.registerPartition("p1", p1)
        eq = factory.makeEquation("A*p1 + B")
        eq.A.setValue(2)
        eq.B.setValue(4)
        self.assertEquals( (2*1+4)+(2*2+4), eq() )


        # Equation with Generator
        g1 = literals.Generator("g1")
        g1.literal = p1
        factory.registerGenerator("g1", g1)
        eq = factory.makeEquation("A*g1 + B")
        eq.A.setValue(2)
        eq.B.setValue(4)
        self.assertEquals( (2*1+4)+(2*2+4), eq() )


        # Partition equation with tags
        p1 = literals.Partition("p1")
        v1, v2 = _makeArgs(2)
        p1.addArgument(v1, "tag", "tag1")
        p1.addArgument(v2, "tag", "tag2")
        factory.registerPartition("p1", p1)
        eq = factory.makeEquation("add(A*p1, B, 'tag1')")
        eq.A.setValue(2)
        eq.B.setValue(4)
        # Addition should only apply to tag1, multiplication applies to both
        self.assertEquals( (2*1+4)+(2*2+0), eq() )

        eq = factory.makeEquation("add( multiply(A, p1, 'tag1'), B, 'tag2')")
        eq.A.setValue(2)
        eq.B.setValue(4)
        # Addition should only apply to tag2, multiplication applies only to
        # tag1
        self.assertEquals( (2*1+0)+(2+4), eq() )

        eq = factory.makeEquation("multiply(A, p1, 'tag1', combine=True) + c")
        eq.A.setValue(2)
        eq.c.setValue(1)
        # Multiplication should only apply to tag2, 'c' should only be added
        # once, since  the partition is combined by the multiplication
        self.assertEquals( ((2*1)+(2))+1, eq() )

        return

    def testBuildEquation(self):

        from numpy import array_equal

        # simple equation
        sin = builder.sin
        a = builder.ArgumentBuilder(name="a", value = 1)
        A = builder.ArgumentBuilder(name="A", value = 2)
        x = numpy.arange(0, numpy.pi, 0.1)

        beq = A*sin(a*x)
        eq = beq.getEquation()

        self.assertTrue("a" in eq.argdict)
        self.assertTrue("A" in eq.argdict)
        self.assertTrue(array_equal(eq(), 2*numpy.sin(x)))


        # custom function
        def _f(a, b):
            return (a-b)*1.0/(a+b)

        f = builder.wrapFunction("f", _f, 2, 1)
        sin = builder.sin
        a = builder.ArgumentBuilder(name="a", value = 2)
        b = builder.ArgumentBuilder(name="b", value = 1)

        beq = sin(f(a,b))
        eq = beq.getEquation()
        self.assertEqual(eq(), numpy.sin(_f(2, 1)))

        # complex function
        sqrt = builder.sqrt
        e = numpy.e
        _x = numpy.arange(0, 1, 0.05)
        x = builder.ArgumentBuilder(name="x", value = _x, const = True)
        sigma = builder.ArgumentBuilder(name="sigma", value = 0.1)
        beq = sqrt(e**(-0.5*(x/sigma)**2))
        eq = beq.getEquation()
        f = lambda x, sigma : sqrt(e**(-0.5*(x/sigma)**2))
        self.assertTrue(array_equal(eq(), numpy.sqrt(e**(-0.5*(_x/0.1)**2))))


        # Equation with partition
        _p1 = literals.Partition("p1")
        v1, v2 = _makeArgs(2)
        _p1.addArgument(v1)
        _p1.addArgument(v2)
        A = builder.ArgumentBuilder(name="A", value = 2)
        B = builder.ArgumentBuilder(name="B", value = 4)
        p1 = builder.wrapPartition("p1", _p1)
        beq = A*p1 + B
        eq = beq.getEquation()
        self.assertEquals( (2*1+4)+(2*2+4), eq() )
        

        # Equation with Generator
        _g1 = literals.Generator("g1")
        _g1.literal = _p1
        g1 = builder.wrapGenerator("g1", _g1)
        geq = A*g1 + B
        eq = geq.getEquation()
        self.assertEquals( (2*1+4)+(2*2+4), eq() )


        # Equation with conditional operator
        add =  builder.add
        _p1 = literals.Partition("p1")
        v1, v2 = _makeArgs(2)
        _p1.addArgument(v1, "tag1")
        _p1.addArgument(v2, "tag2")
        p1 = builder.wrapPartition("p1", _p1)
        beq = add(A*p1, B, "tag1")
        eq = beq.getEquation()
        self.assertEquals( (2*1+4)+(2*2+0), eq() )
        
        return


if __name__ == "__main__":

    unittest.main()

