#!/usr/bin/env python
"""Tests for refinableobj module."""

import diffpy.srfit.equation.visitors as visitors
import diffpy.srfit.equation.literals as literals

import numpy

from utils import _makeArgs

x = numpy.arange(0, 20, 0.05)

def makeLazyEquation():
    """Make a lazy equation and see how fast it is."""

    # Make some variables
    v1, v2, v3, v4, v5, v6, v7 = _makeArgs(7)

    # Make some operations
    mult = literals.MultiplicationOperator()
    plus = literals.AdditionOperator()
    minus = literals.SubtractionOperator()
    pow = literals.ExponentiationOperator()
    exp = literals.ExponentiationOperator()
    mult2 = literals.MultiplicationOperator()
    
    # Create the equation ((v1+v2)*(v3-v4))**v5 * v6^v7
    plus.addLiteral(v1)
    plus.addLiteral(v2)
    minus.addLiteral(v3)
    minus.addLiteral(v4)
    mult.addLiteral(plus)
    mult.addLiteral(minus)
    pow.addLiteral(mult)
    pow.addLiteral(v5)
    exp.addLiteral(v6)
    exp.addLiteral(v7)
    mult2.addLiteral(pow)
    mult2.addLiteral(exp)

    v2.setValue(x)
    v3.setValue(50*x)
    v5.setValue(2.11)
    v6.setValue(numpy.e)

    evaluator = visitors.Evaluator()

    def _f(a, b, c, d, e):
        v1.setValue(a)
        v4.setValue(b)
        v5.setValue(c)
        v6.setValue(d)
        v7.setValue(e)
        mult2.identify(evaluator)
        evaluator.clicker.click()
        return evaluator.value

    return _f

def makeEquation1():
    """Make the same equation as the lazy one."""

    y = 50*x

    def _f(a, b, c, d, e):
        return ((a+x)*(y-b))**c * d**e

    return _f

def timeFunction(f, *args, **kw):
    """Time a function in ms."""
    import time
    t1 = time.time()
    f(*args, **kw)
    t2 = time.time()
    return (t2-t1)*1000

def speedTest1():
    import random
    f1 = makeLazyEquation()
    f2 = makeEquation1()

    args = [3.1, 8.19973123410, 2.1, numpy.e, numpy.pi]

    total1 = 0
    total2 = 0
    for i in xrange(len(args)):
        args[i] = 10*random.random()
        print "Changing argument %i"%(i+1)
        t1 = timeFunction(f1, *args)
        t2 = timeFunction(f2, *args)
        total1 += t1
        total2 += t2
        print "lazy", t1
        print "regular", t2

    print "Totals:"
    print "lazy", total1
    print "regular", total2
    print "Ratio (lazy/regular)", total1/total2

def speedTest2():

    from diffpy.srfit.equation.builder import makeEquation

    x = numpy.arange(0, 20, 0.05)
    qsig = 0.01
    sigma = 0.003


    eqstr = "exp((x*qsig)**2)*(exp(((x-1)/sigma1)**2)+exp(((x-2)/sigma2)**2))"
    eq = makeEquation(eqstr)
    #print eq.args.keys()
    eq.x.setValue(x)
    eq.qsig.setValue(qsig)
    eq.sigma1.setValue(sigma)
    eq.sigma2.setValue(sigma)

    t1 = timeFunction(eq)
    eq.qsig.setValue(0.02)
    t2 = timeFunction(eq)
    eq.sigma2.setValue(0.005)
    t3 = timeFunction(eq)
    eq.sigma1.setValue(0.005)
    t4 = timeFunction(eq)
    print "Mine", t1, t2, t3, t4, t1+t2+t3+t4

    from numpy import exp
    def f(qsig, sigma1, sigma2):
        return exp((x*qsig)**2)*(exp(((x-1)/sigma1)**2)+exp(((x-2)/sigma2)**2))

    t1 = timeFunction(f, qsig, sigma, sigma)
    t2 = timeFunction(f, 0.02, sigma, sigma)
    t3 = timeFunction(f, 0.02, 0.005, sigma)
    t4 = timeFunction(f, 0.02, 0.005, 0.005)
    print "Numpy", t1, t2, t3, t4, t1+t2+t3+t4
    return

if __name__ == "__main__":
    #speedTest1()
    speedTest2()
