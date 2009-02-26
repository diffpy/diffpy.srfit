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

def makeEquation():
    """Make the same equation as the lazy one."""

    y = 50*x

    def _f(a, b, c, d, e):
        return ((a+x)*(y-b))**c * d**e

    return _f

def timeFunction(f, *args):
    """Time a function in ms."""
    import time
    t1 = time.time()
    f(*args)
    t2 = time.time()
    return (t2-t1)*1000

if __name__ == "__main__":

    import random
    f1 = makeLazyEquation()
    f2 = makeEquation()

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

