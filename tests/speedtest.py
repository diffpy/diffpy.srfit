#!/usr/bin/env python
"""Tests for refinableobj module."""

import diffpy.srfit.equation.visitors as visitors
import diffpy.srfit.equation.literals as literals

import numpy

def _makeArgs(num):
    args = []
    for i in xrange(num):
        args.append(literals.Argument())
        args[-1].name = "v%i"%(i+1)
    return args

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

    def _f(a, b, c):
        v1.setValue(a)
        v4.setValue(b)
        v7.setValue(c)
        mult2.identify(evaluator)
        evaluator.clicker.click()
        return evaluator.value

    return _f

def makeEquation():
    """Make the same equation as the lazy one."""

    y = 50*x

    def _f(a, b, c):
        return ((a+x)*(y-b))**2.11 * numpy.exp(c)

    return _f

def timeFunction(f, *args):
    """Time a function in ms."""
    import time
    t1 = time.time()
    f(*args)
    t2 = time.time()
    return (t2-t1)*1000

if __name__ == "__main__":

    print "Initial"
    f1 = makeLazyEquation()
    f2 = makeEquation()

    a, b, c = 3.1, 8.19973123410, 2.1

    t1 = timeFunction(f1, a, b, c)
    t2 = timeFunction(f2, a, b, c)
    print "lazy", t1
    print "regular", t2

    # Change one of the variables
    print "Single change"
    c = 0.1
    t1 = timeFunction(f1, a, b, c)
    t2 = timeFunction(f2, a, b, c)
    print "lazy", t1
    print "regular", t2

    # Change two
    print "Two change"
    a = 2.0
    b = 2.0
    t1 = timeFunction(f1, a, b, c)
    t2 = timeFunction(f2, a, b, c)
    print "lazy", t1
    print "regular", t2

    # Change three
    print "Three change"
    a = 1.0
    b = 1.0
    c = 1.0
    t1 = timeFunction(f1, a, b, c)
    t2 = timeFunction(f2, a, b, c)
    print "lazy", t1
    print "regular", t2

