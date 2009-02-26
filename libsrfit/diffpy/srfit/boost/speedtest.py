#!/usr/bin/env python
"""Tests for refinableobj module."""

from __init__ import *

import numpy

def _makeArgs(num):
    args = []
    for i in xrange(num):
        args.append(Argument())
        args[-1].name = "v%i"%(i+1,)
    return args

x = numpy.arange(0, 20, 0.05)


class AdditionOperator(UFuncOperator):
    """Addition operator."""

    def __init__(self):
        """Initialization."""
        UFuncOperator.__init__(self)
        self.name = "add"
        self.setUFunc(numpy.add, "+")
        return

class SubtractionOperator(UFuncOperator):
    """Subtraction operator."""

    def __init__(self):
        """Initialization."""
        UFuncOperator.__init__(self)
        self.name = "subtract"
        self.setUFunc(numpy.subtract, "-")
        return

class MultiplicationOperator(UFuncOperator):
    """Multiplication operator."""

    def __init__(self):
        """Initialization."""
        UFuncOperator.__init__(self)
        self.name = "multiply"
        self.setUFunc(numpy.multiply, "*")
        return

class DivisionOperator(UFuncOperator):
    """Division operator."""

    def __init__(self):
        """Initialization."""
        UFuncOperator.__init__(self)
        self.name = "divide"
        self.setUFunc(numpy.divide, "/")
        return

class ExponentiationOperator(UFuncOperator):
    """Exponentiation operator."""

    def __init__(self):
        """Initialization."""
        UFuncOperator.__init__(self)
        self.name = "power"
        self.setUFunc(numpy.power, "**")
        return


class RemainderOperator(UFuncOperator):
    """Remainder operator."""

    def __init__(self):
        """Initialization."""
        UFuncOperator.__init__(self)
        self.name = "mod"
        self.setUFunc(numpy.mod, "%")
        return

def testOperator():
    o = MultiplicationOperator()
    print o.callFunction((3,4))


def makeLazyEquation():
    """Make a lazy equation and see how fast it is."""

    # Make some variables
    v1, v2, v3, v4, v5, v6, v7 = _makeArgs(7)

    # Make some operations
    mult = MultiplicationOperator()
    plus = AdditionOperator()
    minus = SubtractionOperator()
    pow = ExponentiationOperator()
    exp = ExponentiationOperator()
    mult2 = MultiplicationOperator()
    
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

    v1.setValue(0)
    v2.setValue(x)
    v3.setValue(50*x)
    v4.setValue(0)
    v5.setValue(2.11)
    v6.setValue(numpy.e)
    v7.setValue(0)

    evaluator = Evaluator()

    v1.identify(evaluator)
    v2.identify(evaluator)

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
        return ((a+x)*(50*x-b))**2.11 * numpy.exp(c)

    return _f

def timeFunction(f, *args):
    """Time a function in ms."""
    import time
    t1 = time.time()
    f(*args)
    t2 = time.time()
    return (t2-t1)*1000

if __name__ == "__main__":

    f1 = makeLazyEquation()
    f2 = makeEquation()

    a, b, c = 3.1, 8.19973123410, 2.1

    t1 = timeFunction(f1, a, b, c)
    t2 = timeFunction(f2, a, b, c)
    print "Initial"
    print "lazy", t1
    print "regular", t2

    import time
    import _purespeed
    t1 = time.time()
    _purespeed.speedy()
    t2 = time.time()
    print "Pure Speed", (t2-t1)*1000


    # Change one of the variables
    b = 10.1
    t1 = timeFunction(f1, a, b, c)
    t2 = timeFunction(f2, a, b, c)
    print "Single change"
    print "lazy", t1
    print "regular", t2

    # Change two
    a = 2.0
    b = 2.0
    t1 = timeFunction(f1, a, b, c)
    t2 = timeFunction(f2, a, b, c)
    print "Two change"
    print "lazy", t1
    print "regular", t2

    # Change three
    a = 1.0
    b = 1.0
    c = 1.0
    t1 = timeFunction(f1, a, b, c)
    t2 = timeFunction(f2, a, b, c)
    print "Three change"
    print "lazy", t1
    print "regular", t2
