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

from __future__ import print_function

import random
import numpy

import diffpy.srfit.equation.visitors as visitors
import diffpy.srfit.equation.literals as literals

from diffpy.srfit.tests.utils import _makeArgs

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
    f1 = makeLazyEquation()
    f2 = makeEquation1()

    args = [3.1, 8.19973123410, 2.1, numpy.e, numpy.pi]

    total1 = 0
    total2 = 0
    for i in xrange(len(args)):
        args[i] = 10*random.random()
        print("Changing argument %i" % (i + 1))
        t1 = timeFunction(f1, *args)
        t2 = timeFunction(f2, *args)
        total1 += t1
        total2 += t2
        print("lazy", t1)
        print("regular", t2)

    print("Totals:")
    print("lazy", total1)
    print("regular", total2)
    print("Ratio (lazy/regular)", total1/total2)

def speedTest2(mutate = 2):

    from diffpy.srfit.equation.builder import EquationFactory
    factory = EquationFactory()

    x = numpy.arange(0, 20, 0.05)
    qsig = 0.01
    sigma = 0.003

    eqstr = """\
    A0*exp(-(x*qsig)**2)*(exp(-((x-1.0)/sigma1)**2)+exp(-((x-2.0)/sigma2)**2))\
    + polyval(list(b1, b2, b3, b4, b5, b6, b7, b8), x)\
    """
    factory.registerConstant("x", x)
    eq = factory.makeEquation(eqstr)
    eq.qsig.setValue(qsig)
    eq.sigma1.setValue(sigma)
    eq.sigma2.setValue(sigma)
    eq.A0.setValue(1.0)
    eq.b1.setValue(0)
    eq.b2.setValue(1)
    eq.b3.setValue(2.0)
    eq.b4.setValue(2.0)
    eq.b5.setValue(2.0)
    eq.b6.setValue(2.0)
    eq.b7.setValue(2.0)
    eq.b8.setValue(2.0)

    from numpy import exp
    from numpy import polyval
    def f(A0, qsig, sigma1, sigma2, b1, b2, b3, b4, b5, b6, b7, b8):
        return A0*exp(-(x*qsig)**2)*(exp(-((x-1.0)/sigma1)**2)+exp(-((x-2.0)/sigma2)**2)) + polyval([b8, b7, b6, b5,b4,b3,b2,b1],x)

    tnpy = 0
    teq = 0
    # Randomly change variables
    numargs = len(eq.args)
    choices = range(numargs)
    args = [0.0]*(len(eq.args))

    # The call-loop
    random.seed()
    numcalls = 1000
    for _i in xrange(numcalls):
        # Mutate values
        n = mutate
        if n == 0:
            n = random.choice(choices)
        c = choices[:]
        for _j in xrange(n):
            idx = random.choice(c)
            c.remove(idx)
            args[idx] = random.random()

        # Time the different functions with these arguments
        tnpy += timeFunction(f, *args)
        teq += timeFunction(eq, *args)

    print("Average call time (%i calls, %i mutations/call):" %
          (numcalls, mutate))
    print("numpy: ", tnpy/numcalls)
    print("equation: ", teq/numcalls)
    print("ratio: ", teq/tnpy)

    return

def speedTest3(mutate = 2):
    """Test wrt sympy.

    Results - sympy is 10 to 24 times faster without using arrays (ouch!).
            - diffpy.srfit.equation is slightly slower when using arrays, but
              not considerably worse than versus numpy alone.

    """

    from diffpy.srfit.equation.builder import EquationFactory
    factory = EquationFactory()

    x = numpy.arange(0, 20, 0.05)
    qsig = 0.01
    sigma = 0.003

    eqstr = """\
    A0*exp(-(x*qsig)**2)*(exp(-((x-1.0)/sigma1)**2)+exp(-((x-2.0)/sigma2)**2))\
    + polyval(list(b1, b2, b3, b4, b5, b6, b7, b8), x)\
    """
    factory.registerConstant("x", x)
    eq = factory.makeEquation(eqstr)
    eq.qsig.setValue(qsig)
    eq.sigma1.setValue(sigma)
    eq.sigma2.setValue(sigma)
    eq.A0.setValue(1.0)
    eq.b1.setValue(0)
    eq.b2.setValue(1)
    eq.b3.setValue(2.0)
    eq.b4.setValue(2.0)
    eq.b5.setValue(2.0)
    eq.b6.setValue(2.0)
    eq.b7.setValue(2.0)
    eq.b8.setValue(2.0)

    from sympy import var, exp, lambdify
    from numpy import polyval
    A0, qsig, sigma1, sigma2, b1, b2, b3, b4, b5, b6, b7, b8, xx = vars = var("A0 qsig sigma1 sigma2 b1 b2 b3 b4 b5 b6 b7 b8 xx")
    f = lambdify(vars, A0*exp(-(xx*qsig)**2)*(exp(-((xx-1.0)/sigma1)**2)+exp(-((xx-2.0)/sigma2)**2)) + polyval([b1, b2, b3, b4, b5, b6, b7, b8], xx), "numpy")

    tnpy = 0
    teq = 0
    # Randomly change variables
    numargs = len(eq.args)
    choices = range(numargs)
    args = [1.0]*(len(eq.args))
    args.append(x)

    # The call-loop
    random.seed()
    numcalls = 1000
    for _i in xrange(numcalls):
        # Mutate values
        n = mutate
        if n == 0:
            n = random.choice(choices)
        c = choices[:]
        for _j in xrange(n):
            idx = random.choice(c)
            c.remove(idx)
            args[idx] = random.random()

        # Time the different functions with these arguments
        teq += timeFunction(eq, *(args[:-1]))
        tnpy += timeFunction(f, *args)

    print("Average call time (%i calls, %i mutations/call):" %
          (numcalls, mutate))
    print("sympy: ", tnpy/numcalls)
    print("equation: ", teq/numcalls)
    print("ratio: ", teq/tnpy)

    return

def speedTest4(mutate = 2):
    """Test wrt sympy.

    Results - sympy is 10 to 24 times faster without using arrays (ouch!).
            - diffpy.srfit.equation is slightly slower when using arrays, but
              not considerably worse than versus numpy alone.

    """

    from diffpy.srfit.equation.builder import EquationFactory
    factory = EquationFactory()

    x = numpy.arange(0, 20, 0.05)

    eqstr = """\
    b1 + b2*x + b3*x**2 + b4*x**3 + b5*x**4 + b6*x**5 + b7*x**6 + b8*x**7\
    """
    factory.registerConstant("x", x)
    eq = factory.makeEquation(eqstr)

    from sympy import var, lambdify
    from numpy import polyval
    b1, b2, b3, b4, b5, b6, b7, b8, xx = vars = var("b1 b2 b3 b4 b5 b6 b7 b8 xx")
    f = lambdify(vars, polyval([b1, b2, b3, b4, b5, b6, b7, b8], xx), "numpy")

    tnpy = 0
    teq = 0
    # Randomly change variables
    numargs = len(eq.args)
    choices = range(numargs)
    args = [1.0]*(len(eq.args))
    args.append(x)

    # The call-loop
    random.seed()
    numcalls = 1000
    for _i in xrange(numcalls):
        # Mutate values
        n = mutate
        if n == 0:
            n = random.choice(choices)
        c = choices[:]
        for _j in xrange(n):
            idx = random.choice(c)
            c.remove(idx)
            args[idx] = random.random()

        # Time the different functions with these arguments
        teq += timeFunction(eq, *(args[:-1]))
        tnpy += timeFunction(f, *args)

    print("Average call time (%i calls, %i mutations/call):" %
          (numcalls, mutate))
    print("sympy: ", tnpy/numcalls)
    print("equation: ", teq/numcalls)
    print("ratio: ", teq/tnpy)

    return

def weightedTest(mutate = 2):
    """Show the benefits of a properly balanced equation tree."""

    from diffpy.srfit.equation.builder import EquationFactory
    factory = EquationFactory()

    x = numpy.arange(0, 10, 0.01)

    eqstr = """\
    b1 + b2*x + b3*x**2 + b4*x**3 + b5*x**4 + b6*x**5 + b7*x**6 + b8*x**7\
    """
    factory.registerConstant("x", x)
    eq = factory.makeEquation(eqstr)

    eq.b1.setValue(0)
    eq.b2.setValue(1)
    eq.b3.setValue(2.0)
    eq.b4.setValue(2.0)
    eq.b5.setValue(2.0)
    eq.b6.setValue(2.0)
    eq.b7.setValue(2.0)
    eq.b8.setValue(2.0)

    #scale = visitors.NodeWeigher()
    #eq.root.identify(scale)
    #print scale.output

    from numpy import polyval
    def f(b1, b2, b3, b4, b5, b6, b7, b8):
        return polyval([b8, b7, b6, b5,b4,b3,b2,b1],x)

    tnpy = 0
    teq = 0
    # Randomly change variables
    numargs = len(eq.args)
    choices = range(numargs)
    args = [0.1]*numargs

    # The call-loop
    random.seed()
    numcalls = 1000
    for _i in xrange(numcalls):
        # Mutate values
        n = mutate
        if n == 0:
            n = random.choice(choices)
        c = choices[:]
        for _j in xrange(n):
            idx = random.choice(c)
            c.remove(idx)
            args[idx] = random.random()

        #print args

        # Time the different functions with these arguments
        teq += timeFunction(eq, *args)
        tnpy += timeFunction(f, *args)

    print("Average call time (%i calls, %i mutations/call):" %
          (numcalls, mutate))
    print("numpy: ", tnpy/numcalls)
    print("equation: ", teq/numcalls)
    print("ratio: ", teq/tnpy)

    return

def profileTest():

    from diffpy.srfit.builder import EquationFactory
    factory = EquationFactory()

    x = numpy.arange(0, 10, 0.001)

    eqstr = """\
    b1 + b2*x + b3*x**2 + b4*x**3 + b5*x**4 + b6*x**5 + b7*x**6 + b8*x**7\
    """
    factory.registerConstant("x", x)
    eq = factory.makeEquation(eqstr)

    eq.b1.setValue(0)
    eq.b2.setValue(1)
    eq.b3.setValue(2.0)
    eq.b4.setValue(2.0)
    eq.b5.setValue(2.0)
    eq.b6.setValue(2.0)
    eq.b7.setValue(2.0)
    eq.b8.setValue(2.0)

    mutate = 8
    numargs = len(eq.args)
    choices = range(numargs)
    args = [0.1]*numargs

    # The call-loop
    random.seed()
    numcalls = 1000
    for _i in xrange(numcalls):
        # Mutate values
        n = mutate
        if n == 0:
            n = random.choice(choices)
        c = choices[:]
        for _j in xrange(n):
            idx = random.choice(c)
            c.remove(idx)
            args[idx] = random.random()

        eq(*args)

    return


if __name__ == "__main__":
    for i in range(1, 13):
        speedTest2(i)
    """
    for i in range(1, 9):
        weightedTest(i)
    """
    """
    from diffpy.srfit.equation.builder import EquationFactory
    import random
    import cProfile
    cProfile.run('profileTest()', 'prof')
    import pstats
    p = pstats.Stats('prof')
    p.strip_dirs()
    p.sort_stats('time')
    p.print_stats(10)
    profileTest()
    """
