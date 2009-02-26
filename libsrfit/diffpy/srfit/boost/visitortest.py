#!/usr/bin/env python

from __init__ import *

def testCallout(CalloutClass, ArgClass):

    o1 = UFuncOperator()
    o1.name = "operator1"

    o2 = UFuncOperator()
    o2.name = "operator2"

    arg1 = ArgClass()
    arg1.name = "arg1"

    arg2 = ArgClass()
    arg2.name = "arg2"

    o1.addLiteral(arg1)
    o1.addLiteral(o2)
    o2.addLiteral(arg2)

    c = CalloutClass()
    o1.identify(c)

    print o1.args
    return


class PyCallout(Visitor):

    def visitArgument(self, a):
        print a.name

    def visitOperator(self, o):
        print o.name
        for l in o.args:
            l.identify(self)


class PyArgument(Argument):
    pass



if __name__ == "__main__":

    print "c++ Visitor, c++ Argument"
    testCallout(Callout, Argument)
    print "python Visitor, c++ Argument"
    testCallout(PyCallout, Argument)
    print "c++ Visitor, python Argument"
    testCallout(Callout, PyArgument)
    print "python Visitor, python Argument"
    testCallout(PyCallout, PyArgument)
