#!/usr/bin/env python

from __init__ import *

if __name__ == "__main__":

    import numpy
    o = UFuncOperator()
    o.setUFunc(numpy.sin, "sin")
    print o.callFunction((3,))

    o.setUFunc(numpy.add, "+")
    print o.callFunction((8,9))

    o.setUFunc(numpy.multiply, "*")
    print o.callFunction((8,numpy.arange(0, 2, 0.1)))

    o.setUFunc(numpy.cos, "cos")
    print o.callFunction((3,))
