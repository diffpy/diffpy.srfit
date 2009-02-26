#!/usr/bin/env python

def callfunction(f, *args):
    """Call a numpy function."""
    import _ufunc
    return _ufunc.callUFunc(f, args)

if __name__ == "__main__":

    import numpy
    import _ufunc

    a = numpy.arange(10)
    b = 0.5

    c = callfunction(numpy.add, a, b)
    print c

    callfunction(numpy.multiply, a, b, c)
    print c

    d = callfunction(numpy.sin, c)
    print d

