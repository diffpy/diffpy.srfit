#!/usr/bin/env python

__all__ = []

import numpy
import sys

from diffpy.srfit.adapters.adaptersmod import UnboundOperator

def __adaptUFuncs():
    """Adapt all ufuncs from numpy."""
    module = sys.modules[__name__]
    for name in dir(numpy):
        op = getattr(numpy, name)
        if isinstance(op, numpy.ufunc):
            module.__all__.append(name)
            setattr(module, name, UnboundOperator(name, op))
    return
__adaptUFuncs()

array = UnboundOperator("array", numpy.array)
__all__.append("array")

dot = UnboundOperator("dot", numpy.dot)
__all__.append("dot")

sum = UnboundOperator("sum", numpy.sum)
__all__.append("sum")

concatenate = UnboundOperator("concatenate", numpy.concatenate)
__all__.append("concatenate")

interp = UnboundOperator("interp", numpy.interp)
__all__.append("interp")
