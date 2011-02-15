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
            symname = name + "_"
            module.__all__.append(symname)
            setattr(module, symname, UnboundOperator(name, op))
    return
__adaptUFuncs()

array_ = UnboundOperator("array", numpy.array)
__all__.append("array_")

dot_ = UnboundOperator("dot", numpy.dot)
__all__.append("dot_")

sum_ = UnboundOperator("sum", numpy.sum)
__all__.append("sum_")

concatenate_ = UnboundOperator("concatenate", numpy.concatenate)
__all__.append("concatenate_")

interp_ = UnboundOperator("interp", numpy.interp)
__all__.append("interp_")
