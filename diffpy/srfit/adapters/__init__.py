########################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2011 Trustees of the Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Chris Farrow
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
########################################################################

"""Base adapters for diffpy.srfit.

"""

__all__ = ["adapt"]

from diffpy.srfit.adapters.adaptersmod import adapt

# Support for method and ufunc pickling.
def _pickle_method(method):
    name = method.im_func.__name__
    obj = method.im_self
    cls = method.im_class
    return _unpickle_method, (name, obj, cls)

def _unpickle_method(name, obj, cls):
    for cls in cls.mro():
        try:
            func = cls.__dict__[name]
        except KeyError:
            pass
        else:
            break
    return func.__get__(obj, cls)

# This requires the ufunc is defined in the numpy module. This is good for most
# of what we're doing.
def _pickle_ufunc(func):
    name = func.__name__
    return _unpickle_ufunc, (name,)

import numpy
def _unpickle_ufunc(name):
    func = getattr(numpy, name)
    return func

import types
import copy_reg
copy_reg.pickle(types.MethodType, _pickle_method, _unpickle_method)
copy_reg.pickle(numpy.ufunc, _pickle_ufunc, _unpickle_ufunc)

__id__ = "$Id$"
# End of file
