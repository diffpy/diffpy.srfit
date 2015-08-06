#!/usr/bin/env python
########################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2008 The Trustees of Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Chris Farrow
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
########################################################################

"""Complex modeling framework for structure refinement and solution.

SrFit is a tool for coherently combining known information about a material to
derive other properties, in particular material structure. SrFit allows the
customization and creation of structure representations, profile calculators,
constraints, restraints and file input parsers. The customized pieces can be
glued together within SrFit to optimize a structure, or other physically
relevant information from one or more experimental profiles. Other known
information about the system of interest can be included with arbitrarily
complex constraints and restraints. In this way, the end user creates a
customized fitting application that suits the problem to the available
information.

The subpackages herein define various pieces of the SrFit framework. Developers
are encouraged to work through the examples described in the documentation to
learn how to use and customize the various parts of SrFit.

"""
from __future__ import print_function
import six
__all__ = ["__version__"]

# package version
from diffpy.srfit.version import __version__

# Support for method and ufunc pickling.
def _pickle_method(method):
    name = method.__func__.__name__
    obj = method.__self__
    cls = method.__self__.__class__
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
try:
    from copy_reg import pickle
except ImportError:
    # python 3 path
    from copyreg import pickle
pickle(types.MethodType, _pickle_method, _unpickle_method)
pickle(numpy.ufunc, _pickle_ufunc, _unpickle_ufunc)


# End of file
