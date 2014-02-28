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

"""Import abc module for python2.6 or before.

This imports the abc module in python > 2.6, or a modified abc module for
earlier python versions. If the modified abc module is used,

isinstance and issubclass are redefined here so that ABCs can be used in python
< 2.6. To assure backwards compatibility, import these methods from this
module rather than using the built-in isinstance and issubclass.

"""

__all__ = ["ABCMeta", "abstractmethod", "abstractproperty", "isinstance",
        "issubclass"]

try:
    # Import abc, and bring isinstance and issubclass into the module namespace
    from abc import ABCMeta, abstractmethod, abstractproperty
    from __builtin__ import isinstance
    from __builtin__ import issubclass

except ImportError:

    # Not python2.6. Use the modifed abc, and redefine isinstance and
    # issubclass so that MetaABC works.
    from _abc import ABCMeta, abstractmethod, abstractproperty

    import __builtin__
    isinstance_builtin = __builtin__.isinstance
    issubclass_builtin = __builtin__.issubclass

    def isinstance(obj, cls):
        try:
            cls = iter(cls)
        except TypeError:
            cls = (cls,)

        ok = 0
        for c in cls:
            if hasattr(c, "__instancecheck__") and\
                    hasattr(c, "__subclasscheck__") and\
                    hasattr(c, "__subclasshook__"):
                ok += c.__instancecheck__(obj)
            else:
                ok += isinstance_builtin(obj, cls)

            if ok: return True

        return False

    def issubclass(obj, cls):
        try:
            cls = iter(cls)
        except TypeError:
            cls = (cls,)

        ok = 0
        for c in cls:
            if hasattr(c, "__subclasscheck__") and\
                    hasattr(c, "__subclasshook__"):
                ok += c.__subclasscheck__(obj)
            else:
                ok += issubclass_builtin(obj, cls)

            if ok: return True

        return False


# End of file
