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
"""Abstract Base Classes for Literals."""

__all__ = ["isinstance", "issubclass", "LiteralABC", "ArgumentABC",
        "OperatorABC"]

from diffpy.srfit.util.abcport import ABCMeta
from diffpy.srfit.util.abcport import abstractmethod, abstractproperty
from diffpy.srfit.util.abcport import isinstance, issubclass

class LiteralABC(object):
    """Abstract Base Class for Literal. See Literal for usage."""

    __metaclass__ = ABCMeta

    @abstractmethod
    def identify(self, visitor): pass

    @abstractmethod
    def getValue(self): pass

    name = abstractproperty(None, None)

# End class LiteralABC

class ArgumentABC(LiteralABC):
    """Abstract Base Class for Argument. See Argument for usage."""

    @abstractmethod
    def setValue(self, value): pass

    const = abstractproperty(None, None)
    value = abstractproperty(None, None)

# End class ArgumentABC

class OperatorABC(LiteralABC):
    """Abstract Base Class for Operator. See Operator for usage."""

    @abstractmethod
    def addLiteral(self, literal): pass

    args = abstractproperty(None, None)
    nin = abstractproperty(None, None)
    nout = abstractproperty(None, None)
    operation = abstractproperty(None, None)
    symbol = abstractproperty(None, None)
    value = abstractproperty(None, None)

# End class OperatorABC
