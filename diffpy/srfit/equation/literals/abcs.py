#!/usr/bin/env python
########################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2008 Trustees of the Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Chris Farrow
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
########################################################################
"""Abstract Base Classes for Literals."""

__all__ = ["LiteralABC", "ArgumentABC", "GeneratorABC", "OperatorABC",
        "PartitionABC", "isinstance", "issubclass"]

from diffpy.srfit.util.abc import *

class LiteralABC(object):
    """Abstract Base Class for Literal. See Literal for usage."""

    __metaclass__ = ABCMeta

    @abstractmethod
    def identify(self, visitor): pass

    name = abstractproperty(None, None)
    clicker = abstractproperty(None, None)

# End class LiteralABC

class ArgumentABC(LiteralABC):
    """Abstract Base Class for Argument. See Argument for usage."""

    @abstractmethod
    def setValue(self, value): pass

    @abstractmethod
    def getValue(self): pass

    const = abstractproperty(None, None)
    value = abstractproperty(None, None)

# End class ArgumentABC


class GeneratorABC(LiteralABC):
    """Abstract Base Class for Generator. See Generator for usage."""

    @abstractmethod
    def addLiteral(self, literal): pass

    @abstractmethod
    def generate(self, clicker): pass

    args = abstractproperty(None, None)
    literal = abstractproperty(None, None)

# End class GeneratorABC


class OperatorABC(LiteralABC):
    """Abstract Base Class for Operator. See Operator for usage."""

    @abstractmethod
    def addLiteral(self, literal): pass

    @abstractmethod
    def setCombine(self, combine): pass

    @abstractmethod
    def addTags(self, *tags): pass

    @abstractmethod
    def clearTags(self): pass

    args = abstractproperty(None, None)
    nin = abstractproperty(None, None)
    nout = abstractproperty(None, None)
    operation = abstractproperty(None, None)
    symbol = abstractproperty(None, None)
    _tags = abstractproperty(None, None)
    _cancombine = abstractproperty(None, None)
    _proxy = abstractproperty(None, None)

# End class OperatorABC


class PartitionABC(LiteralABC):
    """Abstract Base Class for Partition. See Partition for usage."""

    @abstractmethod
    def addArgument(self, arg, *tags): pass

    @abstractmethod
    def combine(self, vals): pass

    @abstractmethod
    def _prepare(self, *tags): pass

    args = abstractproperty(None, None)
    tags = abstractproperty(None, None)
    tagmap = abstractproperty(None, None)
    _partvals = abstractproperty(None, None)

# End class PartitionABC
