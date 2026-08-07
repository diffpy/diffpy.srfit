#!/usr/bin/env python
##############################################################################
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
##############################################################################
"""Abstract Base Classes for Literals."""

__all__ = ["LiteralABC", "ArgumentABC", "OperatorABC"]

from abc import ABC, abstractmethod


class LiteralABC(ABC):
    """Abstract Base Class for Literal.

    See Literal for usage.
    """

    @abstractmethod
    def identify(self, visitor):
        """Identify this literal using a visitor."""
        pass

    @abstractmethod
    def getValue(self):
        """Return the value of the literal."""
        pass

    @property
    @abstractmethod
    def name(self):
        """Name of the literal."""
        pass


# End class LiteralABC


class ArgumentABC(LiteralABC):
    """Abstract Base Class for Argument.

    See Argument for usage.
    """

    @abstractmethod
    def set_value(self, value):
        """Set the value of the argument."""
        pass

    @property
    @abstractmethod
    def const(self):
        """Whether the argument is constant."""
        pass

    @property
    @abstractmethod
    def value(self):
        """Value of the argument."""
        pass


# End class ArgumentABC


class OperatorABC(LiteralABC):
    """Abstract Base Class for Operator.

    See Operator for usage.
    """

    @abstractmethod
    def addLiteral(self, literal):
        """Add a literal argument to the operator."""
        pass

    @property
    @abstractmethod
    def args(self):
        """Arguments of the operator."""
        pass

    @property
    @abstractmethod
    def nin(self):
        """Number of input arguments."""
        pass

    @property
    @abstractmethod
    def nout(self):
        """Number of outputs."""
        pass

    @property
    @abstractmethod
    def operation(self):
        """Callable implementing the operator."""
        pass

    @property
    @abstractmethod
    def symbol(self):
        """Symbol representing the operator."""
        pass

    @property
    @abstractmethod
    def value(self):
        """Value produced by the operator."""
        pass


# End class OperatorABC
