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
"""Argument class.

Arguments are the leaves of an equation tree, in essence a variable or a
constant.
"""

__all__ = ["Argument"]

from numpy import array_equal, ndarray

from diffpy.srfit.equation.literals.abcs import ArgumentABC
from diffpy.srfit.equation.literals.literal import Literal


class Argument(Literal, ArgumentABC):
    """Argument class.

    Attributes
    name    --  A name for this Argument.
    const   --  A flag indicating whether this is considered a constant.
                Constants may be given special treatment by the Visitors.
    _value  --  The value of the Argument. Modified with 'setValue'.
    value   --  Property for 'getValue' and 'setValue'.
    """

    const = None

    def __init__(self, name=None, value=None, const=False):
        """Initialization."""
        Literal.__init__(self, name)
        self.const = const
        self.value = value
        return

    def identify(self, visitor):
        """Identify self to a visitor."""
        return visitor.onArgument(self)

    def getValue(self):
        """Get the value of this Literal."""
        return self._value

    def setValue(self, val):
        """Set the value of the Literal.

        val --  The value to assign
        """
        if isinstance(self._value, ndarray) or isinstance(val, ndarray):
            notequiv = not array_equal(self._value, val)
        else:
            notequiv = self._value != val
        if notequiv:
            self._value = val
            self.notify()
        return

    value = property(
        lambda self: self.getValue(), lambda self, val: self.setValue(val)
    )


# End of file
