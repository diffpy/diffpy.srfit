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

"""Literal base class used to construct equation trees.

Literals are base pieces of the equation hierarchy. The 'identify' method
identifies the Literal to a visitor by calling the identifying method of the
vistior.
"""

__all__ = ["Literal"]

from diffpy.srfit.equation.literals.abcs import LiteralABC
from diffpy.srfit.util.observable import Observable

class Literal(Observable,LiteralABC):
    """Abstract class for equation pieces, such as operators and arguments.

    Literal derives from Observable. See diffpy.srfit.util.observable.

    Attributes
    name    --  A name for this Literal (default None).
    _value  --  The value of the Literal.

    """

    name = None
    _value = None

    def __init__(self, name=None):
        """Initialization."""
        Observable.__init__(self)
        if name is not None:
            self.name = name
        return

    def getValue(self):
        """Get the value of the Literal."""
        raise NotImplementedError("Define in derived class")

    def identify(self, visitor):
        """Identify self to a visitor."""
        m = "'%s' must override 'identify'" % self.__class__.__name__
        raise NotImplementedError(m)

    def _flush(self, other):
        """Invalidate my state and notify observers."""
        if self._value is None:
            return
        self._value = None
        self.notify(other)
        return

    def __str__(self):
        return "%s(%s)"%(self.__class__.__name__, self.name)

# End of file
