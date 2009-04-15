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
"""Argument class. 

Arguments are the leaves of an equation tree, in essense a variable or a
constant.
"""

from .literal import Literal

import numpy

class Argument(Literal):
    """Argument class.
    
    Attributes
    name    --  A name for this Argument.
    clicker --  A Clicker instance for recording change in the value.
    const   --  A flag indicating whether this is considered a constant.
                Constants may be given special treatment by the Visitors.
    value   --  The value of the Argument. Modified with setValue.
    """

    def __init__(self, value = None, name = None, const = False):
        """Initialization."""
        Literal.__init__(self)
        self.name = name
        self.const = const
        self.value = None
        self.setValue(value)
        return

    def identify(self, visitor):
        """Identify self to a visitor."""
        visitor.onArgument(self)
        return

    def setValue(self, val):
        """Set the value of the argument to something.

        This will click the clicker.
        """
        # This faster than using numpy.array_equiv.
        # The most common case is that of comparing scalars so we check that
        # first.
        notequiv = (val != self.value)
        if notequiv is True:
            self.value = val
            self.clicker.click()
        elif notequiv is False:
            return
        elif notequiv.any():
            self.value = val
            self.clicker.click()
        return

    def __str__(self):
        if self.name:
            return "Argument(" + self.name + ")"
        return self.__repr__()

# version
__id__ = "$Id$"

#
# End of file
