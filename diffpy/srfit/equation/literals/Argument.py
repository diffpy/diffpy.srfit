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

Arguments are the leaves of an equation tree.
"""

from .Literal import Literal

import numpy

class Argument(Literal):
    """Argument class.
    
    Attributes
    name    --  A name for this Argument.
    clicker --  A Clicker instance for recording change in the value.
    value   --  The value of the Argument. Modified with setValue.
    const   --  A flag indicating whether this is considered a constant.
                Constants may be given special treatment by the Visitors.
    """

    def __init__(self, value = None, name = None, const = False):
        """Initialization."""
        Literal.__init__(self)
        self.value = value
        self.name = name
        self.const = const
        return

    def identify(self, visitor):
        """Identify self to a visitor."""
        visitor.onArgument(self)
        return

    def setValue(self, val):
        """Set the value of the argument to something.

        This will click the clicker.
        """
        if not numpy.array_equiv(val, self.value):
            #print self.name, "changed"
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
