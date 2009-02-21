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

from diffpy.srfit.equation.literals.Literal import Literal

class Argument(Literal):
    """Argument class."""

    def identify(self, visitor):
        """Identify self to a visitor."""
        visitor.onArgument(self)
        return

    def setValue(self, val):
        """Set the value of the argument to something.

        This will click the clicker.
        """
        self.value = val
        self.clicker.click()
        return

# version
__id__ = "$Id$"

#
# End of file
