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

    This class inherits from Literal. See the Literal documentation.

    """

    def _identify(self, visitor):
        """Identify self to a visitor."""
        visitor.onArgument(self)
        return

    def __str__(self):
        if self.name:
            return "Argument(" + self.name + ")"
        return self.__repr__()

# version
__id__ = "$Id$"

#
# End of file
