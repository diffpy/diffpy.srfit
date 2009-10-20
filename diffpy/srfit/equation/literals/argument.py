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

    Attributes
    const   --  A flag indicating whether this is a constant.

    """

    def __init__(self, value = None, const = False):
        Literal.__init__(self, value)
        self.const = const
        return

# version
__id__ = "$Id$"

#
# End of file
