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
"""Literal class. 

Literals are base pieces of the equation hierarchy. The have a single method,
'identify' that identifies the Literal to a visitor. Literal has a single data
attribute, 'clicker', that records changes in the state of a Literal. Clicker is
found in the diffpy.equation.Clicker module.
"""

from diffpy.srfit.equation.Clicker import Clicker

class Literal(object):
    """Abstract class for equation pieces, such as operators and variables.

    This is the base class for the equation data structure. The structure is
    biased heavily towards the evaluation of the equation tree, and so some
    structure required to make that fast and efficient is contained in the base
    class, breaking the encapsulation of responsibility.
    """

    def __init__(self):
        """Initialization."""
        self.clicker = Clicker()
        return

    def identify(self, visitor):
        """Identify self to a visitor."""
        m = "'%s' must override 'identify'" % self.__class__.__name__
        raise NotImplementedError(m)


# version
__id__ = "$Id$"

#
# End of file
