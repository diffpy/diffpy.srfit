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
"""Visitor for extracting the Argument entries in a Literal tree.
""" 

from .Visitor import Visitor

class ArgFinder(Visitor):
    """ArgFinder extracts the Argument entries from a Literal tree.

    Attributes
    args    --  The set of collected Arguments (a set, not a list!)
    """

    def __init__(self):
        """Initialize."""
        self.args = set()
        return

    def onArgument(self, arg):
        """Process an Argument node."""
        self.args.add(arg)
        return

    def onOperator(self, op):
        """Process an Operator node."""
        for arg in op.args:
            arg.identify(self)
        return


# version
__id__ = "$Id$"

#
# End of file
