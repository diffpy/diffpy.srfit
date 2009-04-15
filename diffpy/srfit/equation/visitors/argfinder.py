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

ArgFinder extracts all Arguments from a literal true, even "hidden" ones that
are contained in Partitions within and below Generators.

""" 

from .visitor import Visitor

class ArgFinder(Visitor):
    """ArgFinder extracts Arguments from a Literal tree.

    Attributes
    args    --  The set of collected Arguments (a set, not a list!)
    getconsts  --  Flag indicating whether to grab constant arguments.
    """

    def __init__(self, getconsts = True):
        """Initialize.

        Arguments
        getconsts   --  Flag indicating whether to grap constant arguments
                        (default True).
        
        """
        self.args = set()
        self.getconsts = getconsts
        return

    def reset(self):
        """Reset the argument set."""
        self.args = set()
        return

    def onArgument(self, arg):
        """Process an Argument node."""
        if self.getconsts or not arg.const:
            self.args.add(arg)
        return

    def onOperator(self, op):
        """Process an Operator node."""
        for arg in op.args:
            arg.identify(self)
        return

    def onPartition(self, part):
        """Process a Partition node."""
        for arg in part.args:
            arg.identify(self)
        return

    def onGenerator(self, gen):
        """Process a Generator node."""
        for arg in gen.args:
            arg.identify(self)
        return


# version
__id__ = "$Id$"

#
# End of file
