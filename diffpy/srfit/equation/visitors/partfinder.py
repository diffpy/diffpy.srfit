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
"""Visitor for determining if Literal tree has any Partitions.

The PartFinder has a parts list that is filled with Partitions in the Literal
tree. It detects Partitions that are created by Generators, provided that the
Generator has it's literal attribute defined as a Partition.

""" 

from ..literals import Partition

from .visitor import Visitor

class PartFinder(Visitor):
    """PartFinder determines if a Literal tree has Partitions.

    Attributes
    parts   --  A list of Partitions from the tree.
    """

    def __init__(self):
        """Initialize."""
        self.parts = []
        return

    def reset(self):
        """Reset the Partition list."""
        self.parts = []
        return

    def onArgument(self, arg):
        """Process an Argument node."""
        return

    def onOperator(self, op):
        """Process an Operator node."""
        for arg in op.args:
            arg.identify(self)
        return

    def onPartition(self, part):
        """Process a Partition node."""
        self.parts.append(part)
        return

    def onGenerator(self, gen):
        """Process a Generator node."""
        if isinstance(gen.literal, Partition):
            self.parts.append(gen.literal)
        return


# version
__id__ = "$Id$"

#
# End of file
