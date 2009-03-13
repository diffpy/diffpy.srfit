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
"""Partition class. 

Partitions are Arguments that have been broken up into logical pieces.
Partitions actually contain a list of Arguments and a 'combine' method that can
un-partition the Arguments.

Operations operate on each Argument of a Partition separately, and perhaps
conditionally based on the tags of the operator. Since multiple conditional
operators can be applied to Partition, the partitioning is preserved throughout
a calculation.  However, once all conditional operations in a Literal tree have
operated on a Partition, it will be combined so subsequent operations will be
quicker. The exception to this rule is when two or more Partitions must act as
arguments in the same Operator. In this case, each Partition in the Operator
will be combined before the operation. In addition, Operators have a 'combine'
flag that can be set to deliberatly combine a Partition after the operation.
"""

from .Literal import Literal

import numpy

class Partition(Literal):
    """Abstract class for Argument partitions.
    
    Attributes
    name    --  A name for this Partition.
    clicker --  A Clicker instance for recording change in any contained
                Argument.
    args    --  List of Arguments, modfied by the 'addArgument' method. 
    tags    --  Set of tags from all arguments.
    tagmap  --  A map of tags to lists of Argument indicies from self.args.
    """ 

    def __init__(self):
        """Initialization."""
        Literal.__init__(self)
        self.name = ""
        self.args = []
        self.tags = set()
        self.tagmap = {}
        # values of arguments, used by Evaluator
        self._partvals = []
        return

    def identify(self, visitor):
        """Identify self to a visitor."""
        visitor.onPartition(self)
        return

    def addArgument(self, arg, *tags):
        """Add an argument to this partition.

        arg     --  Instance of Argument.
        All remaining method arguments are interpreted as tags.
        """
        self.args.append(arg)
        self.tags.update(tags)
        for tag in tags:
            tagset = self.tagmap.get(tag, [])
            tagset.append(len(self.args)-1)
            self.tagmap[tag] = tagset
        self.clicker.addSubject(arg.clicker)
        return

    def combine(self, vals):
        """Combine the values in the partition.

        This is the 'sum' method by default.
        
        values  --  List of values to operate on.

        Returns the combined values
        """
        return sum(vals)

    def _prepare(self):
        """Called by the Evaluator visitor to fill the _partvals list."""
        self._partvals = [arg.value for arg in self.args]
        return


    def __str__(self):
        if self.name:
            return "Argument(" + self.name + ")"
        return self.__repr__()

# version
__id__ = "$Id$"

#
# End of file
