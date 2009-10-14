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

ArgFinder extracts all Arguments from a literal true.

""" 

from .visitor import Visitor

from diffpy.srfit.util.ordereddict import OrderedDict

def getArguments(literal, getconsts = True):
    """Get an OrderedDict of arguments from an evaluation network."""
    f = ArgFinder(getconsts)
    literal.identify(f)
    args = OrderedDict([(a.name, a) for a in f.args])
    return args

class ArgFinder(Visitor):
    """ArgFinder extracts Arguments from a Literal tree.

    Attributes
    args    --  The list of collected Arguments
    getconsts  --  Flag indicating whether to grab constant arguments.

    """

    def __init__(self, getconsts = True):
        """Initialize.

        Arguments
        getconsts   --  Flag indicating whether to grap constant arguments
                        (default True).
        
        """
        self.args = []
        self.getconsts = getconsts
        return

    def reset(self):
        """Reset the argument set."""
        self.args = []
        return

    def onArgument(self, arg):
        """Process an Argument node."""
        print arg, arg.getValue()
        if self.getconsts or not arg.const:
            self.args.append(arg)
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
