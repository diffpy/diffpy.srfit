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

class Visitor(object):
    """Abstract class for all visitors to a Literal tree.

    See implemented visitors for examples of use.
    """

    def onArgument(self, arg):
        """Process an Argument node."""
        return self._abstract("onArgument")

    def onOperator(self, op):
        """Process an Operator node."""
        return self._abstract("onOperator")

    def onPartition(self, part):
        """Process a Partition node."""
        return self._abstract("onPartition")


    # throw an exception
    def _abstract(self, method):
        raise NotImplementedError(
            "class '%s' should override method '%s'" % (self.__class__.__name__, method))


# version
__id__ = "$Id$"

#
# End of file
