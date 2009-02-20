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
    """Abstract class for all visitors to a literal tree.

    The methods are considered self-documented. See implemented visitors for
    examples of use.
    """

    # Variables
    def onVariable(self, var):
        return self._abstract("onVariable")

    # Operators
    def onOperator(self, op):
        return self._abstract("onOperator")


    # throw an exception
    def _abstract(self, method):
        raise NotImplementedError(
            "class '%s' should override method '%s'" % (self.__class__.__name__, method))


# version
__id__ = "$Id$"

#
# End of file
