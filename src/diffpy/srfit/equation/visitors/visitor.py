#!/usr/bin/env python
##############################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2008 The Trustees of Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Chris Farrow
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
##############################################################################
"""Base visitor class.

Visitors work with Literal trees to perform a specified action. They are
designed according to the visitor pattern (
http://en.wikipedia.org/wiki/Visitor_pattern).

Visitors work with the following Literal classes
(diffpy.srfit.equation.literals):
Argument
Operator

See the Visitor class for the required methods that each Visitor must overload.
"""

__all__ = ["Visitor"]


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

    def onEquation(self, eq):
        """Process an Equation node.

        Equations are specialized Operators. They don't need to be
        specifically supported by a Visitor.
        """
        return self.onOperator(eq)

    # throw an exception
    def _abstract(self, method):
        raise NotImplementedError(
            "class '%s' should override method '%s'"
            % (self.__class__.__name__, method)
        )


# End of file
