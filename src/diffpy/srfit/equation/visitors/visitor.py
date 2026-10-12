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

from diffpy.utils._deprecator import build_deprecation_message, deprecated

visitor_base = "diffpy.srfit.equation.visitors.visitor.Visitor"
removal_version = "4.0.0"
onArgument_dep_msg = build_deprecation_message(
    visitor_base, "onArgument", "on_argument", removal_version
)
onOperator_dep_msg = build_deprecation_message(
    visitor_base, "onOperator", "on_operator", removal_version
)
onEquation_dep_msg = build_deprecation_message(
    visitor_base, "onEquation", "on_equation", removal_version
)


class Visitor(object):
    """Abstract class for all visitors to a Literal tree.

    Subclasses must override ``on_argument`` and ``on_operator``. See
    implemented visitors for examples of use.
    """

    def on_argument(self, arg):
        """Process an Argument node.

        Parameters
        ----------
        arg : Argument
            The Argument to process.

        Raises
        ------
        NotImplementedError
            Always, because subclasses must override this method.
        """
        return self._abstract("on_argument")

    def on_operator(self, op):
        """Process an Operator node.

        Parameters
        ----------
        op : Operator
            The Operator to process.

        Raises
        ------
        NotImplementedError
            Always, because subclasses must override this method.
        """
        return self._abstract("on_operator")

    def on_equation(self, eq):
        """Process an Equation node.

        Equations are specialized Operators. They don't need to be
        specifically supported by a Visitor, so by default they are
        passed to ``on_operator``.

        Parameters
        ----------
        eq : Equation
            The Equation to process.

        Returns
        -------
        object
            The result of ``on_operator(eq)``.
        """
        return self.on_operator(eq)

    @deprecated(onArgument_dep_msg)
    def onArgument(self, arg):
        """This function has been deprecated and will be removed in
        version 4.0.0.

        Please use
        diffpy.srfit.equation.visitors.visitor.Visitor.on_argument
        instead.
        """
        return self.on_argument(arg)

    @deprecated(onOperator_dep_msg)
    def onOperator(self, op):
        """This function has been deprecated and will be removed in
        version 4.0.0.

        Please use
        diffpy.srfit.equation.visitors.visitor.Visitor.on_operator
        instead.
        """
        return self.on_operator(op)

    @deprecated(onEquation_dep_msg)
    def onEquation(self, eq):
        """This function has been deprecated and will be removed in
        version 4.0.0.

        Please use
        diffpy.srfit.equation.visitors.visitor.Visitor.on_equation
        instead.
        """
        return self.on_equation(eq)

    # throw an exception
    def _abstract(self, method):
        raise NotImplementedError(
            "class '%s' should override method '%s'"
            % (self.__class__.__name__, method)
        )


# End of file
