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
"""Visitor for extracting the Argument entries in a Literal tree.

ArgFinder extracts all Arguments from a literal tree.
"""

__all__ = ["ArgFinder"]

from diffpy.srfit.equation.visitors.visitor import Visitor


class ArgFinder(Visitor):
    """ArgFinder extracts Arguments from a Literal tree.

    Parameters
    ----------
    getconsts : bool, optional
        The flag that is True (default) to also collect constant
        Arguments.

    Attributes
    ----------
    args : list
        The collected Arguments.
    getconsts : bool
        The flag that is True to also collect constant Arguments.
    """

    def __init__(self, getconsts=True):
        self.args = []
        self.getconsts = getconsts
        return

    def reset(self):
        """Reset the argument set."""
        self.args = []
        return

    def on_argument(self, arg):
        """Collect an Argument node.

        Parameters
        ----------
        arg : Argument
            The Argument to collect.

        Returns
        -------
        list
            The Arguments collected so far.
        """
        if self.getconsts or not arg.const:
            self.args.append(arg)
        return self.args

    def on_operator(self, op):
        """Collect the Arguments of an Operator node.

        Parameters
        ----------
        op : Operator
            The Operator whose Arguments are collected.

        Returns
        -------
        list
            The Arguments collected so far.
        """
        for arg in op.args:
            arg.identify(self)
        return self.args


# End of file
