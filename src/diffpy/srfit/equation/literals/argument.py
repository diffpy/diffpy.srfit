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
"""Argument class.

Arguments are the leaves of an equation tree, in essence a variable or a
constant.
"""

__all__ = ["Argument"]

from numpy import array_equal, ndarray

from diffpy.srfit.equation.literals.abcs import ArgumentABC
from diffpy.srfit.equation.literals.literal import Literal


class Argument(Literal, ArgumentABC):
    """A leaf of an equation tree, holding a variable or a constant.

    Parameters
    ----------
    name : str, optional
        The name of the Argument. Default is None.
    value : object, optional
        The initial value of the Argument. Default is None.
    const : bool, optional
        The flag that is True if the Argument is a constant. Default is
        False.

    Attributes
    ----------
    name : str or None
        The name of the Argument.
    const : bool
        The flag that is True if the Argument is a constant. Visitors
        may treat constants differently.
    value : object
        The value of the Argument, read with ``get_value`` and changed
        with ``set_value``.
    """

    const = None

    def __init__(self, name=None, value=None, const=False):
        Literal.__init__(self, name)
        self.const = const
        self.value = value
        return

    def identify(self, visitor):
        """Identify this Argument to a visitor.

        Parameters
        ----------
        visitor : Visitor
            The visitor to identify to.

        Returns
        -------
        object
            The result of ``visitor.on_argument(self)``.
        """
        return visitor.on_argument(self)

    def get_value(self):
        """Get the value of this Argument.

        Returns
        -------
        object
            The current value of the Argument.
        """
        return self._value

    def set_value(self, val):
        """Set the value of this Argument.

        Observers are notified only if the new value differs from the
        current one.

        Parameters
        ----------
        val : object
            The value to assign.
        """
        if isinstance(self._value, ndarray) or isinstance(val, ndarray):
            notequiv = not array_equal(self._value, val)
        else:
            notequiv = self._value != val
        if notequiv:
            self._value = val
            self.notify()
        return

    value = property(
        lambda self: self.get_value(), lambda self, val: self.set_value(val)
    )


# End of file
