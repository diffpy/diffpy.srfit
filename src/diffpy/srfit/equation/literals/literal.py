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
"""Literal base class used to construct equation trees.

Literals are base pieces of the equation hierarchy. The 'identify'
method identifies the Literal to a visitor by calling the identifying
method of the visitor.
"""

__all__ = ["Literal"]

from diffpy.srfit.equation.literals.abcs import LiteralABC
from diffpy.srfit.util.observable import Observable
from diffpy.utils._deprecator import build_deprecation_message, deprecated

literal_base = "diffpy.srfit.equation.literals.Literal"
removal_version = "4.0.0"
getValue_dep_msg = build_deprecation_message(
    literal_base, "getValue", "get_value", removal_version
)


class Literal(Observable, LiteralABC):
    """Abstract class for equation pieces, such as operators and
    arguments.

    Literal derives from Observable. See diffpy.srfit.util.observable.

    Parameters
    ----------
    name : str, optional
        The name of the Literal. If None (default), the class attribute
        ``name`` is used.

    Attributes
    ----------
    name : str or None
        The name of the Literal.
    _value : object
        The cached value of the Literal, or None if it must be
        recomputed.
    """

    name = None
    _value = None

    def __init__(self, name=None):
        Observable.__init__(self)
        if name is not None:
            self.name = name
        return

    def get_value(self):
        """Get the value of the Literal.

        Raises
        ------
        NotImplementedError
            Always, because derived classes must define this method.
        """
        raise NotImplementedError("Define in derived class")

    @deprecated(getValue_dep_msg)
    def getValue(self):
        """This function has been deprecated and will be removed in
        version 4.0.0.

        Please use diffpy.srfit.equation.literals.Literal.get_value
        instead.
        """
        return self.get_value()

    def identify(self, visitor):
        """Identify this Literal to a visitor.

        Parameters
        ----------
        visitor : Visitor
            The visitor to identify to.

        Raises
        ------
        NotImplementedError
            Always, because derived classes must define this method.
        """
        m = "'%s' must override 'identify'" % self.__class__.__name__
        raise NotImplementedError(m)

    def _flush(self, other):
        """Clear the cached value and notify observers."""
        if self._value is None:
            return
        self._value = None
        self.notify(other)
        return

    def __str__(self):
        return "%s(%s)" % (self.__class__.__name__, self.name)


# End of file
