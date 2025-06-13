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
"""Visitors that perform on Literal networks.

Visitors are designed to traverse and extract information from Literal
networks (diffpy.srfit.equation.literals). Visitors are used to
validate, print and extracting Arguments from Literal networks.

The Literal-Visitor relationship is that described by the Visitor
pattern (

http://en.wikipedia.org/wiki/Visitor_pattern).
"""

from __future__ import print_function

from diffpy.srfit.equation.visitors.argfinder import ArgFinder
from diffpy.srfit.equation.visitors.printer import Printer
from diffpy.srfit.equation.visitors.swapper import Swapper
from diffpy.srfit.equation.visitors.validator import Validator


def getArgs(literal, getconsts=True):
    """Get the Arguments of a Literal tree.

    getconsts   --  If True (default), then Arguments designated as constant
                    are also retrieved.

    Returns a list of Arguments searched for depth-first.
    """
    v = ArgFinder(getconsts)
    return literal.identify(v)


def getExpression(literal, eqskip=None):
    """Get math expression string from the Literal tree object.

    eqskip  --  regular expression pattern for Equation objects that should
                be printed under their own name.  Expand all Equation objects
                when None.

    Return string.
    """
    v = Printer()
    v.eqskip = eqskip
    rv = literal.identify(v)
    return rv


def prettyPrint(literal):
    """Print a Literal tree."""
    print(getExpression(literal))
    return


def validate(literal):
    """Validate a Literal tree.

    Raises ValueError if the tree contains errors.
    """
    v = Validator()
    errors = literal.identify(v)
    if errors:
        m = "Errors found in Literal tree '%s'\n" % literal
        m += "\n".join(errors)
        raise ValueError(m)
    return


def swap(literal, oldlit, newlit):
    """Swap one literal for another in a Literal tree.

    Corrections are done in-place unless literal is oldlit, in which
    case the return value is newlit.

    Returns the literal tree with oldlit swapped for newlit.
    """
    if literal is oldlit:
        return newlit
    v = Swapper(oldlit, newlit)
    literal.identify(v)
    return literal
