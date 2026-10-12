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

from diffpy.srfit.equation.visitors.argfinder import ArgFinder
from diffpy.srfit.equation.visitors.printer import Printer
from diffpy.srfit.equation.visitors.swapper import Swapper
from diffpy.srfit.equation.visitors.validator import Validator
from diffpy.utils._deprecator import build_deprecation_message, deprecated

visitors_base = "diffpy.srfit.equation.visitors"
removal_version = "4.0.0"
getArgs_dep_msg = build_deprecation_message(
    visitors_base, "getArgs", "get_args", removal_version
)
getExpression_dep_msg = build_deprecation_message(
    visitors_base, "getExpression", "get_expression", removal_version
)
prettyPrint_dep_msg = build_deprecation_message(
    visitors_base, "prettyPrint", "pretty_print", removal_version
)


def get_args(literal, getconsts=True):
    """Get the Arguments of a Literal tree.

    Parameters
    ----------
    literal : Literal
        The root of the Literal tree.
    getconsts : bool, optional
        The flag that is True (default) to also get the Arguments that
        are constant.

    Returns
    -------
    list of Argument
        The Arguments, found depth-first.
    """
    v = ArgFinder(getconsts)
    return literal.identify(v)


def get_expression(literal, eqskip=None):
    """Get the math expression string of a Literal tree.

    Parameters
    ----------
    literal : Literal
        The root of the Literal tree.
    eqskip : str, optional
        The regular expression pattern for Equation objects that should
        be written under their own name. If None (default), all
        Equation objects are expanded.

    Returns
    -------
    str
        The expression of the Literal tree.
    """
    v = Printer()
    v.eqskip = eqskip
    rv = literal.identify(v)
    return rv


def pretty_print(literal):
    """Print the math expression of a Literal tree.

    Parameters
    ----------
    literal : Literal
        The root of the Literal tree.
    """
    print(get_expression(literal))
    return


def validate(literal):
    """Validate a Literal tree.

    Parameters
    ----------
    literal : Literal
        The root of the Literal tree.

    Raises
    ------
    ValueError
        If the tree contains errors.
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

    Parameters
    ----------
    literal : Literal
        The root of the Literal tree.
    oldlit : Literal
        The literal to replace.
    newlit : Literal
        The literal to put in its place.

    Returns
    -------
    Literal
        The Literal tree with ``oldlit`` swapped for ``newlit``.
    """
    if literal is oldlit:
        return newlit
    v = Swapper(oldlit, newlit)
    literal.identify(v)
    return literal


@deprecated(getArgs_dep_msg)
def getArgs(literal, getconsts=True):
    """This function has been deprecated and will be removed in version
    4.0.0.

    Please use diffpy.srfit.equation.visitors.get_args instead.
    """
    return get_args(literal, getconsts)


@deprecated(getExpression_dep_msg)
def getExpression(literal, eqskip=None):
    """This function has been deprecated and will be removed in version
    4.0.0.

    Please use diffpy.srfit.equation.visitors.get_expression instead.
    """
    return get_expression(literal, eqskip)


@deprecated(prettyPrint_dep_msg)
def prettyPrint(literal):
    """This function has been deprecated and will be removed in version
    4.0.0.

    Please use diffpy.srfit.equation.visitors.pretty_print instead.
    """
    return pretty_print(literal)
