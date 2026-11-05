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
"""Utilities and constants used throughout SrFit."""

import math
import re

from diffpy.utils._deprecator import build_deprecation_message, deprecated

_DASHEDLINE = 78 * "-"

sortKeyForNumericString_dep_msg = build_deprecation_message(
    "diffpy.srfit.util",
    "sortKeyForNumericString",
    "sort_key_for_numeric_string",
    "4.0.0",
)


def sort_key_for_numeric_string(string):
    """Return a sort key that orders the numbers in a string by value.

    The string is split into text and integer segments, so ``"a2"`` sorts
    before ``"a10"``. Signs, decimal points and exponents are ignored. Use it
    as the ``key`` argument of ``sorted`` or ``list.sort``.

    Parameters
    ----------
    string : str
        The string to build a key for, which may contain numbers, e.g.
        ``"a12b"``.

    Returns
    -------
    tuple
        The text segments of ``string`` with integer values in between,
        e.g. ``("a", 12, "b")``.
    """
    segments = re.split(r"(\d+)", string)
    sort_key = tuple(
        int(segment) if segment.isdecimal() else segment
        for segment in segments
    )
    return sort_key


@deprecated(sortKeyForNumericString_dep_msg)
def sortKeyForNumericString(s):
    """This function has been deprecated and will be removed in version
    4.0.0.

    Please use diffpy.srfit.util.sort_key_for_numeric_string instead.
    """
    return sort_key_for_numeric_string(s)


def format_value_with_uncertainty(value, uncertainty):
    """Return a value and its uncertainty in parenthesis notation.

    The uncertainty is rounded to two significant digits, and the value is
    rounded to the same decimal place, so ``3.19965644 +/- 0.0021228`` gives
    ``"3.1997(21)"``. An uncertainty that is zero, negative, infinite or not
    a number gives the value alone, with six significant digits.

    Parameters
    ----------
    value : float
        The value to format.
    uncertainty : float
        The uncertainty in ``value``.

    Returns
    -------
    str
        The value with the digits of its uncertainty in parentheses, e.g.
        ``"3.1997(21)"``.
    """
    if not math.isfinite(uncertainty) or uncertainty <= 0:
        formatted_value = f"{value:.6g}"
        return formatted_value
    decimals = 1 - math.floor(math.log10(uncertainty))
    if round(uncertainty * 10**decimals) >= 100:
        decimals -= 1
    decimals = max(0, decimals)
    uncertainty_digits = round(uncertainty * 10**decimals)
    formatted_value = f"{value:.{decimals}f}({uncertainty_digits})"
    return formatted_value


# End of file
