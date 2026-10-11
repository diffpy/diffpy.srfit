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


# End of file
