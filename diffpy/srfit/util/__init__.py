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


"""\
Utilities and constants used throughout SrFit.
"""

_DASHEDLINE = 78 * '-'


def sortKeyForNumericString(s):
    """\
    Compute key for sorting strings according to their integer numeric value.

    Each string gets split to string and integer segments to create keys
    for comparison.  Signs, decimal points and exponents are ignored.
    This function is intended as the ``key`` argument for the ``sorted``
    or ``list.sort`` function.

    Parameters
    ----------
    s : str
        String which may have numeric components, e.g., "a12b".

    Returns
    -------
    tuple
        Tuple of non-numeric segments intermixed with integer values.
    """
    if sortKeyForNumericString._rx is None:
        import re
        sortKeyForNumericString._rx = re.compile(r'(\d+)')
    rx = sortKeyForNumericString._rx
    rv = tuple((int(w) if i % 2 else w)
               for i, w in enumerate(rx.split(s)))
    return rv

sortKeyForNumericString._rx = None

# End of file
