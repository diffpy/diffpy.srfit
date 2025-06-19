#!/usr/bin/env python
##############################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2010 The Trustees of Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
##############################################################################
"""Helper routines for testing."""

import sys

import six

import diffpy.srfit.equation.literals as literals
from diffpy.srfit.sas.sasimport import sasimport
from tests import logger

# Helper functions for testing -----------------------------------------------


def _makeArgs(num):
    args = []
    for i in range(num):
        j = i + 1
        args.append(literals.Argument(name="v%i" % j, value=j))
    return args


def noObserversInGlobalBuilders():
    """True if no observer function leaks to global builder objects.

    Ensure objects are not immortal due to a reference from static
    value.
    """
    from diffpy.srfit.equation.builder import _builders

    rv = True
    for n, b in _builders.items():
        if b.literal and b.literal._observers:
            rv = False
            break
    return rv


def capturestdout(f, *args, **kwargs):
    """Capture the standard output from a call of function f."""
    savestdout = sys.stdout
    fp = six.StringIO()
    try:
        sys.stdout = fp
        f(*args, **kwargs)
    finally:
        sys.stdout = savestdout
    return fp.getvalue()


# End of file
