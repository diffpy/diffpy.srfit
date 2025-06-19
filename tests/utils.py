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
