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
from diffpy.srfit.tests import logger


# Resolve availability of optional third-party packages.

# srfit-sasview or sasview

try:
    _msg_nosas = "No module named 'sas.pr.invertor'"
    sasimport('sas.pr.invertor')
    _msg_nosas = "No module named 'sas.models'"
    sasimport('sas.models')
    has_sas = True
except ImportError as e:
    has_sas = False
    logger.warning('%s, SaS tests skipped.', e)

# diffpy.structure

_msg_nostructure = "No module named 'diffpy.structure'"
try:
    import diffpy.structure as m; del m
    has_structure = True
except ImportError:
    has_structure = False
    logger.warning('Cannot import diffpy.structure, Structure tests skipped.')

# pyobjcryst

_msg_nopyobjcryst = "No module named 'pyobjcryst'"
try:
    import pyobjcryst as m; del m
    has_pyobjcryst = True
except ImportError:
    has_pyobjcryst = False
    logger.warning('Cannot import pyobjcryst, pyobjcryst tests skipped.')

# diffpy.srreal

_msg_nosrreal = "No module named 'diffpy.srreal'"
try:
    import diffpy.srreal.pdfcalculator as m; del m
    has_srreal = True
except ImportError:
    has_srreal = False
    logger.warning('Cannot import diffpy.srreal, PDF tests skipped.')

# Helper functions for testing -----------------------------------------------

def _makeArgs(num):
    args = []
    for i in range(num):
        j=i+1
        args.append(literals.Argument(name="v%i"%j, value=j))
    return args


def noObserversInGlobalBuilders():
    """True if no observer function leaks to global builder objects.

    Ensure objects are not immortal due to a reference from static value.
    """
    from diffpy.srfit.equation.builder import _builders
    rv = True
    for n, b in _builders.items():
        if b.literal and b.literal._observers:
            rv = False
            break
    return rv


def datafile(filename):
    from pkg_resources import resource_filename
    rv = resource_filename(__name__, "testdata/" + filename)
    return rv


def capturestdout(f, *args, **kwargs):
    """Capture the standard output from a call of function f.
    """
    savestdout = sys.stdout
    fp = six.StringIO()
    try:
        sys.stdout = fp
        f(*args, **kwargs)
    finally:
        sys.stdout = savestdout
    return fp.getvalue()

# End of file
