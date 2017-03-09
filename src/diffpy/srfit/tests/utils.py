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

import diffpy.srfit.equation.literals as literals
from diffpy.srfit.sas.sasimport import sasimport
from diffpy.srfit.tests import logger
from unittest import TestCase

# Create a singleton and a test for optional test cases
_TestCaseDisabled = object

def _allactive(*testcaseclasses):
    'True if none of the testcaseclasses has been disabled.'
    return not _TestCaseDisabled in testcaseclasses

def testoptional(*testcaseclasses):
    'Return unittest.TestCase only if all testcaseclasses are active.'
    if _allactive(*testcaseclasses):  return TestCase
    return _TestCaseDisabled

# Resolve the TestCase*Optional classes
try:
    sasimport('sas.pr.invertor')
    sasimport('sas.models')
    TestCaseSaS = TestCase
except ImportError, e:
    TestCaseSaS = _TestCaseDisabled
    logger.warning('%s, SaS tests skipped.', e)

try:
    import diffpy.Structure as m; del m
    TestCaseStructure = TestCase
except ImportError:
    TestCaseStructure = _TestCaseDisabled
    logger.warning('Cannot import diffpy.Structure, Structure tests skipped.')

try:
    import pyobjcryst as m; del m
    TestCaseObjCryst = TestCase
except ImportError:
    TestCaseObjCryst = _TestCaseDisabled
    logger.warning('Cannot import pyobjcryst, pyobjcryst tests skipped.')

try:
    import diffpy.srreal.pdfcalculator as m; del m
    TestCasePDF = TestCase
except ImportError:
    TestCasePDF = _TestCaseDisabled
    logger.warning('Cannot import diffpy.srreal, PDF tests skipped.')


def _makeArgs(num):
    args = []
    for i in xrange(num):
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

# End of file
