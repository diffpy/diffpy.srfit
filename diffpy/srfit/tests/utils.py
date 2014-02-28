#!/usr/bin/env python
##############################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2010 Trustees of the Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
##############################################################################
"""Helper routines for testing."""

import logging
import diffpy.srfit.equation.literals as literals
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
    import sans.pr.invertor
    import sans.models
    TestCaseSaS = TestCase
except ImportError, e:
    TestCaseSaS = _TestCaseDisabled
    logging.warning('%s, SaS tests skipped.', e)

try:
    import diffpy.Structure
    import pyobjcryst
    import diffpy.srreal
    TestCaseStructure = TestCase
except ImportError, e:
    TestCaseStructure = _TestCaseDisabled
    logging.warning('%s, Structure tests skipped.', e)

try:
    import diffpy.srreal
    TestCasePDF = TestCase
except ImportError, e:
    TestCasePDF = _TestCaseDisabled
    logging.warning('%s, PDF tests skipped.', e)


def _makeArgs(num):
    args = []
    for i in xrange(num):
        j=i+1
        args.append(literals.Argument(name="v%i"%j, value=j))
    return args


def datafile(filename):
    from pkg_resources import resource_filename
    rv = resource_filename(__name__, "testdata/" + filename)
    return rv

# End of file
