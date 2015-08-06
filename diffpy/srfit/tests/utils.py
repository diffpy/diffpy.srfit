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

from __future__ import print_function
import six

import logging
import diffpy.srfit.equation.literals as literals
from unittest import TestCase

# python2/3 compatible xrange. xrange was renamed to range in python 3 and
# range was removed
try:
    xrange
except NameError:
    xrange = range

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
except ImportError as ie:
    TestCaseSaS = _TestCaseDisabled
    logging.warning('%s, SaS tests skipped.', ie)

try:
    import diffpy.Structure
    TestCaseStructure = TestCase
except ImportError:
    TestCaseStructure = _TestCaseDisabled
    logging.warning('Cannot import diffpy.Structure, Structure tests skipped.')

try:
    import pyobjcryst
    TestCaseObjCryst = TestCase
except ImportError:
    TestCaseObjCryst = _TestCaseDisabled
    logging.warning('Cannot import pyobjcryst, pyobjcryst tests skipped.')

try:
    import diffpy.srreal.pdfcalculator
    TestCasePDF = TestCase
except ImportError:
    TestCasePDF = _TestCaseDisabled
    logging.warning('Cannot import diffpy.srreal, PDF tests skipped.')


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
