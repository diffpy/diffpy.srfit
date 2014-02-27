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
TestCaseOptional = object
def isinstalled(*TestCases):
    'check for optional dependencies'
    return not TestCaseOptional in TestCases
def testcase(*TestCases):
    'get appropriate TestCase object for optional dependencies'
    return TestCase if isinstalled(*TestCases) else TestCaseOptional

# Resolve the TestCase*Optional classes
try:
    import sans.pr.invertor
    import sans.models
    TestCaseSaS = TestCase
except ImportError, e:
    TestCaseSaS = TestCaseOptional
    logging.warning('%s, SaS tests skipped.', e)

try:
    import diffpy.Structure
    import pyobjcryst
    import diffpy.srreal
    TestCaseStructure = TestCase
except ImportError, e:
    TestCaseStructure = TestCaseOptional
    logging.warning('%s, Structure tests skipped.', e)

try:
    import diffpy.srreal
    TestCasePdf = TestCase
except ImportError, e:
    TestCasePdf = TestCaseOptional
    logging.warning('%s, Pdf tests skipped.', e)


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
