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

"""Unit tests for diffpy.srfit.
"""
from __future__ import print_function
import six

def testsuite():
    '''Build a unit tests suite for the diffpy.srfit package.

    Return a unittest.TestSuite object.
    '''
    import unittest
    modulenames = '''
        testbuilder
        testcharacteristicfunctions
        testconstraint
        testcontribution
        testdiffpyparset
        testequation
        testfitrecipe
        testfitresults
        testliterals
        testobjcrystparset
        testordereddict
        testparameter
        testparameterset
        testpdf
        testprofile
        testprofilegenerator
        testrecipeorganizer
        testrestraint
        testsas
        testsgconstriants
        testtagmanager
        testvisitors
    '''.split()
    suite = unittest.TestSuite()
    loader = unittest.defaultTestLoader
    return unittest.defaultTestLoader.loadTestsFromNames(
        modulenames)


def test():
    '''Execute all unit tests for the diffpy.srfit package.
    Return a unittest TestResult object.
    '''
    import unittest
    suite = testsuite()
    runner = unittest.TextTestRunner()
    result = runner.run(suite)
    return result


# End of file
