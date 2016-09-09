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

# create logger instance for the tests subpackage
import logging

logging.basicConfig()
logger = logging.getLogger(__name__)
del logging


def testsuite():
    '''Build a unit tests suite for the diffpy.srfit package.

    Return a unittest.TestSuite object.
    '''
    import unittest
    modulenames = '''
        diffpy.srfit.tests.testbuilder
        diffpy.srfit.tests.testcharacteristicfunctions
        diffpy.srfit.tests.testconstraint
        diffpy.srfit.tests.testcontribution
        diffpy.srfit.tests.testdiffpyparset
        diffpy.srfit.tests.testequation
        diffpy.srfit.tests.testfitrecipe
        diffpy.srfit.tests.testfitresults
        diffpy.srfit.tests.testliterals
        diffpy.srfit.tests.testobjcrystparset
        diffpy.srfit.tests.testparameter
        diffpy.srfit.tests.testparameterset
        diffpy.srfit.tests.testpdf
        diffpy.srfit.tests.testprofile
        diffpy.srfit.tests.testprofilegenerator
        diffpy.srfit.tests.testrecipeorganizer
        diffpy.srfit.tests.testrestraint
        diffpy.srfit.tests.testsas
        diffpy.srfit.tests.testsgconstriants
        diffpy.srfit.tests.testtagmanager
        diffpy.srfit.tests.testvisitors
        diffpy.srfit.tests.testweakrefcallable
    '''.split()
    suite = unittest.TestSuite()
    loader = unittest.defaultTestLoader
    mobj = None
    for mname in modulenames:
        exec ('import %s as mobj' % mname)
        suite.addTests(loader.loadTestsFromModule(mobj))
    return suite


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
