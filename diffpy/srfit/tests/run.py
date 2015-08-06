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

"""Convenience module for executing all unit tests with

python -m diffpy.srfit.tests.run
"""
from __future__ import print_function
import six

if __name__ == '__main__':
    import sys
    from diffpy.srfit.tests import test
    # produce zero exit code for a successful test
    sys.exit(not test().wasSuccessful())

# End of file
