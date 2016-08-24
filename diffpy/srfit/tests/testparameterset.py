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

"""Tests for refinableobj module."""

import unittest

from diffpy.srfit.fitbase.parameter import Parameter
from diffpy.srfit.fitbase.parameterset import ParameterSet

class TestParameterSet(unittest.TestCase):

    def setUp(self):
        self.parset = ParameterSet("test")
        return

    def testAddParameterSet(self):
        """Test the addParameterSet method."""
        parset2 = ParameterSet("parset2")
        p1 = Parameter("parset2", 1)

        self.parset.addParameterSet(parset2)
        self.assertTrue(self.parset.parset2 is parset2)

        self.assertRaises(ValueError, self.parset.addParameterSet, p1)

        p1.name = "p1"
        parset2.addParameter(p1)

        self.assertTrue(self.parset.parset2.p1 is p1)

        return


if __name__ == "__main__":
    unittest.main()
