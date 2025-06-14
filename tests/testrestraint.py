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

from diffpy.srfit.fitbase.restraint import Restraint
from diffpy.srfit.fitbase.recipeorganizer import equationFromString
from diffpy.srfit.fitbase.parameter import Parameter
from diffpy.srfit.equation.builder import EquationFactory


class TestRestraint(unittest.TestCase):

    def testRestraint(self):
        """Test the Restraint class."""

        p1 = Parameter("p1", 1)
        p2 = Parameter("p2", 2)

        factory = EquationFactory()

        factory.registerArgument("p1", p1)
        factory.registerArgument("p2", p2)

        # Restrain 1 <  p1 + p2 < 5
        eq = equationFromString("p1 + p2", factory)
        r = Restraint(eq, 1, 5)

        # This should have no penalty
        p1.setValue(1)
        p2.setValue(1)
        self.assertEqual(0, r.penalty())

        # Make p1 + p2 = 0
        # This should have a penalty of 1*(1 - 0)**2 = 1
        p1.setValue(-1)
        p2.setValue(1)
        self.assertEqual(1, r.penalty())

        # Make p1 + p2 = 8
        # This should have a penalty of 1*(8 - 5)**2 = 9
        p1.setValue(4)
        p2.setValue(4)
        self.assertEqual(9, r.penalty())

        # Set the chi^2 to get a dynamic penalty
        r.scaled = True
        self.assertEqual(13.5, r.penalty(1.5))

        # Make a really large number to check the upper bound
        import numpy
        r.ub = numpy.inf
        p1.setValue(1e100)
        self.assertEqual(0, r.penalty())

        return


if __name__ == "__main__":
    unittest.main()
