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

from diffpy.srfit.fitbase.constraint import Constraint
from diffpy.srfit.fitbase.recipeorganizer import equationFromString
from diffpy.srfit.fitbase.parameter import Parameter
from diffpy.srfit.equation.builder import EquationFactory


class TestConstraint(unittest.TestCase):

    def testConstraint(self):
        """Test the Constraint class."""

        p1 = Parameter("p1", 1)
        p2 = Parameter("p2", 2)

        factory = EquationFactory()

        factory.registerArgument("p1", p1)
        factory.registerArgument("p2", p2)

        c = Constraint()
        # Constrain p1 = 2*p2
        eq = equationFromString("2*p2", factory)
        c.constrain(p1, eq)

        self.assertTrue(p1.constrained)
        self.assertFalse(p2.constrained)

        eq2 = equationFromString("2*p2+1", factory)
        c2 = Constraint()
        self.assertRaises(ValueError, c2.constrain, p1, eq2)
        p2.setConst()
        eq3 = equationFromString("p1", factory)
        self.assertRaises(ValueError, c2.constrain, p2, eq3)

        p2.setValue(2.5)
        c.update()
        self.assertEqual(5.0, p1.getValue())

        p2.setValue(8.1)
        self.assertEqual(5.0, p1.getValue())
        c.update()
        self.assertEqual(16.2, p1.getValue())
        return


if __name__ == "__main__":
    unittest.main()
