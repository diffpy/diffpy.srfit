#!/usr/bin/env python
"""Tests for refinableobj module."""

import unittest

from diffpy.srfit.fitbase.restraint import Restraint, BoundsRestraint
from diffpy.srfit.fitbase.modelorganizer import equationFromString
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

        r = Restraint()
        # Restrain 1 <  p1 + p2 < 5
        eq = equationFromString("p1 + p2", factory)
        r.restrain(eq, 1, 5)

        # This should have no penalty
        p1.setValue(1)
        p2.setValue(1)
        self.assertEquals(0, r.penalty())

        # Make p1 + p2 = 0
        # This should have a penalty of 1*(1 - 0)**2 = 1
        p1.setValue(-1)
        p2.setValue(1)
        self.assertEquals(1, r.penalty())

        # Make p1 + p2 = 8
        # This should have a penalty of 1*(8 - 5)**2 = 9
        p1.setValue(4)
        p2.setValue(4)
        self.assertEquals(9, r.penalty())

        # Try a different penalty function
        from numpy import inf
        eq = equationFromString("p1", factory)
        r.restrain(eq, 0, inf, 2, 4, True)

        # Make p1 = -2
        # This should give a penalty of 2*2**4 = 32
        p1.setValue(-2)
        self.assertEquals(32, r.penalty())

        # Set the chi^2 to get a dynamic penalty
        self.assertEquals(48, r.penalty(1.5))

        # Make a really large number to check the upper bound
        p1.setValue(1e100)
        self.assertEquals(0, r.penalty())

        return

    def testBoundsRestraint(self):
        """Test the BoundsRestraint class."""

        p1 = Parameter("p1", 1)
        p2 = Parameter("p2", 2)

        factory = EquationFactory()

        factory.registerArgument("p1", p1)
        factory.registerArgument("p2", p2)

        r = BoundsRestraint()
        # Restrain 1 <  p1 + p2 < 5
        eq = equationFromString("p1 + p2", factory)
        r.confine(eq, 1, 5)

        # This should have no penalty
        p1.setValue(1)
        p2.setValue(1)
        self.assertEquals(0, r.penalty())

        # This should have an infinite penalty
        p1.setValue(9)
        from numpy import inf
        self.assertEquals(inf, r.penalty())
        p1.setValue(-2)
        self.assertEquals(inf, r.penalty())

        return

if __name__ == "__main__":

    unittest.main()

