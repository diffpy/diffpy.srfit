#!/usr/bin/env python
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

        p2.setValue(2.5)
        c.update()
        self.assertEquals(5.0, p1.getValue())

        p2.setValue(8.1)
        self.assertEquals(5.0, p1.getValue())
        c.update()
        self.assertEquals(16.2, p1.getValue())
        return

if __name__ == "__main__":

    unittest.main()

