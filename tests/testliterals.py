#!/usr/bin/env python
"""Tests for refinableobj module."""

import unittest

from diffpy.srfit.util.clicker import Clicker
import diffpy.srfit.equation as equation
import diffpy.srfit.equation.literals as literals


class TestArgument(unittest.TestCase):

    def testSetValue(self):
        """Test initialization."""
        l = literals.Argument()
        c1 = Clicker()

        l.setValue(3.14)
        self.assertAlmostEqual(l.getValue(), 3.14)
        self.assertTrue(c1 < l.clicker)

        # Try again
        c1.click()
        l.setValue(3.14)
        self.assertTrue(c1 > l.clicker)

        # Try array
        import numpy
        x = numpy.arange(0, 10, 0.1)
        l.setValue(x)
        self.assertTrue( l.getValue() is x )

        # Change the array
        y = numpy.arange(0, 10, 0.5)
        l.setValue(y)
        self.assertTrue( l.getValue() is y )

        # Back to scalar
        l.setValue(1.01)
        self.assertAlmostEqual(l.getValue(), 1.01)
        return


class TestOperator(unittest.TestCase):

    def testAddLiteral(self):
        """Test adding a literal."""
        a1 = literals.Argument()
        op1 = literals.Operator()
        op2 = literals.Operator()

        op1.addLiteral(a1)
        self.assertEquals(op1.args, [a1])

        self.assertTrue(a1.clicker == op1.clicker)
        return
        

class TestPartition(unittest.TestCase):

    def testAddArgument(self):
        """Test adding an argument."""
        a1 = literals.Argument(value = 1.0)
        a2 = literals.Argument(value = 2.0)

        p1 = literals.Partition()
        p1.addArgument(a1, "tag1")
        self.assertTrue(a1.clicker, p1.clicker)
        p1.addArgument(a2, "tag1", "tag2")
        self.assertTrue(a2.clicker, p1.clicker)

        self.assertEquals([a1, a2], p1.args)
        self.assertEquals(set(["tag1", "tag2"]), p1.tags)
        self.assertEquals([0, 1], p1.tagmap["tag1"])
        self.assertEquals([1], p1.tagmap["tag2"])

        p1._prepare()
        self.assertEquals( 3, p1.combine(p1._partvals) )
        return





if __name__ == "__main__":

    unittest.main()

