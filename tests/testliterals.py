#!/usr/bin/env python
"""Tests for refinableobj module."""

import diffpy.srfit.equation.literals as literals
import diffpy.srfit.equation as equation
import unittest


class TestArgument(unittest.TestCase):

    def testSetValue(self):
        """Test initialization."""
        l = literals.Argument()
        c1 = equation.Clicker()

        l.setValue(3.14)
        self.assertAlmostEqual(l.value, 3.14)
        self.assertTrue(c1 < l.clicker)
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

