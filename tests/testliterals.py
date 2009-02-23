#!/usr/bin/env python
"""Tests for refinableobj module."""

import diffpy.srfit.equation.literals as literals
import unittest


class TestArgument(unittest.TestCase):

    def testSetValue(self):
        """Test initialization."""
        l = literals.Argument()
        c1 = literals.Clicker()

        l.setValue(3.14)
        self.assertAlmostEqual(l.value, 3.14)
        self.assertTrue(c1 < l.clicker)
        return


class TestOperator(unittest.TestCase):

    def testAddLiteral(self):
        """Test adding a literal."""
        a1 = literals.Argument()
        op1 = literals.Operator()
        op1.nin = 1
        op2 = literals.Operator()

        op1.addLiteral(a1)
        self.assertEquals(op1.args, [a1])

        self.assertTrue(a1.clicker == op1.clicker)
        return
        



if __name__ == "__main__":

    unittest.main()

