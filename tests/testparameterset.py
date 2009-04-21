#!/usr/bin/env python
"""Tests for refinableobj module."""

import unittest

from diffpy.srfit.fitbase.parameter import Parameter, ParameterSet

class TestParameterSet(unittest.TestCase):

    def setUp(self):
        self.parset = ParameterSet("test")
        return

    def testInit(self):
        """Test the initialization."""
        self.assertEquals("test", self.parset.name)
        self.assertTrue(hasattr(self.parset, "clicker"))
        self.assertTrue(hasattr(self.parset, "constraints"))
        self.assertTrue(hasattr(self.parset, "restraints"))
        self.assertTrue(hasattr(self.parset, "pardict"))
        self.assertTrue(hasattr(self.parset, "suborganizers"))
        self.assertTrue(hasattr(self.parset, "_eqfactory"))
        return

    def testInsert(self):
        """Test the insert method."""

        p1 = Parameter("p1", 1)
        p2 = Parameter("p1", 2)
        p3 = Parameter("", 3)

        # Check normal insert
        self.parset.insert(p1)
        self.assertTrue("p1" in self.parset.pardict.keys())
        self.assertTrue(p1 in self.parset.pardict.values())

        # Try to insert another parameter with the same name
        self.assertRaises(ValueError, self.parset.insert, p2)

        # Try to insert a parameter without a name
        self.assertRaises(ValueError, self.parset.insert, p3)

        return

    def testClickers(self):
        """Test make sure that parameters are observed by parameter sets."""
        p1 = Parameter("p1", 1)
        self.parset.insert(p1)

        self.assertTrue(self.parset.clicker >= p1.clicker)

        p1.setValue(1.234)
        self.assertTrue(self.parset.clicker >= p1.clicker)

        return

    def testAccess(self):
        """Check attribute access."""
        p1 = Parameter("p1", 1)
        self.parset.insert(p1)

        self.assertTrue(p1 is self.parset.p1)

        self.assertRaises(AttributeError, self.parset.__getattr__, "p2")

        return




if __name__ == "__main__":

    unittest.main()

