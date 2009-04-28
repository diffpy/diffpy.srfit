#!/usr/bin/env python
"""Tests for diffpy.srfit.structure package."""

import unittest

from diffpy.srfit.structure.wrappers import ParameterWrapper
from diffpy.srfit.fitbase.parameter import Parameter


class TestParameterWrapper(unittest.TestCase):

    def testWrapper(self):
        """Test the adapter.

        This adapts a Parameter to the Parameter interface. :)
        """
        l = Parameter("l", 3.14)

        # Try Accessor adaptation
        la = ParameterWrapper(l, "l", getter = Parameter.getValue, setter =
                Parameter.setValue)

        self.assertEqual(l.name, la.name)
        self.assertEqual(l.getValue(), la.getValue())

        # Change the parameter
        l.setValue(2.3)
        self.assertEqual(l.getValue(), la.getValue())

        # Change the adapter
        la.setValue(3.2)
        self.assertEqual(l.getValue(), la.getValue())

        # Try Attribute adaptation
        la = ParameterWrapper(l, "l", attr = "value")

        self.assertEqual(l.name, la.name)
        self.assertEqual("value", la.attr)
        self.assertEqual(l.getValue(), la.getValue())

        # Change the parameter
        l.setValue(2.3)
        self.assertEqual(l.getValue(), la.getValue())

        # Change the adapter
        la.setValue(3.2)
        self.assertEqual(l.getValue(), la.getValue())

        return

if __name__ == "__main__":

    unittest.main()

