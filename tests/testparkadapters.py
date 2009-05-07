#!/usr/bin/env python
"""Tests for diffpy.srfit.park.adapters module."""

import unittest

from diffpy.srfit.fitbase.parameter import Parameter
from diffpy.srfit.park.adapters import ParkParameterProxy


class TestParameterProxy(unittest.TestCase):

    def testParameterProxy(self):
        """Test ParamterProxy methods."""
        par = Parameter("test", 3.14)
        proxy = ParkParameterProxy(par)

        self.assertEquals("fitted", proxy.status)

        self.assertEquals(par.name, proxy.name)
        self.assertEquals(par.getValue(), proxy.value)

        proxy.value = 2.0
        self.assertEquals(par.getValue(), proxy.value)

        par.setValue(8.1)
        self.assertEquals(par.getValue(), proxy.value)
        return

if __name__ == "__main__":

    unittest.main()

