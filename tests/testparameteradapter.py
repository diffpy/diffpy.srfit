#!/usr/bin/env python

import unittest
from pickle import dumps, loads

from diffpy.srfit.fit.parameters import Parameter
from diffpy.srfit.adapters.adaptersmod import ParameterAdapter
from diffpy.srfit.adapters.adaptersmod import funcgetter, nosetter
from diffpy.srfit.adapters.adaptersmod import UnboundOperator

class TestParameterAdapter(unittest.TestCase):

    def testAdapter(self):
        """Test the adapter."""
        l = Adapted()
        l.set(3.14)
        la = ParameterAdapter("l", l, getter = Adapted.get, 
                setter = Adapted.set)

        dl = dumps((l, la))
        lp, lap = loads(dl)

        self._testAdapter(l, la)
        self._testAdapter(lp, lap)
        return

    def _testAdapter(self, l, la):
        """Test the adapter."""

        self.assertEqual(l.get(), la.value)

        # Change the parameter. This will not register with the adapter because
        # it is not done through the adapter.
        l.set(2.3)
        self.assertNotEqual(l.get(), la.value)
        self.assertEqual(l.get(), la._get())

        # Change the adapter
        la.set(3.2)
        self.assertEqual(l.get(), la.value)

        # Constrain the adapter
        m = Parameter("m", 4)
        la.constrain(m)
        # The adapted object should be updated every time the adapter value is
        # called.
        self.assertEqual(m.value, la.value)
        self.assertEqual(l.get(), la.value)

        return

class Adapted(object):
    def get(self):
        return self._value
    def set(self, val):
        self._value = val
        return

if __name__ == "__main__":

    unittest.main()

