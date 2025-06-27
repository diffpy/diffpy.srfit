#!/usr/bin/env python
##############################################################################
#
# diffpy.srfit      Complex Modeling Initiative
#                   (c) 2016 Brookhaven Science Associates,
#                   Brookhaven National Laboratory.
#                   All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
##############################################################################
"""Unit tests for the weakrefcallable module."""


import pickle
import unittest

from diffpy.srfit.fitbase import FitContribution
from diffpy.srfit.fitbase.parameter import Parameter
from diffpy.srfit.util.weakrefcallable import WeakBoundMethod, weak_ref

# ----------------------------------------------------------------------------


class TestWeakBoundMethod(unittest.TestCase):

    def setUp(self):
        self.f = FitContribution("f")
        self.f.setEquation("7")
        self.w = weak_ref(self.f._eq._flush, fallback=_fallback_example)
        return

    def tearDown(self):
        self.f = None
        self.assertTrue(None is self.w._wref())
        obj, args, kw = self.w("any", "argument", foo=37)
        self.assertTrue(obj is self.w)
        self.assertEqual(("any", "argument"), args)
        self.assertEqual({"foo": 37}, kw)
        return

    def test___init__(self):
        """Check WeakBoundMethod.__init__()"""
        self.assertTrue(self.w.fallback is _fallback_example)
        wf = weak_ref(self.f._flush)
        self.assertTrue(None is wf.fallback)
        return

    def test___call__(self):
        """Check WeakBoundMethod.__call__()"""
        f = self.f
        self.assertEqual(7, f.evaluate())
        self.assertEqual(7, f._eq._value)
        # verify f has the same effect as f._eq._flush
        self.w(())
        self.assertTrue(None is f._eq._value)
        # check WeakBoundMethod behavior with no fallback
        x = Parameter("x", value=3)
        wgetx = weak_ref(x.getValue)
        self.assertEqual(3, wgetx())
        del x
        self.assertRaises(ReferenceError, wgetx)
        return

    def test___hash__(self):
        """Check WeakBoundMethod.__hash__()"""
        f1 = FitContribution("f1")
        w1 = weak_ref(f1._flush)
        h0 = hash(w1)
        del f1
        self.assertTrue(None is w1._wref())
        self.assertEqual(h0, hash(w1))
        w1c1 = pickle.loads(pickle.dumps(w1))
        w1c2 = pickle.loads(pickle.dumps(w1))
        self.assertEqual(hash(w1c1), hash(w1c2))
        return

    def test___eq__(self):
        """Check WeakBoundMethod.__eq__()"""
        f1 = FitContribution("f1")
        w1 = weak_ref(f1._flush)
        w2 = weak_ref(f1._flush)
        self.assertEqual(w1, w2)
        w1c = pickle.loads(pickle.dumps(w1))
        # pickle-copied objects should have empty reference
        self.assertTrue(None is w1c._wref())
        self.assertNotEqual(w1, w1c)
        del f1
        self.assertTrue(None is w1._wref())
        self.assertEqual(w1, w1c)
        w1cc = pickle.loads(pickle.dumps(w1c))
        self.assertTrue(None is w1cc._wref())
        self.assertEqual(w1c, w1cc)
        self.assertEqual(w1, w1cc)
        return

    def test_pickling(self):
        """Verify unpickling works when it involves __hash__ call."""
        holder = set([self.w])
        objs = [holder, self.f._eq, self.w]
        data = pickle.dumps(objs)
        objs2 = pickle.loads(data)
        h2, feq2, w2 = objs2
        self.assertTrue(w2 in h2)
        self.assertTrue(feq2 is w2._wref())
        return

    def test_observable_deregistration(self):
        """Check if Observable drops dead Observer."""
        f = self.f
        x = f.newParameter("x", 5)
        f.setEquation("3 * x")
        self.assertEqual(15, f.evaluate())
        self.assertEqual(15, f._eq._value)
        # get one of the observer callables that are associated with f
        xof = next(iter(x._observers))
        self.assertTrue(isinstance(xof, WeakBoundMethod))
        # changing value of x should reset f._eq
        x.setValue(x.value + 1)
        self.assertTrue(None is f._eq._value)
        self.assertEqual(18, f.evaluate())
        # deallocate f now
        self.f = f = None
        self.assertTrue(xof in x._observers)
        # since f does not exist anymore, the next notification call
        # should drop the associated observer.
        x.setValue(x.value + 1)
        self.assertEqual(0, len(x._observers))
        return


# End of class TestWeakBoundMethod

# Local Routines -------------------------------------------------------------


def _fallback_example(wbm, *args, **kwargs):
    return (wbm, args, kwargs)


if __name__ == "__main__":
    unittest.main()

# End of file
