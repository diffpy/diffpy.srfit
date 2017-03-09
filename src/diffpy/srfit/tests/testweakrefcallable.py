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

"""\
Unit tests for the weakrefcallable module.
"""


import unittest
import pickle

from diffpy.srfit.fitbase import FitContribution
from diffpy.srfit.util.weakrefcallable import weak_ref

# ----------------------------------------------------------------------------

class TestWeakBoundMethod(unittest.TestCase):

    def setUp(self):
        self.f = FitContribution('f')
        self.holder = set()
        self.w = weak_ref(self.f._flush, holder=self.holder)
        self.holder.add(self.w)
        return


    def tearDown(self):
        self.assertTrue(self.w in self.holder)
        self.f = None
        self.assertTrue(None is self.w('any', 'argument'))
        self.assertFalse(self.w in self.holder)
        return


    def test___init__(self):
        """check WeakBoundMethod.__init__()
        """
        # make sure there is no 'weak object has gone away' warning.
        wf = weak_ref(self.f._flush, holder=self.holder)
        del wf
        return


    def test___call__(self):
        """check WeakBoundMethod.__call__()
        """
        f = self.f
        f.newParameter('x', 5)
        f.setEquation('3 * x')
        wx = weak_ref(f._eq._flush, holder=self.holder)
        self.holder.add(wx)
        self.assertEqual(15, f.evaluate())
        self.assertEqual(15, f._eq._value)
        wx(())
        self.assertTrue(None is f._eq._value)
        self.assertEqual(15, f.evaluate())
        self.assertTrue(wx in self.holder)
        f.setEquation('2 * t')
        self.assertFalse(wx in self.holder)
        return


    def test___hash__(self):
        """check WeakBoundMethod.__hash__()
        """
        f1 = FitContribution('f1')
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
        """check WeakBoundMethod.__eq__()
        """
        f1 = FitContribution('f1')
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

# End of class TestWeakBoundMethod

if __name__ == '__main__':
    unittest.main()

# End of file
