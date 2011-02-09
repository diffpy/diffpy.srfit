#!/usr/bin/env python
import unittest

from diffpy.srfit.adapters import nodes
from diffpy.srfit.adapters.adaptersmod import UnboundOperator

class TestUnboundOperator(unittest.TestCase):
    """Test UnboundOperator nodes."""

    def testUnboundOperator(self):
        """Test unbound operator."""

        uop = UnboundOperator("f", op)

        p1 = nodes.Parameter("a", 1.0)
        p2 = nodes.Parameter("b", 2.0)

        op1 = uop(p1, p2)
        self.assertTrue(uop is op1._obj)
        self.assertTrue( p1 in op1._args )
        self.assertTrue( p2 in op1._args )
        self.assertTrue( {} == op1._kw )
        self.assertTrue( op1 in p1._viewers )
        self.assertTrue( op1 in p2._viewers )

        self.assertAlmostEqual(3, op1.value)

        # try this again
        op2 = uop(p1, p2)
        self.assertTrue( op1 is op2 )

        # Verify we can pass non-adapted objects
        op3 = uop(1.0, 4.0)
        self.assertAlmostEqual(5, op3.value)

        # Try passing too much or too little or wrong key-words
        self.assertRaises(TypeError, uop, p1)
        self.assertRaises(TypeError, uop, p1, p2, p2)
        self.assertRaises(TypeError, uop, p1, p2, a = p2)
        self.assertRaises(TypeError, uop, p1, p2, a = p1)
        self.assertRaises(TypeError, uop, p1, p2, c = p1)
        self.assertRaises(TypeError, uop, a = p1)

        self._testPickle(uop)
        return

    def testUnboundOperatorKW(self):
        """Test unbound operator with keywords."""

        uop = UnboundOperator("f", op)

        p1 = nodes.Parameter("a", 1.0)
        p2 = nodes.Parameter("b", 2.0)

        op1 = uop(p1, b = p2)
        self.assertTrue(uop is op1._obj)
        self.assertTrue( p1 in op1._args )
        self.assertTrue( p2 in op1._kw.values() )
        self.assertTrue( op1 in p1._viewers )
        self.assertTrue( op1 in p2._viewers )

        self.assertAlmostEqual(3, op1.value)

        # Try this again. We should be able to get the same object even if we
        # don't use the keywords.
        op2 = uop(p1, p2)
        self.assertTrue( op1 is op2 )

        # Verify we can pass non-adapted objects
        op3 = uop(1.0, 4.0)
        self.assertAlmostEqual(5, op3.value)

        self._testPickle(uop)
        return

    def testMethod(self):
        """Test unbound operator on a method."""

        uop = UnboundOperator("f", mtest.f)

        p1 = nodes.Parameter("a", 1.0)
        p2 = nodes.Parameter("b", 2.0)

        op1 = uop(p1, p2)
        self.assertTrue(uop is op1._obj)
        self.assertTrue( p1 in op1._args )
        self.assertTrue( p2 in op1._args )
        self.assertTrue( {} == op1._kw )
        self.assertTrue( op1 in p1._viewers )
        self.assertTrue( op1 in p2._viewers )

        self.assertAlmostEqual(3, op1.value)

        # try this again
        op2 = uop(p1, p2)
        self.assertTrue( op1 is op2 )

        # Verify we can pass non-adapted objects
        op3 = uop(1.0, 4.0)
        self.assertAlmostEqual(5, op3.value)

        self.assertRaises(TypeError, uop, p1)
        self.assertRaises(TypeError, uop, p1, p2, p2)
        self.assertRaises(TypeError, uop, p1, p2, a = p2)
        self.assertRaises(TypeError, uop, p1, p2, a = p1)
        self.assertRaises(TypeError, uop, p1, p2, c = p1)
        self.assertRaises(TypeError, uop, a = p1)

        self._testPickle(uop)
        return

    def testBuiltin(self):
        """Test unbound operator on a built-in function."""

        uop = UnboundOperator("f", max)

        p1 = nodes.Parameter("a", 1.0)
        p2 = nodes.Parameter("b", 2.0)

        op1 = uop(p1, p2)
        self.assertTrue(uop is op1._obj)
        self.assertTrue( p1 in op1._args )
        self.assertTrue( p2 in op1._args )
        self.assertTrue( {} == op1._kw )
        self.assertTrue( op1 in p1._viewers )
        self.assertTrue( op1 in p2._viewers )

        self.assertAlmostEqual(2, op1.value)

        # try this again
        op2 = uop(p1, p2)
        self.assertTrue( op1 is op2 )

        # Verify we can pass non-adapted objects
        op3 = uop(1.0, 4.0)
        self.assertAlmostEqual(4, op3.value)

        # Skip these. We know they fail for built-in functions. The only way to
        # check these is to evaluate the function. However, the arguments may
        # not be ready for evaluation when we create a function, so we really
        # can't do anything.
        return
        self.assertRaises(TypeError, uop, p1)
        self.assertRaises(TypeError, uop, p1, p2, p2)
        self.assertRaises(TypeError, uop, p1, p2, a = p2)
        self.assertRaises(TypeError, uop, p1, p2, a = p1)
        self.assertRaises(TypeError, uop, p1, p2, c = p1)
        self.assertRaises(TypeError, uop, a = p1)
        return

        self._testPickle(uop)
        return

    def testFunctor(self):
        """Test unbound operator on a built-in function."""

        uop = UnboundOperator("f", ftest)

        p1 = nodes.Parameter("a", 1.0)
        p2 = nodes.Parameter("b", 2.0)

        op1 = uop(p1, p2)
        self.assertTrue(uop is op1._obj)
        self.assertTrue( p1 in op1._args )
        self.assertTrue( p2 in op1._args )
        self.assertTrue( {} == op1._kw )
        self.assertTrue( op1 in p1._viewers )
        self.assertTrue( op1 in p2._viewers )

        self.assertAlmostEqual(3, op1.value)

        # Check lazy update.
        # FIXME - there has to be a better way to check for this.
        p1.value = 4.0
        self.assertTrue(None is op1._value)
        self.assertAlmostEqual(4.5, op1.value)
        self.assertAlmostEqual(4.5, op1._value)

        # try this again
        op2 = uop(p1, p2)
        self.assertTrue( op1 is op2 )

        # Verify we can pass non-adapted objects
        op3 = uop(1.0, 4.0)
        self.assertAlmostEqual(5, op3.value)

        self.assertRaises(TypeError, uop, p1)
        self.assertRaises(TypeError, uop, p1, p2, p2)
        self.assertRaises(TypeError, uop, p1, p2, a = p2)
        self.assertRaises(TypeError, uop, p1, p2, a = p1)
        self.assertRaises(TypeError, uop, p1, p2, c = p1)
        self.assertRaises(TypeError, uop, a = p1)

        self._testPickle(uop)
        return

    def testUFunc(self):
        """Test unbound operator on a built-in function."""

        import numpy
        uop = UnboundOperator("f", numpy.add)

        p1 = nodes.Parameter("a", 1.0)
        p2 = nodes.Parameter("b", 2.0)

        op1 = uop(p1, p2)
        self.assertTrue(uop is op1._obj)
        self.assertTrue( p1 in op1._args )
        self.assertTrue( p2 in op1._args )
        self.assertTrue( {} == op1._kw )
        self.assertTrue( op1 in p1._viewers )
        self.assertTrue( op1 in p2._viewers )

        self.assertAlmostEqual(3, op1.value)

        # try this again
        op2 = uop(p1, p2)
        self.assertTrue( op1 is op2 )

        # Verify we can pass non-adapted objects
        op3 = uop(1.0, 4.0)
        self.assertAlmostEqual(5, op3.value)

        self.assertRaises(TypeError, uop, p1)
        self.assertRaises(TypeError, uop, p1, p2, p2)
        self.assertRaises(TypeError, uop, p1, p2, a = p2)
        self.assertRaises(TypeError, uop, p1, p2, a = p1)
        self.assertRaises(TypeError, uop, p1, p2, c = p1)
        self.assertRaises(TypeError, uop, a = p1)

        self._testPickle(uop)
        return

    def _testPickle(self, uop):
        """Test pickling the various operations."""
        from pickle import dumps, loads
        duop = dumps(uop)
        uop2 = loads(duop)
        self.assertEqual(uop2.name, uop.name)
        self.assertTrue(uop2._getop() is uop2._getop())
        return

def op(a, b):
    return a + b/a

class _Functor(object):
    def __call__(self, a, b): return a + b/a
ftest = _Functor()

class _Method(object):
    def f(self, a, b): return a + b/a
mtest = _Method()

if __name__ == "__main__":

    unittest.main()

