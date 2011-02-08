#!/usr/bin/env python
import unittest

from diffpy.srfit.adapters import nodes
from diffpy.srfit.adapters.adaptersmod import MethodAdapter
from diffpy.srfit.adapters.adaptersmod import UnboundOperator

class TestMethodAdapter(unittest.TestCase):
    """Test MethodAdapter nodes."""

    def testMethodAdapter(self):
        """Test unbound operator."""

        uop = MethodAdapter("f", op, getter, setter)
        self.assertTrue(uop._getop() is op)
        self.assertTrue(uop._setop(3) is 3)

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

        # Test pickling
        from pickle import dumps, loads
        dop3 = dumps(op3)
        op4 = loads(dop3)
        self.assertEqual(op3.name, op4.name)
        self.assertEqual(op3.value, op4.value)

        return

def getter(op):
    return op
def setter(op, val):
    return val
def op(a,b):
    return a + b/a

if __name__ == "__main__":

    unittest.main()

