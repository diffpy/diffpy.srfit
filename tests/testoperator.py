#!/usr/bin/env python
import unittest

from diffpy.srfit.adapters import nodes

class TestOperator(unittest.TestCase):
    """Test Operator nodes."""

    def testOperator(self):
        """Test operator."""

        op1 = nodes.Operator("f")
        op1._setop(op)
        self.assertTrue(op is op1._getop())

        p1 = nodes.Parameter("a", 1.0)
        p2 = nodes.Parameter("b", 2.0)
        p1._addViewer(op1)
        p2._addViewer(op1)
        op1._args = [p1, p2]

        self.assertAlmostEqual(3, op1.value)
        self.assertAlmostEqual(3, op1._value)

        # Check lazy update.
        # FIXME - there has to be a better way to check for this.
        p1.value = 4.0
        self.assertTrue(None is op1._value)
        self.assertAlmostEqual(4.5, op1.value)
        self.assertAlmostEqual(4.5, op1._value)

        # Check pickle
        from pickle import dumps, loads
        dop1 = dumps(op1)
        op2 = loads(dop1)
        self.assertEqual(op2.value, op1.value)
        self.assertEqual(op2._args[0].name, p1.name)
        self.assertEqual(op2._args[0].value, p1.value)
        self.assertEqual(op2._args[1].name, p2.name)
        self.assertEqual(op2._args[1].value, p2.value)
        self.assertEqual(op2.name, op1.name)
        self.assertTrue(op2._getop() is op2._getop())
        return

    def testOperatorKW(self):
        """Test operator with key words."""

        op1 = nodes.Operator("f")
        op1._setop(op)
        self.assertTrue(op is op1._getop())

        p1 = nodes.Parameter("a", 1.0)
        p2 = nodes.Parameter("b", 2.0)
        p1._addViewer(op1)
        p2._addViewer(op1)
        op1._kw = { "a" : p1, "b" : p2 }

        self.assertAlmostEqual(3, op1.value)
        self.assertAlmostEqual(3, op1._value)

        p1.value = 4.0
        self.assertAlmostEqual(4.5, op1.value)

        return

def op(a,b):
    return a + b/a

if __name__ == "__main__":

    unittest.main()

