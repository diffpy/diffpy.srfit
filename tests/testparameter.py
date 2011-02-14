#!/usr/bin/env python

import unittest

from pickle import dumps, loads

from diffpy.srfit.fit.parameters import *
from diffpy.srfit.util import messages

class TestParameter(unittest.TestCase):

    def testPar(self):
        """Test the Par factory."""
        l = Par("l", 3.14)

        # pickle
        dl = dumps(l)
        lp = loads(dl)

        self._testPar(l)
        self._testPar(lp)
        return

    def testPars(self):
        """Test the Pars factory."""
        l, m = Pars(("l", 3.14), ("m", 3.14))
        self.assertEqual("l", l.name)
        self.assertEqual("m", m.name)

        # pickle
        dl = dumps(l)
        lp = loads(dl)
        self._testPar(l)
        self._testPar(lp)
        self._testPar(m)
        return

    def _testPar(self, l):
        """Test the l factory."""
        self.assertFalse(l.isVaried())
        l.vary()
        self.assertTrue(l.isVaried())
        self.assertTrue( isinstance(l, Parameter) )
        l.fix()
        self.assertFalse(l.isVaried())
        return

    def testVar(self):
        """Test the Var factory."""
        l = Var("l", 3.14)

        # pickle
        dl = dumps(l)
        lp = loads(dl)

        self._testVar(l)
        self._testVar(lp)
        return

    def testPars(self):
        """Test the Vars factory."""
        l, m = Vars(("l", 3.14), ("m", 3.14))
        self.assertEqual("l", l.name)
        self.assertEqual("m", m.name)

        # pickle
        dl = dumps(l)
        lp = loads(dl)
        self._testVar(l)
        self._testVar(lp)
        self._testVar(m)
        return

    def _testVar(self, l):
        """Test the Var factory."""
        self.assertTrue(l.isVaried())
        l.vary()
        self.assertTrue(l.isVaried())
        self.assertTrue( isinstance(l, Parameter) )
        l.fix()
        self.assertFalse(l.isVaried())
        l.vary()
        return

    def testFixed(self):
        """Test the Fixed factory."""
        l = Fixed("l", 3.14)

        # pickle
        dl = dumps(l)
        lp = loads(dl)

        self._testFixed(l)
        self._testFixed(lp)
        return

    def _testFixed(self, l):
        """Test the Fixed factory."""
        self.assertFalse(l.isVaried())
        self.assertRaises(AttributeError, l.vary)
        self.assertFalse(l.isVaried())
        l.fix()
        self.assertFalse(l.isVaried())
        self.assertTrue( isinstance(l, Parameter) )
        return

    def testValue(self):
        """Test value setting and getting."""
        l = Parameter("l")

        dl = dumps(l)
        lp = loads(dl)

        self._testValue(l)
        self._testValue(lp)
        return

    def _testValue(self, l):
        """Test value setting and getting."""

        l.set(3.14)
        self.assertAlmostEqual(3.14, l.get())

        # Try array
        import numpy
        x = numpy.arange(0, 10, 0.1)
        l.set(x)
        self.assertTrue( l.get() is x )
        self.assertTrue( l.value is x )

        # Change the array
        y = numpy.arange(0, 10, 0.5)
        l.value = y
        self.assertTrue( l.get() is y )

        # Back to scalar
        l.set(1.01)
        self.assertAlmostEqual(1.01, l.value)
        return

    def testConstrain(self):
        """Test constraints."""
        l = Parameter("l", 3)
        m = Parameter("m", 4)
        self.assertFalse( l.isConstrained() )
        self.assertFalse( m.isConstrained() )
        self.assertAlmostEqual(3, l.value)
        l.constrain(m)

        # Pickle
        dl = dumps((l, m))
        lp, mp = loads(dl)

        self._testConstrain(l, m)
        self._testConstrain(lp, mp)
        return

    def _testConstrain(self, l, m):
        """Test constraints."""
        self.assertTrue( l.isConstrained() )
        self.assertTrue(l._constraint is m)
        self.assertEqual(4, l.value)

        m.value = 5
        self.assertAlmostEqual(5, m.value)
        self.assertAlmostEqual(5, l.value)

        self.assertRaises(ValueError, l.set, 1)

        l.unconstrain()
        self.assertFalse( l.isConstrained() )
        self.assertAlmostEqual(5, l.value)

        l.value = 9
        self.assertAlmostEqual(9, l.value)

        # check for circular constraint
        p = Parameter("p", 5)
        p.constrain(l)
        self.assertRaises(ValueError, l.constrain, p * m)
        pass

    def testVary(self):
        """Test parameter varying."""
        l = Parameter("l")

        # Pickle
        dl = dumps(l)
        lp = loads(dl)

        self._testVary(l)
        self._testVary(lp)
        return

    def _testVary(self, l):
        """Test parameter varying."""
        self.assertFalse(l.isVaried())

        l.vary()
        self.assertTrue(l.isVaried())

        l.vary()
        self.assertTrue(l.isVaried())

        l.vary(3)
        self.assertAlmostEqual(3, l.value)

        l.fix()
        self.assertFalse(l.isVaried())

        l.fix(4)
        self.assertAlmostEqual(4, l.value)
        return


    def _test_arithmetic(self):
        """Test operations on parameters.

        Behavior:
        Arithmetic operations (+, -, *, /, **, -(unary), +(unary), abs) return
        a node that evaluates the operation. The fact that the operands need
        not be nodes is not part of the public interface, so is not tested
        here. Passing a numpy array to an operand will return an array of nodes
        if the array is on the left of the operation. This is not desired
        behavior, but it is unavoidable, therefore we only operate on nodes.

        We only test the node type, not the operations.
        
        """
        n1 = Parameter("n1", 1)
        n2 = Parameter("n2", 2)

        out = n1 + n2
        self.assertTrue( functions.add is out._obj )
        self.assertTrue( n1 in out._args )
        self.assertTrue( n2 in out._args )
        self.assertTrue( {} == out._kw )
        self.assertTrue( out in n1._viewers )
        self.assertTrue( out in n2._viewers )
        self.assertEqual(3, out.value)

        out = n1 - n2
        self.assertTrue( functions.subtract is out._obj )
        self.assertTrue( n1 in out._args )
        self.assertTrue( n2 in out._args )
        self.assertTrue( {} == out._kw )
        self.assertTrue( out in n1._viewers )
        self.assertTrue( out in n2._viewers )
        self.assertEqual(-1, out.value)

        out = n1 * n2
        self.assertTrue( functions.multiply is out._obj )
        self.assertTrue( n1 in out._args )
        self.assertTrue( n2 in out._args )
        self.assertTrue( {} == out._kw )
        self.assertTrue( out in n1._viewers )
        self.assertTrue( out in n2._viewers )
        self.assertEqual(2, out.value)

        out = n1 / n2
        self.assertTrue( functions.divide is out._obj )
        self.assertTrue( n1 in out._args )
        self.assertTrue( n2 in out._args )
        self.assertTrue( {} == out._kw )
        self.assertTrue( out in n1._viewers )
        self.assertTrue( out in n2._viewers )
        self.assertEqual(0.5, out.value)

        out = n1 ** n2
        self.assertTrue( functions.power is out._obj )
        self.assertTrue( n1 in out._args )
        self.assertTrue( n2 in out._args )
        self.assertTrue( {} == out._kw )
        self.assertTrue( out in n1._viewers )
        self.assertTrue( out in n2._viewers )
        self.assertEqual(1, out.value)

        out = n1 ** 2
        self.assertTrue( functions.power is out._obj )
        self.assertTrue( n1 in out._args )
        self.assertTrue( {} == out._kw )
        self.assertTrue( out in n1._viewers )
        self.assertEqual(1, out.value)

        out = -n1
        self.assertTrue( functions.negative is out._obj )
        self.assertTrue( n1 in out._args )
        self.assertTrue( {} == out._kw )
        self.assertTrue( out in n1._viewers )
        self.assertEqual(-1, out.value)

        out = +n1
        self.assertTrue( out is n1 )

        out = abs(n1)
        self.assertTrue( functions.abs is out._obj )
        self.assertTrue( n1 in out._args )
        self.assertTrue( out in n1._viewers )
        self.assertEqual(1, out.value)

        # In-place operations transfer the name, but that is all
        out = n1 
        out += n2
        self.assertTrue( functions.add is out._obj )
        self.assertTrue( n1 in out._args )
        self.assertTrue( n2 in out._args )
        self.assertTrue( {} == out._kw )
        self.assertTrue( out in n1._viewers )
        self.assertTrue( out in n2._viewers )
        self.assertTrue( out.name is n1.name )
        self.assertEqual(3, out.value)

        out = n1
        out -= n2
        self.assertTrue( functions.subtract is out._obj )
        self.assertTrue( n1 in out._args )
        self.assertTrue( n2 in out._args )
        self.assertTrue( {} == out._kw )
        self.assertTrue( out in n1._viewers )
        self.assertTrue( out in n2._viewers )
        self.assertTrue( out.name is n1.name )
        self.assertEqual(-1, out.value)

        out = n1
        out *= n2
        self.assertTrue( functions.multiply is out._obj )
        self.assertTrue( n1 in out._args )
        self.assertTrue( n2 in out._args )
        self.assertTrue( {} == out._kw )
        self.assertTrue( out in n1._viewers )
        self.assertTrue( out in n2._viewers )
        self.assertTrue( out.name is n1.name )
        self.assertEqual(2, out.value)

        out = n1
        out /= n2
        self.assertTrue( functions.divide is out._obj )
        self.assertTrue( n1 in out._args )
        self.assertTrue( n2 in out._args )
        self.assertTrue( {} == out._kw )
        self.assertTrue( out in n1._viewers )
        self.assertTrue( out in n2._viewers )
        self.assertTrue( out.name is n1.name )
        self.assertEqual(0.5, out.value)

        out = n1
        out **= n2
        self.assertTrue( functions.power is out._obj )
        self.assertTrue( n1 in out._args )
        self.assertTrue( n2 in out._args )
        self.assertTrue( {} == out._kw )
        self.assertTrue( out in n1._viewers )
        self.assertTrue( out in n2._viewers )
        self.assertTrue( out.name is n1.name )
        self.assertEqual(1, out.value)

        return

if __name__ == "__main__":

    unittest.main()

