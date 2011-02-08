#!/usr/bin/env python

import unittest

from utils import TestViewer
from diffpy.srfit.adapters import nodes
from diffpy.srfit.fit import functions

class TestNode(unittest.TestCase):
    """Test methods common to all nodes."""

    def testName(self):
        """Test aspects of naming.

        Behavior:
        Nodes have a name. It is accessible through the 'name' attribute.
        Alternately, the 'rename' method assigns a new name and returns the
        node for chaining.

        """
        n1 = nodes.Node("n1")
        self.assertEquals(n1.name, "n1")
        retval = n1.rename("test")
        self.assertTrue(n1 is retval)
        self.assertEquals(n1.name, "test")
        return

    def test_notify(self):
        """Test notification.
        
        Behavior:
        The _notify method asks all viewers to respond to a message.
        
        """

        n1 = nodes.Node("n1")
        n2 = nodes.Node("n2")
        v1 = TestViewer()
        v2 = TestViewer()

        n1._addViewer(n2)
        n1._addViewer(v1)
        n2._addViewer(v2)

        self.assertTrue(n2 in n1._viewers)
        self.assertTrue(v1 in n1._viewers)
        self.assertTrue(v2 in n2._viewers)

        msg1 = "test1"
        n1._notify(msg1)
        self.assertTrue( msg1 is v1.msg )
        self.assertTrue( msg1 is v2.msg )

        msg2 = "test2"
        n2._notify(msg2)
        self.assertTrue( msg1 is v1.msg )
        self.assertTrue( msg2 is v2.msg )
        return

    def test_arithmetic(self):
        """Test operations on a node.

        Behavior:
        Arithmetic operations (+, -, *, /, **, -(unary), +(unary), abs) return
        a node that evaluates the operation. The fact that the operands need
        not be nodes is not part of the public interface, so is not tested
        here. Passing a numpy array to an operand will return an array of nodes
        if the array is on the left of the operation. This is not desired
        behavior, but it is unavoidable, therefore we only operate on nodes.

        We only test the node type, not the operations.
        
        """
        n1 = nodes.Node("n1")
        n2 = nodes.Node("n2")

        out = n1 + n2
        self.assertTrue( functions.add is out._obj )
        self.assertTrue( n1 in out._args )
        self.assertTrue( n2 in out._args )
        self.assertTrue( {} == out._kw )
        self.assertTrue( out in n1._viewers )
        self.assertTrue( out in n2._viewers )

        out = n1 - n2
        self.assertTrue( functions.subtract is out._obj )
        self.assertTrue( n1 in out._args )
        self.assertTrue( n2 in out._args )
        self.assertTrue( {} == out._kw )
        self.assertTrue( out in n1._viewers )
        self.assertTrue( out in n2._viewers )

        out = n1 * n2
        self.assertTrue( functions.multiply is out._obj )
        self.assertTrue( n1 in out._args )
        self.assertTrue( n2 in out._args )
        self.assertTrue( {} == out._kw )
        self.assertTrue( out in n1._viewers )
        self.assertTrue( out in n2._viewers )

        out = n1 / n2
        self.assertTrue( functions.divide is out._obj )
        self.assertTrue( n1 in out._args )
        self.assertTrue( n2 in out._args )
        self.assertTrue( {} == out._kw )
        self.assertTrue( out in n1._viewers )
        self.assertTrue( out in n2._viewers )

        out = n1 ** n2
        self.assertTrue( functions.power is out._obj )
        self.assertTrue( n1 in out._args )
        self.assertTrue( n2 in out._args )
        self.assertTrue( {} == out._kw )
        self.assertTrue( out in n1._viewers )
        self.assertTrue( out in n2._viewers )

        out = n1 ** 2
        self.assertTrue( functions.power is out._obj )
        self.assertTrue( n1 in out._args )
        self.assertTrue( {} == out._kw )
        self.assertTrue( out in n1._viewers )

        out = -n1
        self.assertTrue( functions.negative is out._obj )
        self.assertTrue( n1 in out._args )
        self.assertTrue( {} == out._kw )
        self.assertTrue( out in n1._viewers )

        out = +n1
        self.assertTrue( out is n1 )

        out = abs(n1)
        self.assertTrue( functions.abs is out._obj )
        self.assertTrue( n1 in out._args )
        self.assertTrue( out in n1._viewers )

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

        out = n1
        out -= n2
        self.assertTrue( functions.subtract is out._obj )
        self.assertTrue( n1 in out._args )
        self.assertTrue( n2 in out._args )
        self.assertTrue( {} == out._kw )
        self.assertTrue( out in n1._viewers )
        self.assertTrue( out in n2._viewers )
        self.assertTrue( out.name is n1.name )

        out = n1
        out *= n2
        self.assertTrue( functions.multiply is out._obj )
        self.assertTrue( n1 in out._args )
        self.assertTrue( n2 in out._args )
        self.assertTrue( {} == out._kw )
        self.assertTrue( out in n1._viewers )
        self.assertTrue( out in n2._viewers )
        self.assertTrue( out.name is n1.name )

        out = n1
        out /= n2
        self.assertTrue( functions.divide is out._obj )
        self.assertTrue( n1 in out._args )
        self.assertTrue( n2 in out._args )
        self.assertTrue( {} == out._kw )
        self.assertTrue( out in n1._viewers )
        self.assertTrue( out in n2._viewers )
        self.assertTrue( out.name is n1.name )

        out = n1
        out **= n2
        self.assertTrue( functions.power is out._obj )
        self.assertTrue( n1 in out._args )
        self.assertTrue( n2 in out._args )
        self.assertTrue( {} == out._kw )
        self.assertTrue( out in n1._viewers )
        self.assertTrue( out in n2._viewers )
        self.assertTrue( out.name is n1.name )

        return

    def testPickle(self):
        from pickle import dumps, loads
        n1 = nodes.Node("n1")
        n2 = nodes.Node("n2")
        v1 = TestViewer()
        v2 = TestViewer()
        eq = n1 + n2
        deq = dumps(eq)
        out = loads(deq)
        # Make sure we have the same viewer structure, etc.
        n1p, n2p = out._args
        self.assertEqual(n1p.name, n1.name)
        self.assertEqual(n2p.name, n2.name)
        self.assertTrue(out in n1p._viewers)
        self.assertTrue(out in n2p._viewers)
        return

if __name__ == "__main__":

    unittest.main()

