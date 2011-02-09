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

    def testPickle(self):
        from pickle import dumps, loads
        eq = nodes.Node("n1")
        deq = dumps(eq)
        eqp = loads(deq)
        self.assertEqual(eqp.name, eq.name)
        return

if __name__ == "__main__":

    unittest.main()

