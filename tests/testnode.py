#!/usr/bin/env python

import unittest

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

    def testPickle(self):
        from pickle import dumps, loads
        eq = nodes.Node("n1")
        cache = eq._cache
        deq = dumps((eq, cache))
        eqp, cachep = loads(deq)
        self.assertEqual(eqp.name, eq.name)
        self.assertTrue(eqp._cache is cachep)
        return

if __name__ == "__main__":

    unittest.main()

