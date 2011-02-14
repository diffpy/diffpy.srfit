#!/usr/bin/env python

import unittest

from diffpy.srfit.util.cachemanager import CacheManager

class TestCacheManager(unittest.TestCase):
    """Test CacheManager class."""

    def test__init__(self):
        """Test __init__."""

        c1 = CacheManager()
        self.assertTrue(c1._nodes == set())
        self.assertTrue(c1._valid == set())
        self.assertTrue(c1._neighbors == {})
        self.assertFalse(c1._locked)
        return

    def testAddNode(self):
        """Test adding nodes to a network."""
        c1 = CacheManager()
        c2 = CacheManager()
        n1 = makeNode(c1)
        n2 = makeNode(c1)
        n3 = makeNode(c2)

        self.assertTrue(c1.isMember(n1))
        self.assertTrue(c1.isMember(n2))
        self.assertTrue(c2.isMember(n3))
        self.assertFalse(c1.isValid(n1))
        self.assertFalse(c1.isValid(n2))
        self.assertFalse(c2.isValid(n3))

        c2.validate(n3)
        self.assertTrue(c2.isValid(n3))
        c1.addNode(n3)
        self.assertFalse(c1.isMember(n3))
        self.assertTrue(c1.isNeighbor(n3))
        self.assertFalse(c2.isValid(n3))
        return

    def testRemoveNode(self):
        """Test removing nodes from a network."""
        c1 = CacheManager()
        c2 = CacheManager()
        n1 = makeNode(c1)
        n2 = makeNode(c1)
        n3 = makeNode(c2)
        c1.addNode(n3)

        # Cannot remove node from own manager
        self.assertRaises(ValueError, c1.removeNode, n1)
        # Cannot remove node from non-neighboring manager
        self.assertRaises(ValueError, c2.removeNode, n1)

        # Test valid removal of node from neighboring manager
        self.assertFalse(c1.isMember(n3))
        self.assertTrue(c2.isMember(n3))
        self.assertTrue(c1.isNeighbor(n3))
        c2.validate(n3)
        self.assertTrue(c2.isValid(n3))
        self.assertTrue(c2.isValid(n3))
        c1.removeNode(n3)
        self.assertTrue(c2.isMember(n3))
        self.assertFalse(c1.isNeighbor(n3))
        self.assertFalse(c2.isValid(n3))

        # Put two nodes in a neighboring network and make sure that the
        # reference count is correct.
        c2.addNode(n1)
        c2.addNode(n2)
        self.assertFalse(c2.isMember(n1))
        self.assertFalse(c2.isMember(n2))
        self.assertTrue(c1.isMember(n1))
        self.assertTrue(c1.isMember(n2))
        self.assertTrue(c2.isNeighbor(n1))
        self.assertTrue(c2.isNeighbor(n2))
        c2.removeNode(n1)
        self.assertFalse(c2.isMember(n1))
        self.assertFalse(c2.isMember(n2))
        self.assertTrue(c1.isMember(n1))
        self.assertTrue(c1.isMember(n2))
        self.assertTrue(c2.isNeighbor(n2))
        # This knowingly fails. isNeighbor does not care about the node, but
        # rather its cache manager. If one node from a network is a neighbor,
        # then they all are.
        # self.assertFalse(c2.isNeighbor(n1))
        c2.removeNode(n2)
        self.assertFalse(c2.isMember(n1))
        self.assertFalse(c2.isMember(n2))
        self.assertTrue(c1.isMember(n1))
        self.assertTrue(c1.isMember(n2))
        self.assertFalse(c2.isNeighbor(n2))
        self.assertFalse(c2.isNeighbor(n1))

        return

    def testValidation(self):
        """Test validation methods."""
        c1 = CacheManager()
        c2 = CacheManager()
        n1 = makeNode(c1)
        n2 = makeNode(c1)
        n3 = makeNode(c2)

        self.assertFalse(c1.isValid(n1))
        self.assertFalse(c1.isValid(n2))

        c1.validate(n1)
        c1.validate(n2)
        self.assertTrue(c1.isValid(n1))
        self.assertTrue(c1.isValid(n2))

        # Check within network
        c1.invalidateNetwork()
        self.assertFalse(c1.isValid(n1))
        self.assertFalse(c1.isValid(n2))

        # Check within neighboring network
        c1.addNode(n3)
        c1.validate(n1)
        c1.validate(n2)
        c2.validate(n3)
        self.assertTrue(c1.isValid(n1))
        self.assertTrue(c1.isValid(n2))
        self.assertTrue(c2.isValid(n3))
        c2.invalidateNetwork()
        self.assertTrue(c1.isValid(n1))
        self.assertTrue(c1.isValid(n2))
        self.assertFalse(c2.isValid(n3))
        c2.validate(n3)
        self.assertTrue(c2.isValid(n3))
        c1.invalidateNetwork()
        self.assertFalse(c1.isValid(n1))
        self.assertFalse(c1.isValid(n2))
        self.assertFalse(c2.isValid(n3))

        # Check within next neighbor networks
        c3 = CacheManager()
        n4 = makeNode(c3)
        c2.addNode(n4)
        c1.validate(n1)
        c1.validate(n2)
        c2.validate(n3)
        c3.validate(n4)
        self.assertTrue(c1.isValid(n1))
        self.assertTrue(c1.isValid(n2))
        self.assertTrue(c2.isValid(n3))
        self.assertTrue(c3.isValid(n4))
        c1.invalidateNetwork()
        self.assertFalse(c1.isValid(n1))
        self.assertFalse(c1.isValid(n2))
        self.assertFalse(c2.isValid(n3))
        self.assertFalse(c3.isValid(n4))

        # Check within cyclic networks
        c2.addNode(n1)
        c1.validate(n1)
        c1.validate(n2)
        c2.validate(n3)
        c3.validate(n4)
        self.assertTrue(c1.isValid(n1))
        self.assertTrue(c1.isValid(n2))
        self.assertTrue(c2.isValid(n3))
        self.assertTrue(c3.isValid(n4))
        c1.invalidateNetwork()
        self.assertFalse(c1.isValid(n1))
        self.assertFalse(c1.isValid(n2))
        self.assertFalse(c2.isValid(n3))
        self.assertFalse(c3.isValid(n4))

        return

    def testConstrain(self):
        """Test the constrain methods."""
        c1 = CacheManager()
        c2 = CacheManager()
        n1 = makeNode(c1)
        n2 = makeNode(c1)
        n3 = makeNode(c2)

        # Test various constraint errors
        self.assertRaises(ValueError, c1.constrain, n3, n1)
        self.assertFalse(c1.isConstrained(n2))
        c1.constrain(n2, n3)
        self.assertRaises(ValueError, c1.constrain, n2, n1)
        self.assertTrue(c1.isConstrained(n2))
        self.assertTrue(c2.isNeighbor(n2))
        c1.unconstrain(n2)
        self.assertFalse(c1.isConstrained(n2))
        self.assertRaises(ValueError, c1.unconstrain, n2)
        self.assertRaises(ValueError, c1.unconstrain, n3)
        self.assertFalse(c2.isNeighbor(n2))

        # Test actual constraints
        n1.value = 1
        n2.value = 2
        n3.value = 3
        self.assertEqual(1, n1.value)
        self.assertEqual(2, n2.value)
        self.assertEqual(3, n3.value)
        eq = n1 + n3
        c1.constrain(n2, n1+n3)
        self.assertTrue(c1.isConstrained(n2))
        c1.updateConstraints()
        self.assertEqual(1, n1.value)
        self.assertEqual(4, n2.value)
        self.assertEqual(3, n3.value)
        n1.value = 0
        c1.updateConstraints()
        self.assertEqual(3, n2.value)

        # This is in the same network as n2. It should not be valid
        self.assertFalse(c1.isValid(n1))
        # n2 should be valid, we just got its value
        self.assertTrue(c1.isValid(n2))
        # n3 is in another network. It should be valid as well.
        self.assertTrue(c2.isValid(n3))
        # The equation should be valid, we just got its value
        self.assertTrue(eq._cache.isValid(eq))

        # Try to make a cyclic constraint
        self.assertRaises(ValueError, eq.constrain, n2)

        return

    def testVary(self):
        """Test the vary methods."""
        c1 = CacheManager()
        c2 = CacheManager()
        n1 = makeNode(c1)
        n2 = makeNode(c1)
        n3 = makeNode(c2)

        self.assertFalse(n1._testvaried)
        self.assertFalse(n2._testvaried)
        self.assertFalse(n3._testvaried)

        c1.vary(n1)
        # Check that _onVary got called
        self.assertTrue(n1._testvaried)
        self.assertTrue(n2._testvaried)
        self.assertFalse(n3._testvaried)
        self.assertTrue(c1.isVaried(n1))
        self.assertFalse(c1.isVaried(n2))
        self.assertFalse(c2.isVaried(n3))

        c1.vary(n1, False)
        n1._testvaried = n2._testvaried = n3._testvaried = False
        self.assertFalse(c1.isVaried(n1))
        self.assertFalse(c1.isVaried(n2))
        self.assertFalse(c2.isVaried(n3))

        self.assertRaises(ValueError, c1.vary, n3)
        c1.addNode(n3)
        self.assertRaises(ValueError, c1.vary, n3)
        c2.vary(n3)
        self.assertFalse(n1._testvaried)
        self.assertFalse(n2._testvaried)
        self.assertTrue(n3._testvaried)
        self.assertFalse(c1.isVaried(n1))
        self.assertFalse(c1.isVaried(n2))
        self.assertTrue(c2.isVaried(n3))

        c2.vary(n3, False)
        n1._testvaried = n2._testvaried = n3._testvaried = False
        n1.vary(c1)
        self.assertTrue(n1._testvaried)
        self.assertTrue(n2._testvaried)
        self.assertTrue(n3._testvaried)
        self.assertTrue(c1.isVaried(n1))
        self.assertFalse(c1.isVaried(n2))
        self.assertFalse(c2.isVaried(n3))
        return


def makeNode(cache):
    """Make a node with a cache."""
    node = TestNode("test")
    node._cache = cache
    cache.addNode(node)
    return node

from diffpy.srfit.fit import Parameter
class TestNode(Parameter):
    def __init__(self, *args, **kw):
        Parameter.__init__(self, *args, **kw)
        self._testvaried = False
    def _onVary(self):
        self._testvaried = True


# End TestNode class

if __name__ == "__main__":

    unittest.main()

