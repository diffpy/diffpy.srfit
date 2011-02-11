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

        self.assertTrue(n1 in c1._nodes)
        self.assertTrue(n2 in c1._nodes)
        self.assertTrue(n3 in c2._nodes)
        self.assertFalse(c1.isValid(n1))
        self.assertFalse(c1.isValid(n2))
        self.assertFalse(c2.isValid(n3))

        c2.validate(n3)
        self.assertTrue(c2.isValid(n3))
        c1.addNode(n3)
        self.assertFalse(n3 in c1._nodes)
        self.assertTrue(c1._neighbors[n3._cache] == 1)
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
        self.assertFalse(n3 in c1._nodes)
        self.assertTrue(n3 in c2._nodes)
        self.assertTrue(c1._neighbors[n3._cache] == 1)
        c2.validate(n3)
        self.assertTrue(c2.isValid(n3))
        self.assertTrue(n3 in c2._valid)
        c1.removeNode(n3)
        self.assertTrue(n3 in c2._nodes)
        self.assertFalse(n3._cache in c1._neighbors)
        self.assertFalse(c2.isValid(n3))
        self.assertFalse(n3 in c2._valid)

        # Put two nodes in a neighboring network and make sure that the
        # reference count is correct.
        c2.addNode(n1)
        c2.addNode(n2)
        self.assertFalse(n1 in c2._nodes)
        self.assertFalse(n2 in c2._nodes)
        self.assertTrue(n1 in c1._nodes)
        self.assertTrue(n2 in c1._nodes)
        self.assertTrue(c2._neighbors[c1] == 2)
        c2.removeNode(n1)
        self.assertFalse(n1 in c2._nodes)
        self.assertFalse(n2 in c2._nodes)
        self.assertTrue(n1 in c1._nodes)
        self.assertTrue(n2 in c1._nodes)
        self.assertTrue(c2._neighbors[c1] == 1)
        c2.removeNode(n2)
        self.assertFalse(n1 in c2._nodes)
        self.assertFalse(n2 in c2._nodes)
        self.assertTrue(n1 in c1._nodes)
        self.assertTrue(n2 in c1._nodes)
        self.assertFalse(c1 in c2._neighbors)

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


def makeNode(cache):
    """Make a node with a cache."""
    node = TestNode()
    node._cache = cache
    cache.addNode(node)
    return node


class TestNode(object):
    """Node for testing."""
    pass

# End TestNode class

if __name__ == "__main__":

    unittest.main()

