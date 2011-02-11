#!/usr/bin/env python
"""CacheManager class for managing value caches in a network."""

class CacheManager(object):
    """Object that manages cache states in a network of nodes.

    The purpose of this class is to organize nodes (or other objects) in a
    messaging network so that changes to the nodes can be properly
    communicated. The network is defined by a collection of nodes whose values
    may be co-dependent. If one value in the network changes, all other nodes
    should invalidate their cache under the assumption that their actual value
    may have changed. This pattern is used in functional relationships, e.g.
    f(a,b) depends on a and b. It is also used in a collection of adapters that
    represent attributes of a single python object, where any one of the
    adapters may modify other attributes.
    
    A node may exist in multiple networks, but should only have one cache
    manager. The manager should be stored as the '_cache' attribute. This is
    used to detect whether a cache manager needs to associate with a
    neighboring one, since cache invalidation can cross the boundary between
    these networks.

    Typical Use
    Node's value changes:   Node updates local cache, 'invalidateNetwork',
                            'validate(node)'.
    Node constrained:       Add node to constraint's network ('addNode'). This
                            invalidates the network of the node's cache
                            manager.
    Node unconstrained:     Remove node from constraint's network
                            ('removeNode'). This invalidates the network of the
                            node's cache manager.
    
    

    Attributes
    _nodes      --  Set of all nodes in the network.
    _valid      --  Set of nodes with valid cache values.
    _neighbors  --  Neighboring networks (CacheManager instances). These are
                    stored in a dictionary, where the key is the neighbor, and
                    the value is how many links that neighbor has to this
                    network.

    """

    def __init__(self):
        """Initialize the cache object."""
        # XXX Would like to use weakref.WeakSet, but this is python2.7. Maybe
        # backport later.
        self._nodes = set()
        self._valid = set()
        # XXX collections.Counter would be nice here, but again python2.7. :(
        self._neighbors = {}
        # Lock is set when updating to avoid recursion between co-connected
        # neighbors.
        self._locked = False
        return

    def isValid(self, node):
        """Indicate if a node has a valid state."""
        return node in self._valid

    def addNode(self, node):
        """Add a node to the network.

        node    --  Node to add to the network.

        This marks the node as invalid.  If node's '_cache' attribute is not
        self, store that cache manager as a neighbor. 
        
        When adding a new node to a network, the client must set the node's
        '_cache' to the network '_cache' before calling addNode.

        """
        # Store the node (or its cache manager) and invalidate the network
        cache = node._cache
        if cache is self:
            self._nodes.add(node)
            self.invalidateNetwork()
        else:
            self._neighbors.setdefault(cache, 0)
            self._neighbors[cache] += 1
            cache.invalidateNetwork()
        return

    def removeNode(self, node):
        """Remove a node from this network.

        node    --  Node to remove from the network.

        This removes the node from the nodes and valid sets. A node cannot be
        removed from its own cache manager. If the node is in another network,
        remove this network from the neighbors unless some other node is also
        in that network.

        Raises ValueError if self is the node's cache manager.
        Raises ValueError if the node's cache manager is not a neighbor.

        """
        cache = node._cache
        if cache is self:
            raise ValueError("Cannot remove node from own network")
        else:
            m = self._neighbors.get(cache)
            if m is None:
                msg = "Node '%s' not neighboring '%s'"%(node, self)
                raise ValueError(msg)
            elif m == 1:
                self._neighbors.pop(cache)
            else:
                self._neighbors[cache] -= 1
            cache.invalidateNetwork()
        return

    def validate(self, node):
        """Mark node as having valid cache.

        It is up to a node to call this whenever its cache becomes valid. This
        does not notify the network.

        """
        self._valid.add(node)
        return

    def invalidateNetwork(self):
        """Invalidate the cache of every node in the network.

        A node may call this while changing its value (possibly to something
        valid).
        
        """
        # No more valid nodes
        self._valid.clear()

        # Lock self so neighbors cannot cause recursion.
        if self._locked: return
        self._locked = True

        # Notify neighbors
        for neighbor in self._neighbors:
            neighbor.invalidateNetwork()

        # Unlock self
        self._locked = False
        return

# End class CacheManager
