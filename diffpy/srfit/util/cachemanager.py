#!/usr/bin/env python
"""CacheManager class for managing value caches in a network."""

from diffpy.srfit.util.utilmod import hasNode

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
    these networks. Nodes also need a local cache named '_value', a local
    handle to a constraint, called '_constraint' and an '_onVary' method that
    handles the response to a change in variable state of a network node (see
    'vary' method.)

    Attributes
    _nodes       -- Set of all nodes in the network.
    _valid       -- Set of nodes with valid cache values.
    _neighbors   -- Neighboring networks (CacheManager instances). These are
                    stored in a dictionary, where the key is the neighbor, and
                    the value is how many links that neighbor has to this
                    network.
    _connodes    -- Set of constrained nodes.
    _varied      -- Set of varied nodes.

    """
    
    def __init__(self):
        """Initialize the cache object."""
        # XXX Would like to use weakref.WeakSet, but this is python2.7. Maybe
        # backport later.
        self._nodes = set()
        self._valid = set()
        # XXX collections.Counter would be nice here, but again python2.7. :(
        self._neighbors = {}
        # Constrained parameters
        self._connodes = set()
        # Varied parameters
        self._varied = set()
        # Lock is set when updating to avoid recursion between co-connected
        # neighbors.
        self._locked = False
        return

    def isMember(self, node):
        """Indicate if a node is a member of this network."""
        return node in self._nodes

    def isNeighbor(self, node):
        """Indicate if a node is a member of a neighboring network."""
        return node._cache in self._neighbors

    def isValid(self, node):
        """Indicate if a node has a valid cache state."""
        return node in self._valid

    def addNode(self, node):
        """Add a node to the network.

        node    --  Node to add to the network.

        If the node's '_cache' is self, the node becomes a member. Otherwise,
        it becomes a neighbor.

        When adding a new node to a network, the client must set the node's
        '_cache' to the network '_cache' before calling addNode, otherwise it
        will not become a member.

        """
        # Store the node (or its cache manager). The node's cache is considered
        # invalid, but this does not invalidate the rest of the network.
        cache = node._cache
        if cache is self:
            self._nodes.add(node)
        else:
            self._neighbors.setdefault(cache, 0)
            self._neighbors[cache] += 1
            cache._valid.discard(node)
        return

    def removeNode(self, node):
        """Remove a neighbor node from this network.

        node    --  Node to remove from the network.

        This removes the node from the collection of neighbors. (Specifically,
        the node's CacheManager reference is decremented in the '_neighbors'
        dictionary.) A node cannot be removed from its own cache manager.

        Note that this does not check the node as a neighbor, only its cache
        manager. Thus, neighbor nodes can be added and different ones can be
        removed for the same effect. In other words, it is up to the node to
        make proper use of this method.

        Raises ValueError if self is the node's cache manager.
        Raises ValueError if the node's cache manager is not a neighbor.

        """
        if self.isMember(node):
            raise ValueError("Cannot remove a member node")
        else:
            cache = node._cache
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

        It is up to a node to call this whenever its cache becomes valid. The
        node should be a member of this CacheManager, but for the sake of speed
        this condition is not verified.

        """
        self._valid.add(node)
        return

    # Constraint handling

    def constrain(self, node, eq):
        """Constrain a node to another node.

        node    --  The node to constrain
        eq      --  The constraint equation.

        Raises ValueError if node is not a member node.
        Raises ValueError if node is already constrained.
        Raises ValueError if eq refrences the node.
        Raises TypeError if node cannot be set.

        """
        if not self.isMember(node):
            msg = "Node '%s' is not in '%s'" % (node, self)
            raise ValueError(msg)
        if self.isConstrained(node):
            raise ValueError("Node is already constrained")
        ## Make sure the constraint doesn't cause self-reference
        if hasNode(eq, node):
            raise ValueError("Constraint causes self-reference")

        # We evaluate the constraint equation to make sure that the node can be
        # set. So we don't waste the effort, we set the value of the node with
        # it so it has a valid cache. 
        val = eq.get()
        node.fix(val)
        eq._cache.addNode(node)

        # Store the constraint.
        node._constraint = eq
        self._connodes.add(node)
        return

    def unconstrain(self, node):
        """Unconstrain a node.

        Raises ValueError if node is not a member node.
        Raises ValueError if the node is not constrained.

        """
        if not self.isMember(node):
            msg = "Node '%s' is not in '%s'" % (node, self)
            raise ValueError(msg)
        if not self.isConstrained(node):
            raise ValueError("Node is not constrained")

        self._connodes.remove(node)
        eq = node._constraint
        eq._cache.removeNode(node)
        # The validity of the cache is not changed by this - we do not need to
        # notify the network.
        return

    def updateConstraints(self):
        """Update all constrained nodes.

        A node should call this before getting its value. This will make sure
        that constraints are properly calculated, and committed before getting
        any other network values.
        
        """
        if self._locked: return

        # We can get out of here if all the constraints have valid cache state.
        if self._connodes <= self._valid: return

        # Get the constrained nodes that are invalid.
        nodes = [node for node in self._connodes if node not in self._valid]

        # Lock up. We'll call 'get' below, which may call 'updateConstraints'.
        self._locked = True

        # Pull all the constraint values.  Keep track of whether any of the
        # constraint values actually get applied.
        invalidate = False
        for node in nodes:
            val = node._constraint.get()
            invalidate |= node._tryset(val)

        # We can unlock now.
        self._locked = False

        # Now we validate the constrained nodes and constraint equations. If
        # any of the constraint values were applied, then invalidate all other
        # nodes.  If none of the constraint values were applied, we do not
        # invalidate the network. 

        if invalidate:
            self.invalidateNetwork()

        self._valid.update(nodes)
        for node in nodes:
            eq = node._constraint
            eq._cache.validate(eq)

        return

    def isConstrained(self, node):
        """Indicate if the node is constrained."""
        return node in self._connodes

    # Variable handling

    def vary(self, node, dovary = True):
        """Set the node to be varied, or not.

        dovary  --  Vary the parameter if True, fix it otherwise.

        Raises ValueError if node is not a member node.

        """
        if not self.isMember(node):
            msg = "Node '%s' is not in '%s'" % (node, self)
            raise ValueError(msg)

        if dovary and node not in self._varied:
            self._varied.add(node)
            self._notifyVary()
        elif node in self._varied and not dovary:
            self._varied.remove(node)
            self._notifyVary()

        return

    def isVaried(self, node):
        """Indicate if the node is varied."""
        return node in self._varied

    # Network messaging

    def invalidateNetwork(self):
        """Invalidate the cache of every node in the network."""
        # Lock self so neighbors cannot cause recursion.
        if self._locked: return
        self._locked = True

        # No more valid nodes
        self._valid.clear()

        # Notify neighbors
        for neighbor in self._neighbors:
            neighbor.invalidateNetwork()

        # Unlock self
        self._locked = False
        return

    def _notifyVary(self):
        """Notify the network of a change in vary state."""
        if self._locked: return
        self._locked = True

        for node in self._nodes:
            node._onVary()

        for neighbor in self._neighbors:
            neighbor._notifyVary()

        self._locked = False
        return


# End class CacheManager
