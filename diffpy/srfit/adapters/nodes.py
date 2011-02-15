#!/usr/bin/env python
########################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2011 Trustees of the Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Chris Farrow
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
########################################################################
"""Base nodes in an evaluation network and tools.

Nodes are nodes in an evaluation network. Their purpose is to hold or compute a
value when asked. Nodes are parameters and variables of a fit, as well as
containers of other nodes and operations on those nodes. The value of an
evaluation network can be refined in a generic way using external optimizes.
Through its nodes, an evaluation network has the power of lazy evaluation:
evaluation of a value is only performed if it may have changed since the last
evaluation.  Nodes can be composed using python arithmetic and function syntax
to construct an evaluation network. Nodes also provide an interface that allow
them to be operated on by visitors via double dispatch.

The functionality of nodes is easily adapted to. Objects that can serve as
nodes need only interface with the '_get' and '_set' methods of the proper node
type. Additional functionality, such as lazy evaluation and node arithmetic
comes for free.  Many object types can be adapted automatically. See the
adaptersmod module for more information.

Node Types

Parameter   
A parameter is a leaf node. It's value can be retrieved directly without
reference to other nodes. Parameters can be varied, with the 'vary' method, in
which case they will be tuned during optimization.  In addition, a parameter
can be constrained to another node, in which case it defers to that node for
its value. This is done in such a way that if the parameter is adapting another
object, the object's value is also controlled by the constraint.

ParameterAdapter
A ParameterAdapter adapts a parameter-like object to the Parameter interface.
It can also adapt a function. The ParameterAdapter class is defined in the
adaptersmod module.

ContainerAdapter
ContainerAdapters are adapters for generic python objects and generate adapters
for the objects' attributes, items and methods.  ContainerAdapters can act like
leaf nodes, but are aware of their sub-objects, and can adapt functions as
well.  Whether a ParameterAdapter or ContainerAdapter is used to adapt a function
depends on the function output.  The ContainerAdapter class is defined in the
adaptersmod module.

"""

from diffpy.srfit.util import hasNode, absName
from diffpy.srfit.util.cachemanager import CacheManager

class Node(object):
    """Nodes in an evaluation network.

    The 'get' and 'set' methods provide the public interface for retrieving a
    node's value. The 'value' property does the same. These methods do the work
    of lazy evaluation, and might also defer to other entities, such as
    constraints. See the 'get' and 'set' methods of specific nodes for more
    information. When the node's value must be recalculated, 'get' and 'set'
    defer to '_get' and '_set'. These methods are provided so that nodes can
    easily be subclassed to serve as adapters to other python objects. In this
    scheme, a subclass overloads _get and _set, while "fancy" operations such
    as lazy evaluation and constraint deferral can be handled by the base
    class.

    Nodes can be composed into networks using python arithmetic syntax.
    In-place arithmetic is available, but it does not preserve the identity of
    the augmented node; it only preserves the name of the node.

    Attributes
    name    --  Name of the node

    Properties
    value   --  The value of the node. Defers to 'get' and 'set' methods.

    """

    def __init__(self, name = None):
        """Initialize this node, optionally with a name."""
        self.name = name
        # Network cache manager
        self._cache = CacheManager()
        self._cache.addNode(self)
        # Local cache value
        self._value = None
        return

    value = property( lambda self: self.get(),
            lambda self, val: self.set(val))

    def _get(self):
        """Evaluate this node."""
        raise ValueError("Cannot get value")

    def get(self):
        """Get the value of this node.

        This is provided as a public interface for _get. Subclasses should
        overload _get rather than this method.

        """
        return self._get()

    def _set(self, val):
        """Set the value of the node, if possible."""
        raise ValueError("Cannot set value")

    def set(self, val):
        """Set the value of this node.

        This is provided as a public interface for _set. Subclasses should
        overload _set rather than this method.

        Returns self so that mutator methods can be chained.

        """
        self._set(val)
        return self

    def rename(self, name):
        """Rename the node.

        Returns self so that mutator methods can be chained.
        
        """
        self.name = name
        return self

    # Method needed by CacheManager

    def _onVary(self):
        """Respond to change in variable state of a network node."""
        return

    # Method for visitors

    def _identify(self, visitor):
        """Identify self to a visitor."""
        raise NotImplementedError

    # Arithmetic operators. These allow nodes to be composed using normal
    # python arithmetic syntax. E.g., newnode = node1 + node2.
    def __add__(self, other): 
        from diffpy.srfit.fit.functions import add_
        return add_(self, other)
    def __sub__(self, other): 
        from diffpy.srfit.fit.functions import subtract_
        return subtract_(self, other)
    def __mul__(self, other): 
        from diffpy.srfit.fit.functions import multiply_
        return multiply_(self, other)
    def __div__(self, other): 
        from diffpy.srfit.fit.functions import divide_
        return divide_(self, other)
    def __pow__(self, other): 
        from diffpy.srfit.fit.functions import power_
        return power_(self, other)
    def __radd__(self, other): 
        from diffpy.srfit.fit.functions import add_
        return add_(other, self)
    def __rsub__(self, other): 
        from diffpy.srfit.fit.functions import subtract_
        return subtract_(other, self)
    def __rmul__(self, other): 
        from diffpy.srfit.fit.functions import multiply_
        return multiply_(other, self)
    def __rdiv__(self, other): 
        from diffpy.srfit.fit.functions import divide_
        return divide_(other, self)
    def __rpow__(self, other): 
        from diffpy.srfit.fit.functions import power_
        return power_(other, self)
    # In-place operations transfer the name, but not identity
    def __iadd__(self, other): 
        from diffpy.srfit.fit.functions import add_
        return add_(self, other).rename(self.name)
    def __isub__(self, other): 
        from diffpy.srfit.fit.functions import subtract_
        return subtract_(self, other).rename(self.name)
    def __imul__(self, other): 
        from diffpy.srfit.fit.functions import multiply_
        return multiply_(self, other).rename(self.name)
    def __idiv__(self, other): 
        from diffpy.srfit.fit.functions import divide_
        return divide_(self, other).rename(self.name)
    def __ipow__(self, other): 
        from diffpy.srfit.fit.functions import power_
        return power_(self, other).rename(self.name)
    # Unary operations
    def __neg__(self): 
        from diffpy.srfit.fit.functions import negative_
        return negative_(self)
    def __pos__(self): return self
    def __abs__(self): 
        from diffpy.srfit.fit.functions import abs_
        return abs_(self)

    def __str__(self):
        return str(self.name)

    def _show(self):
        """Get a detailed description of the node."""
        return str(self)

    def show(self):
        """Print a detailed description of the node."""
        print self._show()

# End class Node

class Parameter(Node):
    """Parameter container for a reference value.

    Parameters are references to a value that may be changed during a
    refinement.

    """

    def __init__(self, name, value = None):
        """Initialize the parameter.

        name    --  Name of the parameter
        value   --  Value of the parameter
        
        """
        Node.__init__(self, name)
        # The local cache. Required by CacheManager.
        self._value = value
        # The local constraint. Required by CacheManager.
        self._constraint = None
        return

    def _get(self):
        """Get the parameter's value."""
        return self._value

    def get(self):
        """Get the parameter's value.

        If the parameter is constrained, get the value of the constraint
        instead.

        """
        if not self._cache.isValid(self):
            # Required to assure that constraints are updated.
            self._cache.updateConstraints()
            # Update our cache. We only need to do this if we're not
            # constrained. 'updateConstraints' will do it for us otherwise. We
            # make this additional check just in case '_get' is slow.
            if not self.isConstrained():
                self._value = self._get()
            # Notify cache manager that our cache is valid.
            self._cache.validate(self)
        return self._value

    def _set(self, val):
        """Set the parameter's value."""
        self._value = val
        return

    def _tryset(self, val):
        """Try to set the value.

        This calls _set only if val is not equal to the cache.  The method is
        provided so that other entities, namely the cache manager, can set the
        value of a node without discussing it with the network.

        val     --  The value to set.

        Returns True if the value was set, False otherwise.

        """
        notequiv = (val != self._value)
        if notequiv is False:
            return False
        if notequiv is True or notequiv.any():
            # Set the value.
            self._set(val)
            # Cache the value. We call '_get' in case the value was somehow
            # rejected or modified by underlying code.
            self._value = self._get()
            return True
        # if not notequiv.any(): falls through
        return False

    def set(self, val):
        """Set the parameter's value.

        Raises ValueError if the parameter is constrained.

        Returns self so that mutator methods can be chained.
        
        """
        if self.isConstrained():
            raise ValueError("parameter is constrained")

        # No matter what happens here, the cache is valid after '_tryset'
        if self._tryset(val):
            self._cache.invalidateNetwork()
        self._cache.validate(self)

        return self

    def constrain(self, eq):
        """Constrain a node with an equation.

        eq  --  The constraint equation.

        Raises ValueError if the node is already constrained.
        Raises ValueError if the constraint refrences the node.
        Raises ValueError if the node cannot be set.

        Returns self so that mutator methods can be chained.
        
        """
        self._cache.constrain(self, eq)
        return self

    def unconstrain(self):
        """Unconstrain this node.

        Raises ValueError if the node is not constrained.

        Returns self so that mutator methods can be chained.

        """
        self._cache.unconstrain(self)
        return self

    def vary(self, val = None, dovary = True):
        """Set this as varied during a fit.

        val     --  New value for the parameter. If this is None, the parameter
                    value will not change.
        dovary  --  Vary the parameter if true, fix it otherwise.

        Returns self so that mutator methods can be chained.
        
        """
        self._cache.vary(self, dovary)
        if val is not None:
            self.set(val)
        return self

    def fix(self, val = None):
        """Set this as fixed during a fit. This is the default state.

        val --  New value for the parameter. If this is None, the parameter
                value will not change.

        Returns self so that mutator methods can be chained.
        
        """
        self.vary(val, False)
        return self

    def isVaried(self):
        """Indicate if this Parameter is varied."""
        return self._cache.isVaried(self)

    def isConstrained(self):
        """Indicate if this Parameter is constrained."""
        return self._cache.isConstrained(self)

    def _identify(self, visitor):
        """Identify self to a visitor."""
        return visitor.onParameter(self)

    def _show(self):
        """Get a string representation of the parameter.

        fullname ~ value (variable)
        fullname = value (fixed)
        fullname : constraint (constrained)

        """
        out = absName(self)
        if self.isVaried():
            out += " ~ " + str(self.value)
        elif self.isConstrained():
            out += " : " + self._constraint._show()
        else:
            out += " = " + str(self.value)
        return out

# End class Parameter
