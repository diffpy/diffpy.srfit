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

ObjectAdapter
ObjectAdapters are adapters for generic python objects and generate adapters
for the objects' attributes, items and methods.  ObjectAdapters can act like
leaf nodes, but are aware of their sub-objects, and can adapt functions as
well.  Whether a ParameterAdapter or ObjectAdapter is used to adapt a function
depends on the function output.  The ObjectAdapter class is defined in the
adaptersmod module.

"""


from diffpy.srfit.util import messages, hasNode, absName, formatEq

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
        # Cache for the value
        self._value = None
        self._viewers = set()

        # Flag indicating that we are locked from responding to notifications.
        # This is here to prevent cycles in cyclic viewers.
        self._nlocked = False
        return

    value = property( lambda self: self.get(),
            lambda self, val: self.set(val))

    def _get(self):
        """Evaluate this node."""
        raise AttributeError("Cannot get value")

    def get(self):
        """Get the value of this node.

        This is provided as a public interface for _get. Subclasses should
        overload _get rather than this method.

        """
        return self._get()

    def _set(self, val):
        """Set the value of the node, if possible."""
        raise AttributeError("Cannot set value")

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

    # Viewable methods

    def _addViewer(self, other):
        """Add a viewer.

        Viewers get notified of changes when '_notify' is called.
        
        """
        self._viewers.add(other)
        return

    def _removeViewer(self, other):
        """Remove a viewer."""
        self._viewers.discard(other)
        return

    def _notify(self, msg):
        """Notify viewers of a change.

        This unconditionally calls the _respond method of all viewers.

        msg --  The message to send to viewers. Standard messages are defined
                in the messages module. The response of the viewer is defined
                by its _respond method.

        """
        [viewer._respond(msg) for viewer in tuple(self._viewers)]
        return

    def _respond(self, msg):
        """Respond to a notification.

        The behavior of _respond is dependent on the message.

        VALUE_CHANGED   --  Notify viewers.
        VARY_CHANGED    --  Notify viewers.
        
        """
        if self._nlocked: return
        self._nlocked = True
        self._notify(msg)
        self._nlocked = False
        return

    # Method for visitors

    def _identify(self, visitor):
        """Identify self to a visitor."""
        raise NotImplementedError

    # Arithmetic operators. These allow nodes to be composed using normal
    # python arithmetic syntax. E.g., newnode = node1 + node2.
    def __add__(self, other): 
        from diffpy.srfit.fit.functions import add
        return add(self, other)
    def __sub__(self, other): 
        from diffpy.srfit.fit.functions import subtract
        return subtract(self, other)
    def __mul__(self, other): 
        from diffpy.srfit.fit.functions import multiply
        return multiply(self, other)
    def __div__(self, other): 
        from diffpy.srfit.fit.functions import divide
        return divide(self, other)
    def __pow__(self, other): 
        from diffpy.srfit.fit.functions import power
        return power(self, other)
    def __radd__(self, other): 
        from diffpy.srfit.fit.functions import add
        return add(other, self)
    def __rsub__(self, other): 
        from diffpy.srfit.fit.functions import subtract
        return subtract(other, self)
    def __rmul__(self, other): 
        from diffpy.srfit.fit.functions import multiply
        return multiply(other, self)
    def __rdiv__(self, other): 
        from diffpy.srfit.fit.functions import divide
        return divide(other, self)
    def __rpow__(self, other): 
        from diffpy.srfit.fit.functions import power
        return power(other, self)
    # In-place operations transfer the name, but not identity
    def __iadd__(self, other): 
        from diffpy.srfit.fit.functions import add
        return add(self, other).rename(self.name)
    def __isub__(self, other): 
        from diffpy.srfit.fit.functions import subtract
        return subtract(self, other).rename(self.name)
    def __imul__(self, other): 
        from diffpy.srfit.fit.functions import multiply
        return multiply(self, other).rename(self.name)
    def __idiv__(self, other): 
        from diffpy.srfit.fit.functions import divide
        return divide(self, other).rename(self.name)
    def __ipow__(self, other): 
        from diffpy.srfit.fit.functions import power
        return power(self, other).rename(self.name)
    # Unary operations
    def __neg__(self): 
        from diffpy.srfit.fit.functions import negative
        return negative(self)
    def __pos__(self): return self
    def __abs__(self): 
        from diffpy.srfit.fit.functions import abs
        return abs(self)

    def __str__(self):
        return str(self.name)

    def _labelself(self):
        """Provide a label for self, if one is not provided."""
        obj = self.get()
        return obj.__class__.__name__

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
        value   --  Value of the parameter (default None)
        
        """
        Node.__init__(self, name)
        self._value = value
        self._varied = False
        self._constraint = None
        # The message sent by 'set'. By default this is
        # VALUE_CHANGED
        self._valmsg = messages.VALUE_CHANGED
        # The message sent by 'vary' and 'fix'. By default this is
        # VARY_CHANGED
        self._varmsg = messages.VARY_CHANGED
        return

    def _get(self):
        """Get the parameter's value."""
        return self._value

    def get(self):
        """Get the parameter's value.

        If the parameter is constrained, get the value of the constraint
        instead.

        """
        if self._value is None:
            # This will set _value to the constraint value
            if self.isConstrained():
                val = self._constraint.get()
                self._set(val)
            self._value = self._get()
        return self._value

    def _set(self, val):
        """Set the parameter's value."""
        self._value = val
        return

    def _tryset(self, val):
        """Try to set the value.

        This calls _set only if val is not the current value returned by get.

        val     --  The value to try and set

        Raises AttributeError if the parameter is constrained.

        Returns True if the value was set, False otherwise.

        """
        if self.isConstrained():
            raise AttributeError("parameter is constrained")

        notequiv = (val != self.get())
        if notequiv is False:
            return False
        if notequiv is True or notequiv.any():
            self._set(val)
            # We must call _get in case we are an adapter. The value we pass to
            # _set might get altered by the adaptee.
            self._value = self._get()
            return True
        # if not notequiv.any(): falls through
        return False

    def set(self, val):
        """Set the parameter's value.

        This calls _set and then notifies viewers of any change in the
        parameter's value.

        Raises AttributeError if the parameter is constrained.

        Returns self so that mutator methods can be chained.
        
        """
        if self._tryset(val):
            self._notify(self._valmsg)
        return self

    def constrain(self, eq):
        """Constrain this parameter with an equation.

        eq  --  The constraint equation.

        Raises AttributeError if the parameter is already constrained.
        Raises AttributeError if the constraint refrences the parameter.
        Raises AttributeError if the parameter cannot be set.

        Returns self so that mutator methods can be chained.
        
        """
        if self.isConstrained():
            raise AttributeError("parameter is already constrained")
        # Make sure the constraint doesn't cause self-reference
        if hasNode(self, eq):
            raise AttributeError("constraint causes self-reference")
        # Make sure we can set the value of the Parameter. If not, then we
        # can't constrain to it. Our best estimate is the value of the
        # equation.
        self._set(eq.get())
        # We want to be fixed, with no value and then notify viewers of the
        # change.
        self.fix()
        self._constraint = eq
        self._constraint._addViewer(self)
        self._value = None
        # XXX we assume the value changed, this might not be the case
        self._notify(self._valmsg)
        return self

    def unconstrain(self):
        """Unconstrain this parameter.

        Raises AttributeError if the parameter is not constrained.

        Returns self so that mutator methods can be chained.

        """
        if not self.isConstrained():
            raise AttributeError("parameter is not constrained")

        self._constraint._removeViewer(self)
        self._constraint = None
        # XXX we assume the value changed, this might not be the case
        self._notify(self._valmsg)
        return self

    def vary(self, val = None, dovary = True):
        """Set this as varied during a fit.

        val     --  New value for the parameter. If this is None (default), the
                    parameter value will not change.
        dovary  --  Actually vary the parameter (default True). If this is
                    False, the parameter will be fixed instead.

        Returns self so that mutator methods can be chained.
        
        """
        if dovary and self.isConstrained():
            raise AttributeError("parameter is constrained")

        msg = 0
        oldvary = self._varied
        self._varied = bool(dovary)

        # Record whether the variability changed
        if oldvary is not self._varied:
            msg |= self._varmsg

        # Try to set the value and record whether the value changed
        if val is not None and self._tryset(val):
            msg |= self._valmsg

        # Send notification if we have a message to send
        if msg:
            self._notify(msg)

        return self

    def fix(self, val = None):
        """Set this as fixed during a fit. This is the default state.

        val --  New value for the parameter. If this is None (default), the
                parameter value will not change.

        Returns self so that mutator methods can be chained.
        
        """
        self.vary(val, False)
        return self

    def isVaried(self):
        """Indicate if this Parameter is varied."""
        return self._varied

    def isConstrained(self):
        """Indicate if this Parameter is constrained."""
        return (self._constraint is not None)

    def _respond(self, msg):
        """Respond to a notification.

        The behavior of _respond is dependent on the message.

        VALUE_CHANGED   --  Set _value to None and notify viewers.
        VARY_CHANGED    --  Notify viewers.
        
        """
        if self._nlocked: return
        self._nlocked = True
        # If we get a VALUE_CHANGED message and we are constrained, then we
        # invalidate our value so it can be recomputed later.
        if (msg & messages.VALUE_CHANGED):
            self._value = None
        self._notify(msg)
        self._nlocked = False
        return

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
