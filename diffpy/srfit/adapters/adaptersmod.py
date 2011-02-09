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
"""Adapters tools.

This module defines adapters that interface between SrFit Residual objects and
client code. The primary goal of these adapters is to faithfully adapt objects
to nodes that are to be used in optimization. A secondary and related goal is
to make it easy to adapt client code for use in optimization. A third goal is
to ensure that the resulting optimization is flexible and fast as much as
possible without intimate knowledge of the optimizer.

adapt
registry
others...

"""

import types
from itertools import chain, izip
from inspect import isfunction, ismethod, isclass, isbuiltin

import numpy

from diffpy.srfit.util.getcallargs import getcallargs
from diffpy.srfit.util import messages, absName
from diffpy.srfit.util import isAdapter, isViewable, getParameters
from diffpy.srfit.util import isConstrained
from diffpy.srfit.adapters.nodes import Node, Parameter

__all__ = ["adapt"]

# Helper getters and setters for adapters. Some of these are classes rather
# than factories for the sake of picklability.
class attrgetter(object):
    """Gets named attribute value."""
    def __init__(self, name):
        self.name = name
    def __call__(self, obj):
        return getattr(obj, self.name)

class attrsetter(object):
    """Sets named attribute value."""
    def __init__(self, name):
        self.name = name
    def __call__(self, obj, val):
        return setattr(obj, self.name, val)

class itemgetter(object):
    """Gets named item value."""
    def __init__(self, key):
        self.key = key
    def __call__(self, obj):
        return obj.__getitem__(self.key)

class itemsetter(object):
    """Set named item value."""
    def __init__(self, key):
        self.key = key
    def __call__(self, obj, val):
        return obj.__setitem__(self.key, val)

def selfgetter(obj):
    """Identity getter."""
    return obj

def nosetter(obj, val):
    """Raises AttributeError. Use when adapter is read-only."""
    msg = "%r cannot be set"%obj
    raise AttributeError(msg)

class funcgetter(object):
    """Getter that evaluates a function.

    This evaluates a function adapted by an UnboundOperator with adapter
    arguments.

    """

    def __init__(self, args, kw):
        """Create the operator getter for a ContainerAdapter.

        args    --  Adapter arguments for the function.
        kw      --  Adapter keyword arguments for the function.

        """
        self._args = args
        self._kw = kw
        return

    def __call__(self, obj):
        """Get the operator's value.

        obj  --  The UnboundOperator that adapts the function.
        
        """
        op = obj._getop()
        adargs = self._args
        adkw = self._kw
        args = [a.get() for a in adargs]
        kw = dict((key, item.get()) for (key,item) in adkw.items())
        retval = op(*args, **kw)
        return retval

# End class funcgetter

class livegetter(object):
    """Getter for "live" adapters.

    This is used by ContainerAdapters to create "live" adapters. This getter holds
    another getter that determines what is retrieved. The '__call__' method
    accepts a function that, when called without arguments, retrieves an object
    that can be passed to the held getter. By passing an ContainerAdapter's 'get'
    method to '__call__', the most current value of the adapted object is used
    by the held getter. This allows ContainerAdapters to switch out the objects
    they adapt without necessarily invalidating their own adapters.

    """

    def __init__(self, getter):
        """Create a live getter.

        getter  --  getter that will retrieve a value using obj(), where obj is
                    passed via __call__.

        """
        self._getter = getter
        return

    def __call__(self, obj):
        """Get the live value.

        obj     --  Function to call for live value of referent.
                    self._getter(obj()) gets the node value.
        
        """
        liveobj = obj()
        return self._getter(liveobj)

# End class livegetter

class livesetter(object):
    """Setter for "live" adapters.

    See an explanation of live adapters in livegetter.

    """

    def __init__(self, setter):
        """Create a live setter.

        setter  --  setter that will set a value using obj(), where obj is
                    passed via __call__.

        """
        self._setter = setter
        return

    def __call__(self, obj, val):
        """Set the live value.

        obj     --  Function to call for live value of referent.
                    self._setter(obj(), val) sets the node value.
        val     --  The value to set.
        
        """
        liveobj = obj()
        return self._setter(liveobj, val)

# End class livesetter

def adapt(obj, name = None, getter = selfgetter, setter = nosetter, ignore =
        []):
    """Adapt an object.

    If obj is already wrapped, it is returned and other arguments are ignored.
    If obj is the only argument then the adapted object will be read-only. Use
    the parameters module for creating variable parameters.

    obj     --  The object being adapted or an object that can be used with
                getter and setter to retrieve and modify the object being
                wrapped.
    name    --  Name for the adapter. If this is None (default), a variation of
                the id of the wrapped object will be used for the name.
    getter  --  getter(obj) gets the adapted object (default selfgetter).
    setter  --  setter(obj, val) sets the adapted object (default None).
    ignore  --  List of attribute names that will not be adapted in a
                ContainerAdapter (default []). See ContainerAdapter.

    Returns the adapter.
    
    """
    if isAdapter(obj):
        return obj
    # Get the adapted object. We need this to get the right adapter.
    val = getter(obj)
    adapter = getAdapter(val)
    # Note juxtaposition of name and obj
    aobj = adapter(name, obj, getter, setter)
    # Check the name. Go for '_labelself'.
    if aobj.name is None:
        aobj.name = aobj._labelself()
    # Assign ignore
    if isinstance(aobj, ContainerAdapter):
        aobj.ignore = set(ignore).union(aobj.ignore)
    return aobj

def getAdapter(obj):
    """Get an adapter for an object."""
    for mro in type.mro(type(obj)):
        adapter = registry.get(mro)
        if adapter is not None:
            return adapter
    return ContainerAdapter

class ParameterAdapter(Parameter):
    """Class for adapting parameter-like objects."""

    def __init__(self, name, obj, getter, setter):
        """See the adapt function."""
        # Object used with getter to get the adapted object
        self._obj = obj
        # Getter for adapted object
        self._getter = getter
        # Setter for adapted object
        self._setter = setter
        # Hub of object network. Used to dispatch notifications efficiently.
        # This will be set by a container.
        self._hub = self
        # Flag indicating if this is to be treated as a function
        self._isfunction = True
        # The function arguments
        self._args = []
        self._kw = {}
        # If this is contained, then keep a reference to the container.
        self._container = None
        Parameter.__init__(self, name = name)
        return

    def _get(self):
        """Get the parameter's value."""
        return self._getter(self._obj)

    def _set(self, val):
        """Set the parameter's value."""
        self._setter(self._obj, val)
        return

    def _makeFunction(self, args, kw):
        """Make this into a function.

        The adapter must already have the proper getter to be a function. This
        makes sure that the function arguments are properly viewed.

        args    --  List of adapted arguments
        kw      --  Dictionary of adapted keyword arguments

        """
        self._isfunction = True
        self._args = args
        self._kw = kw
        # Make sure that I am viewing the arguments.
        [arg._addViewer(self) for arg in args]
        [arg._addViewer(self) for arg in kw.values()]
        return

    def _notify(self, msg):
        """Notify viewers of a change.

        This calls the _respond method of the hub (if self is not its own hub)
        and then asks all viewers to respond.

        msg --  The message to send to viewers. Standard messages are defined
                in the messages module. The response of the viewer is defined
                by its _respond method.

        """
        # Since adapter nodes may be in a network centered around a hub node ,
        # we ask the hub to respond to the message, and it will notify the rest
        # of the network. 
        # XXX We only ask the hub to respond if we are not the hub. Otherwise,
        # the first call to _respond would get past this point and then the
        # flow would continue from this point for a second pass.
        if self is not self._hub: 
            self._hub._respond(msg)
        # Now we can perform normal notification. Since we've already visited
        # the hub, or in the process of doing so, we can temporarily lock the
        # hub. We do this just in case one of our viewers is also in our node
        # network, in which case it will be notified by the hub.
        lval = self._hub._nlocked
        self._hub._nlocked = True
        [viewer._respond(msg) for viewer in tuple(self._viewers)]
        self._hub._nlocked = lval
        return

    def get(self):
        """Get the nodes's value.

        If the parameter is constrained, get the value of the constraint
        instead. This tells the hub to update constraints first.

        """
        if self._value is None:
            if self is self._hub:
                self._updateConstraints()
            self._value = self._get()
        return self._value

    def _updateConstraints(self):
        """Update constraints.

        Update the constrains in the network. This need only be called from
        the hub node in a network. Any change in the
        network invalidates the entire network. When 'get' is called on
        any network node, the chain of livegetters between that node and
        the hub will propagate the call to the hub's 'get' method.
        
        """
        if self._nlocked: return
        self._nlocked = True
        # Update the constraint for every invalidated and constrained node.
        for viewer in filter(isConstrained, self._viewers):
            if viewer._value is None:
                val = viewer._constraint.get()
                viewer._set(val)

        self._nlocked = False
        return

# End class ParameterAdapter

# XXX - Any method or property may change internal information, and therefore
# indirectly modify other attributes. The way we handle this is to treat the
# adaped container as a hub. All adapted containees are viewed by the hub and
# notifications get dispatched from the hub. Thus, when any node changes, they
# are all notified.
# FIXME - change this back to ContainerAdapter
class ContainerAdapter(ParameterAdapter):
    """Adapter for generic python objects.

    The main purpose of the ContainerAdapter is to create attributes on the fly
    that mimic the attributes of the adapted object.  Immutable attributes are
    adapted as ParameterAdapter objects, sub-objects are adapted as
    ContainerAdapters. Adapted attributes refer to the adapted object for their
    value. Methods calls are also adapted and return adapters that refer back
    to the function for their value.
    
    Class Attributes
    ignore      --  A set of method names. The named methods are not adapted.
                    When a call to one of these methods is made, the call is
                    forwarded to the adapted object.

    Properties
    adapters    --  Iterator for all adapted objects created by this object.

    """

    ignore = set()

    adapters = property( lambda self: \
            chain(self._cachedattrs.itervalues(),
                self._cacheditems.itervalues(),
                self._cachedfunctions))

    def __init__(self, name, obj, getter, setter):
        """See the adapt function."""
        ParameterAdapter.__init__(self, name, obj, getter, setter)
        # For caching attribute adapters
        self._cachedattrs = {}
        # For caching item adapters
        self._cacheditems = {}
        # For caching unbound method adapters
        self._cachedunbound = {}
        # For caching spawned functions
        self._cachedfunctions = set()
        return

    def _labelitem(self, idx):
        """Provide a label for adapted items.

        This can be overloaded for list-containers to give names to adapted
        items. This takes prescedent over _labelself, but defers to it when
        returning None (default).
        
        """
        return None

    def __call__(self, *args, **kw):
        """Support for callable containers."""
        unbound = self.__getattr__("__call__")
        adapter = unbound(*args, **kw)
        adapter.name = self.name
        return adapter

    def __getattr__(self, name):
        """Get an attribute, but wrap it first."""
        if name in self._cachedattrs:
            return self._cachedattrs[name]
        if name in self._cachedunbound:
            return self._cachedunbound[name]

        # Check ignore list
        if name in self.ignore:
            obj = self._get()
            return getattr(obj, name)

        # Write a getter and setter for the object we're trying to retrieve. We
        # don't know anything about that object, so we default to getattr and
        # setattr and make them live.
        getter = livegetter(attrgetter(name))
        setter = livesetter(attrsetter(name))

        adapter = self._addattr(name, self.get, getter, setter)

        return adapter

    def _addattr(self, name, obj, getter, setter):
        """Add an object as an attribute and adapt it.

        name    --  Name for the adapter and name by which it is accessed.
        obj     --  The obj required by getter and setter.
        getter  --  Getter for the new attribute.
        setter  --  Setter for the new attribute.
        
        """
        adapter = adapt(obj, name, getter, setter)
        adapter._container = self
        adapter._hub = self._hub

        if isinstance(adapter, MethodAdapter):
            self._cachedunbound[name] = adapter
        # Add the adapter to the network.
        elif isViewable(adapter):
            self._cachedattrs[name] = adapter
            self._addNode(adapter)

        return adapter

    def _addNode(self, node):
        """Add a node to the network.

        Allow the node view the hub, which is in charge of dispatching messages
        to the network.
        
        """
        self._hub._addViewer(node)
        return

    ## Indexing methods

    def __getitem__(self, key):
        """Retrieve and adapt a contained item."""

        # Check for slice
        if isinstance(key, slice):
            indices = key.indices(len(self))
            return [self[i] for i in range(*indices)]

        if key in self._cacheditems:
            return self._cacheditems[key]

        # Create a suitable name. If the key is a string, then keep it.
        # Otherwise, go for _labelitem.
        name = key
        if not isinstance(name, str):
            name = self._labelitem(key)

        # Write a getter and setter for the object we're trying to retrieve. We
        # don't know anything about that object, so we default to itemgetter
        # and itemsetter and make them live.
        getter = livegetter(itemgetter(key))
        setter = livesetter(itemsetter(key))

        # Adapt it
        adapter = adapt(self.get, name, getter, setter)
        adapter._container = self
        adapter._hub = self._hub
        self._cacheditems[key] = adapter

        # Watch it
        self._addNode(adapter)

        return adapter

    def __iter__(self):
        """Get a dictionary-like or list-like iterator."""
        obj = self.get()

        if not hasattr(obj, "__iter__"):
            msg = "'%s' object is not iterable"%obj.__class__.__name__
            raise TypeError(msg)

        # Need both dictionary and list-like. We also need to trigger
        # __getitem__ so that objects are created if necessary.
        if hasattr(obj, "iterkeys"):
            return obj.iterkeys()
        elif hasattr(obj, "keys"):
            return iter(obj.keys())
        else:
            return (self[idx] for idx in xrange(len(obj)))
        return

    def __len__(self):
        obj = self.get()
        return len(obj)

    def keys(self):
        obj = self.get()
        return obj.keys()

    def iterkeys(self):
        obj = self.get()
        return obj.iterkeys()

    def values(self):
        return list(self.itervalues())

    def itervalues(self):
        obj = self.get()
        return (self[key] for key in obj.iterkeys())

    def items(self):
        return zip(self.keys(), self.values())

    def iteritems(self):
        return izip(self.iterkeys(), self.itervalues())


    # Needed for double-dispatch. It is up to the visitor to decide if it wants
    # to treat the container like anything other than a container.
    def _identify(self, visitor):
        """Identify self to a visitor."""
        return visitor.onContainer(self)

    # For pretty-printing.
    def _show(self):
        """Get a detailed description of the container."""
        pars = getParameters(self)
        spars = [par._show() for par in pars]
        spars.sort()
        out = absName(self) + "\n"
        out += "\n  ".join(spars)
        return out

    # For pickling.
    def __setstate__(self, state):
        """Set the state after pickling.

        This is provided to avoid calling __getattr__ during unpickling.

        """
        self.__dict__.update(state)
        return


# End class ContainerAdapter

# FIXME - It is possible that an unbound function could be passed as an
# argument to a function (e.g. map). Do we want to support this? If so, then
# UnboundOperator must be a Node.
class UnboundOperator(object):
    """Factory that generates Operator nodes.

    The __call__ method creates a new ContainerAdapter with a funcgetter for
    retrieving its value. __call__ accepts arguments and kewords, which should
    be nodes.  Multiple calls with the same arguments will return the
    configured ContainerAdapter.

    Attributes
    name    --  Name of the operator

    """

    def __init__(self, name, op, getter = None, setter = None):
        """Create an unbound operator.

        name    --  Name of the operator
        op      --  The operation the Operator is to perform.
        getter  --  Ignored, needed for adapter interface.
        setter  --  Ignored, needed for adapter interface.

        getter and setter are ignored so that this provides a direct interface
        to create UnboundOperators directly. The MethodAdapter class has the
        same functionality, but sticks to the normal adapter interface.

        """
        self.name = name
        self._op = op
        # Cache for objects
        self._cache = {}
        # Generators for cache keys, indexed the same as _cache
        self._keygen = {}
        self._container = None
        self._hub = None
        return

    def _labelself(self):
        """Provide a label for self, if one is not provided."""
        return self.__class__.__name__

    def _getop(self):
        """Get the operation."""
        return self._op

    def _setop(self, op):
        """Set the operator."""
        self._op = op
        return

    def __call__(self, *args, **kw):
        """Call the operation symbolically.

        This sets the arguments for the operation.

        args    --  Arguments for the function.
        kw      --  Keyword arguments for the function.

        The args and kw need not be adapted. These are used to determine a key
        that is used to cache the adapter. Calling this multiple times with the
        same inputs will return the same adapter.

        Returns a new adapter that peforms the operation on the values of args
        and kw.

        Raises TypeError if args and kw are not compatible with the
        requirements of the operation. (This error checking is not performed
        for built-in and non-python function.)
        Raises TypeError if the operation is not callable.
        
        """
        from diffpy.srfit.adapters.adaptersmod import adapt
        # Check and see if we've already created this operator
        adargs = [adapt(arg, "_%s_arg_%i" % (self.name, i)) for i, arg in \
                enumerate(args)]
        kwt = ( (key, adapt(val, key)) for key, val in kw.items() )
        adkw = dict(kwt)

        # Create keys for storing the input for later retrieval.
        operation = self._getop()
        key = self._makeKey(operation, args, kw)
        try: return self._cache[key]
        except KeyError: pass

        # Test the function
        self._checkFunction(args, kw)

        # Spawn an adapter storing the value. The getter for the container will
        # call the function.
        # FIXME - adapt gets the value of the object, which may not be
        # available at the time. Furthermore, it may be time consuming.
        getter = funcgetter(adargs, adkw)
        adapter = adapt(self, self.name, getter, nosetter)
        # Let the adapter know it is a function
        adapter._makeFunction(adargs, adkw)

        ## Store the adapter so we don't have to create it more than once.
        self._keygen[key] = (args, kw)
        self._cache[key] = adapter

        # If I have a container, make sure that the new adapter is properly
        # viewed by it.
        if self._container is not None:
            self._container._addNode(adapter)
            adapter._container = self._container
            adapter._hub = self._hub
            self._container._cachedfunctions.add(adapter)

        return adapter

    def _checkFunction(self, args, kw):
        """Check the function call for errors.

        args  --  The list of arguments
        kw    --  The dictionary of keyword arguments

        """
        callargs = None
        operation = self._getop()

        if not callable(operation):
            msg = "'%s' object is not callable"%type(operation)
            raise TypeError(msg)

        if isinstance(operation, numpy.ufunc):
            nin = operation.nin
            if kw: 
                msg = "'%s' does not accept keywords"%self.name
                raise TypeError(msg)
            nargs = len(args)
            if nargs != nin:
                msg = "%s() takes exactly %i arguments (%i given)" % \
                    (self.name, nin, nargs)
                raise TypeError(msg)

        elif isfunction(operation):
            callargs = getcallargs(operation, *args, **kw)

        # We pull this out to exclude Boost.Python.function objects
        elif ismethod(operation) and isfunction(operation.im_func):
            callargs = getcallargs(operation, *args, **kw)

        elif isclass(operation):
            callargs = getcallargs(operation.__init__, *args, **kw)

        #elif not isbuiltin(operation):
        #    callargs = getcallargs(operation.__call__, *args, **kw)

        # FIXME
        # built-ins, etc. fall through. We don't test them since we don't know
        # the nature of their arguments.

        return

    def _makeKey(self, operation, args, kw):
        """Make a key for the passed args and kw.

        This accepts operation as an argument so we don't have to call _getop
        while unpickling.
        
        """
        try: 
            callargs = getcallargs(operation, *args, **kw)
            items = sorted(callargs.items())
        except TypeError:
            items = chain(args, sorted(kw.items()))
        key = hash("".join(repr(item) for item in items))
        return key

    def __getstate__(self):
        """Get state for pickling.

        This stores the value of _getop so we don't have to call _getop when
        unpickling, which causes recursion.

        """
        state = dict(self.__dict__)
        state["_op"] = self._getop()
        return state

    def __setstate__(self, state):
        """Recache the reconstituted objects since their keys may change."""
        self.__dict__.update(state)
        op = self._op
        for key in self._keygen.keys():
            args, kw = self._keygen.pop(key)
            newkey = self._makeKey(op, args, kw)
            adapter = self._cache.pop(key)
            self._cache[newkey] = adapter
        return

# End class UnboundOperator

class MethodAdapter(UnboundOperator):
    """Class for adapting methods.

    This an UnboundOperator that aids in adapting methods. It is used by
    ContainerAdapter.
    
    """

    def __init__(self, name, obj, getter, setter):
        """See the adapt function."""
        self._obj = obj
        # Getter and setter for operation, not for value
        self._getter = getter
        self._setter = setter
        UnboundOperator.__init__(self, name, None)
        return

    def _getop(self):
        """Get the operation.

        This defers to the getter for the operation.
        
        """
        return self._getter(self._obj)

    def _setop(self, val):
        """Get the operation.
        
        This defers to the setter.
        
        """
        return self._setter(self._obj, val)

# End class MethodAdapter

# Registry of adapters indexed by type of object to be adapted.
registry = {}
registry[types.FunctionType] = UnboundOperator
registry[types.LambdaType] = UnboundOperator
registry[numpy.ufunc] = UnboundOperator
registry[numpy.lib.function_base.vectorize] = UnboundOperator
registry[types.MethodType] = MethodAdapter
# Parameters
registry[types.BooleanType] = ParameterAdapter
registry[types.FloatType] = ParameterAdapter
registry[types.IntType] = ParameterAdapter
registry[types.LongType] = ParameterAdapter
registry[numpy.floating] = ParameterAdapter
# Everything else is handled by ContainerAdapter, which is the default.
#registry[types.ComplexType] = ContainerAdapter
#registry[types.NoneType] = ContainerAdapter
#registry[types.StringType] = ContainerAdapter
#registry[types.UnicodeType] = ContainerAdapter
#registry[numpy.ndarray] = ContainerAdapter

__id__ = "$Id$"
