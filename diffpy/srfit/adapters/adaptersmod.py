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
    return obj

class funcgetter(object):
    """Evaluates a function.

    This evaluates a function held by and UnboundOperator using the arguments
    set by an UnboundOperator.

    """

    def __init__(self, args, kw):
        """Create the operator getter for a ObjectAdapter.

        args    --  Adapter arguments for the function.
        kw      --  Adapter keyword arguments for the function.

        """
        self._args = args
        self._kw = kw
        return

    def __call__(self, obj):
        """Get the operator's value.

        obj     --  The UnboundOperator that created this getter.
        
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

    This is used by ObjectAdapters to create a live getter. This calls the
    obj function passed to __call__ with no arguments to get a current object
    to pass to another getter to obtain the value for a node.  This is used
    with functions that return a ObjectAdapter so that the attributes and
    items of that container refer back to the function.

    """

    def __init__(self, getter):
        """Create a live getter.

        getter  --  getter that will retrieve a value using the
                    return value of obj.get(), where obj is passed via
                    __call__.

        """
        self._getter = getter
        return

    def __call__(self, obj):
        """Get the live value.

        obj     --  Function to call for live value of referent.
                    self._getter(obj()) gets the node value.
        
        """
        liveobj = obj()
        val = self._getter(liveobj)
        return val

# End class livegetter

class livesetter(object):
    """Setter for "live" adapters.

    This is used by ObjectAdapters to create a live setter. This calls the
    obj function passed to __call__ with no arguments to get a current object
    to pass to another setter to set the value for a node.  This is used
    with functions that return a ObjectAdapter so that the attributes and
    items of that container refer back to the function. With this pattern, it
    is assumed that the object returned by the getter is an internal object,
    and thus its attributes can be modified. If this assumption is not true,
    then setting the value will have no actual effect.

    """

    def __init__(self, setter):
        """Create a live setter.

        setter  --  setter that will set a value using the
                    return value of obj.get(), where obj is passed via
                    __call__.

        """
        self._setter = setter
        return

    def __call__(self, obj, val):
        """Set the live value.

        obj     --  Function to call for live value of referent.
                    self._getter(obj()) gets the node value.
        val     --  The value to set.
        
        """
        liveobj = obj()
        val = self._setter(liveobj, val)
        return val

# End class livesetter


def nosetter(obj, val):
    msg = "%r cannot be set"%obj
    raise AttributeError(msg)

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
                ObjectAdapter (default []). See ObjectAdapter.

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
    if isinstance(aobj, ObjectAdapter):
        aobj.ignore = set(ignore).union(aobj.ignore)
    return aobj

def getAdapter(obj):
    """Get an adapter for an object."""
    for mro in type.mro(type(obj)):
        adapter = registry.get(mro)
        if adapter is not None:
            return adapter
    return ObjectAdapter

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
        return self

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

    def _respond(self, msg):
        """Respond to a notification.

        The behavior of _respond is dependent on the message.

        VALUE_CHANGED       --  Set _value to None and notify viewers.
        VARY_CHANGED        --  Notify viewers.

        """
        if self._nlocked: return
        # If we get a VALUE_CHANGED message, then we invalidate our value so it
        # can be recomputed later.
        if (msg & messages.VALUE_CHANGED):
            self._value = None
        self._notify(msg)
        return

    def _updateConstraints(self):
        """Update constraints within the container network."""
        if self._nlocked: return
        self._nlocked = True
        if self._container is not None:
            # If we're in a container, send the message to update the
            # constraints.
            self._container._updateConstraints()
        if self.isConstrained():
            val = self._constraint.get()
            self._set(val)
        self._nlocked = False
        return


# End class ParameterAdapter


# XXX - perhaps implement __hash__ for this and other nodes
# XXX - Any method or property may change internal information, and therefore
# indirectly modify other attributes. The way we handle this is to treat
# the adaped container as a hub. All adapted containees view the hub, and the
# hub views all containees. When any of these is modified, all get notified of
# the change.
class ObjectAdapter(ParameterAdapter):
    """Adapter for generic python objects.

    The main purpose of the ObjectAdapter is to create attributes on the fly
    that mimic the attributes of the adapted object.  Immutable attributes are
    adapted as ParameterAdapter objects, sub-objects are adapted as
    ObjectAdapters. Adapted attributes refer to the adapted object for their
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

        # Make sure we are mutual observers
        if isinstance(adapter, MethodAdapter):
            self._cachedunbound[name] = adapter
        elif isViewable(adapter):
            self._cachedattrs[name] = adapter
            self._addNode(adapter)

        return adapter

    def _addNode(self, node):
        """Add a node to the network."""
        self._addViewer(node)
        node._addViewer(self)
        return

    def _updateConstraints(self):
        """Update constraints within the container network."""
        if self._nlocked: return
        self._nlocked = True
        for adapter in self.adapters:
            adapter._updateConstraints()
        self._nlocked = False
        ParameterAdapter._updateConstraints(self)
        return

    # XXX - Special adapters might have auxillary data that could be lost when
    # a new accessor is created. This method tries to cache accessors, but
    # cannot do so when the arguments to the accessor are not hashable.
    # Furthermore, when attributes are retrieved by accessor, it cannot be
    # guaranteed that the same adapter is returned for the same reason.
    # IDEA - maybe allow some sort indication that a method is an accessor. For
    # example, container.access.getIt() would identify 'getIt' as an
    # accessor and wrap it up.
    def _wrapAccessor(self, name):
        """Retrieve attributes that are hiding behind methods.

        This wraps an accessor method of the specified name. A call to this
        container's accessor returns an adapted object that in turn calls the
        underlying accessor for its value. The adapter cannot be set unless a
        client specifically assigns the setter. The resulting adapter is linked
        to the container as an attribute would be. The adapter is cached (if
        possible) so multiple calls to the container's accessor give the same
        adapter. This is not possible if the arguments of the accessor are not
        hashable. It is currently not possible to sync adapters created from
        an accessor with their attribute equivalents (if any).

        """
        obj = self.get()
        return _AccessorWrapper(self, obj, name)

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
        self._cacheditems[key] = adapter

        # Watch it
        self._addNode(adapter)

        return adapter

    def _identify(self, visitor):
        """Identify self to a visitor."""
        return visitor.onObject(self)

    def __setstate__(self, state):
        """Set the state after pickling.

        This is provided to avoid calling __getattr__ during unpickling.

        """
        self.__dict__.update(state)
        return

    def _show(self):
        """Get a detailed description of the container."""
        pars = getParameters(self)
        spars = [par._show() for par in pars]
        spars.sort()
        out = absName(self) + "\n"
        out += "\n  ".join(spars)
        return out

# End class ObjectAdapter

# FIXME - It is possible that an unbound function could be passed as an
# argument to a function (e.g. map). Do we want to support this? If so, then
# UnboundOperator must be a Node.
class UnboundOperator(object):
    """Factory that generates Operator nodes.

    The __call__ method creates a new ObjectAdapter with a funcgetter for
    retrieving its value. __call__ accepts arguments and kewords, which should
    be nodes.  Multiple calls with the same arguments will return the
    configured ObjectAdapter.

    Attributes
    name    --  Name of the operator

    """

    def __init__(self, name, op, getter = None, setter = None):
        """Create an unbound operator.

        name    --  Name of the operator
        op      --  The operation the Operator is to perform.
        getter  --  Ignored, needed for adapter interface.
        setter  --  Ignored, needed for adapter interface.

        """
        self.name = name
        self._op = op
        # Cache for objects
        self._cache = {}
        # Generators for cache keys, indexed the same as _cache
        self._keygen = {}
        self._container = None
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

        *args   --  Node arguments.
        **kw    --  Named node arguments.

        Returns a new Operator that peforms the operation on the values of args
        and kw.

        Raises TypeError if args and kw are not compatible with the
        requirements of the operation.
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
        # available at the time. Furthermore, it is time consuming.
        getter = funcgetter(adargs, adkw)
        adapter = adapt(self, self.name, getter, nosetter)
        # Make the adapter a function
        adapter._makeFunction(adargs, adkw)

        ## Store the adapter so we don't have to create it more than once.
        self._keygen[key] = (args, kw)
        self._cache[key] = adapter

        # If I have a container, make sure that the new adapter is properly
        # viewed by it.
        if self._container is not None:
            self._container._addNode(adapter)
            adapter._container = self._container
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
        # built-ins, etc. fall through. We don't introspect them.

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
        """Recache the reconstituted objects since their keys will change."""
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
    ObjectAdapter.
    
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
# Everything else is handled by ObjectAdapter, which is the default.
#registry[types.ComplexType] = ObjectAdapter
#registry[types.NoneType] = ObjectAdapter
#registry[types.StringType] = ObjectAdapter
#registry[types.UnicodeType] = ObjectAdapter
#registry[numpy.ndarray] = ObjectAdapter

__id__ = "$Id$"
