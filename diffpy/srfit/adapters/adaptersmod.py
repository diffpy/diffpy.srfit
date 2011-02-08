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
from diffpy.srfit.util import isAdapted, isViewable, getParameters
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
    held by a ContainerAdapter.

    """

    def __init__(self):
        """Create the operator getter for a ContainerAdapter.

        container   --  ContainerAdapter using this getter.

        """
        self._container = None
        return

    def __call__(self, obj):
        """Get the operator's value.

        obj     --  The UnboundOperator that created this getter.
        
        """
        op = obj._getop()
        adargs = self._container._args
        adkw = self._container._kw
        args = [a.get() for a in adargs]
        kw = dict((key, item.get()) for (key,item) in adkw.items())
        retval = op(*args, **kw)
        return retval

# End class funcgetter

class livegetter(object):
    """Getter for "live" adapters.

    This is used by ContainerAdapters to create a live getter. This calls the
    obj function passed to __call__ with no arguments to get a current object
    to pass to another getter to obtain the value for a node.  This is used
    with functions that return a ContainerAdapter so that the attributes and
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

        obj     --  The obj from the container.
        
        """
        liveobj = obj()
        val = self._getter(liveobj)
        return val

# End class livegetter

def nosetter(obj, val):
    msg = "%r cannot be set"%obj
    raise AttributeError(msg)

class abandonedgetter(object):
    """Getter for objects that have been abandoned by a container."""
    def __init__(self, abandoner):
        self.abandoner = abandoner
    def __call__(self, obj):
        msg = '"%s" has been detached from "%s"'%(obj, self.abandoner)
        raise AttributeError(msg)

class abandonedsetter(object):
    """Setter for objects that have been abandoned by a container."""
    def __init__(self, abandoner):
        self.abandoner = abandoner
    def __call__(self, obj, val):
        msg = '"%s" has been detatched from "%s"'%(obj, abandoner)
        raise AttributeError(msg)

def adapt(obj, name = None, getter = selfgetter, setter = nosetter, 
        accessors = [], ignore = []):
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
    accessors   -- List of accessor names (default []). See ContainerAdapter.
    ignore  --  List of attribute names that will not be adapted in a
                ContainerAdapter (default []). See ContainerAdapter.

    Returns the adapter.
    
    """
    if isAdapted(obj):
        return obj
    val = getter(obj)
    # Get the appropriate adapter
    adapter = getAdapter(val)
    # Note juxtaposition of name and obj
    aobj = adapter(name, obj, getter, setter)
    # Check the name. Go for '_labelself'.
    if aobj.name is None and hasattr(aobj, "_labelself"):
        aobj.name = aobj._labelself()
    # Assign ignore and accessors
    if isinstance(aobj, ContainerAdapter):
        aobj.ignore = set(ignore).union(aobj.ignore)
        aobj.accessors = set(accessors).union(aobj.accessors)
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
        self._obj = obj
        self._getter = getter
        self._setter = setter
        Parameter.__init__(self, name = name)
        return

    def _get(self):
        """Get the parameter's value."""
        return self._getter(self._obj)

    def _set(self, val):
        """Set the parameter's value."""
        self._setter(self._obj, val)
        return self

# End class ParameterAdapter

# FIXME - Once we have ContainerAdapter working as a function, we can refine
# the accessor machinery.
class _AccessorWrapper(object):
    """Wrapper for accessor methods.

    This is provided to support picklability.
    
    """
    def __init__(self, container, obj, name):
        self.container = container
        self.obj = obj
        self.name = name
        return

    def __call__(self, *args, **kw):
        """Call the accessor and link the output."""
        accessor = getattr(self.obj.__class__, self.name)
        # Get a key for retrieval
        keygen = (self.name, args, kw)
        key = self.container._accessorkey(*keygen)
        # Try to retrieve the adapter.
        adapter = self.container._registeredaccessors.get(key)
        if adapter is not None:
            return adapter

        # If we didn't find the adapter, then we make it. Create a getter
        # and setter. The getter calls accessor with the arguments. The
        # setter is nosetter.
        getter = _AccessorGetter(accessor, args, kw)
        setter = nosetter
        # Make a name
        sargs = (str(a) for a in args)
        skw = ((str(name), str(val)) for name, val in kw.items())
        fsargs = ", ".join(sargs)
        fsargs += ", ".join("%s = %s" % item for item in skw)
        name = "%s(%s)"%(self.name, fsargs)
        adapter = adapt(self.obj, name, getter, setter)
        adapter._keygen = keygen
        adapter._container = self.container
        self.container._addNode(adapter)
        self.container._registeredaccessors[key] = adapter
        return adapter

# End class _AccessorWrapper

class _AccessorGetter(object):
    """Getter for accessors.

    This is provided to support picklability. The following are equivalent.
    getter = lambda obj: accessor(obj, *args, **kw)
    getter = _AccessorGetter(accessor, *args, **kw)
    
    """
    def __init__(self, accessor, args, kw):
        self.accessor = accessor
        self.args = args
        self.kw = kw
        return

    def __call__(self, obj):
        """Call the accessor with obj and pre-determined arguments."""
        return self.accessor(obj, *self.args, **self.kw)

# End class _AccessorGetter


# XXX - Any method or property may change internal information, and therefore
# indirectly modify other attributes. The way we handle this is to treat
# the adaped container as a hub. All adapted containees view the hub, and the
# hub views all containees. When any of these is modified, all get notified of
# the change.
class ContainerAdapter(Node):
    """Adapter for containers.

    The main purpose of the ContainerAdapter is to create attributes on
    the fly that mimic the attributes of the adapted container.  Numerical
    attributes are adapted as Parameter objects, sub-objects are adapted as
    ContainerAdapters. Adapted attributes are "live"; they refer to the adapted
    container for their value. Methods calls are also adapted and either
    return Operators that operate on the (adapted) method arguments, or they
    are treated as accessors for "hidden" attributes. See the class attributes
    for a further explanation of how methods are handled.
    
    Class Attributes
    accessors   --  A set of method names. The named methods are accessors,
                    which give access to sub-objects of the container.
                    Accessors are adapted to return an adapted output. The
                    value of the adapted output is "live", in that it calls the
                    underlying accessor for its value. The adapters returned
                    form accessors cannot be 'set' unless instructed how to do
                    so by the client by assigning the '_getter' attribute. (See
                    the 'adapt' function. The 'obj' held by the adapter is the
                    adapted container.) See the '_wrapAccessor' method for more
                    information.
    ignore      --  A set of method names. The named methods are not adapted.
                    When a call to one of these methods is made, the call is
                    forwarded to the container.

    Properties
    adapters    --  All adapted objects tied to this adapter.

    """

    ignore = set()
    accessors = set()

    adapters = property( lambda self: \
            chain(self._registeredattrs.values(),
                self._registereditems.values(),
                self._registeredaccessors.values()))

    def __init__(self, name, obj, getter, setter):
        """See the adapt function."""
        Node.__init__(self, name)
        self._obj = obj
        self._getter = getter
        self._setter = setter
        self._registeredattrs = {}
        self._registereditems = {}
        self._registeredaccessors = {}
        self._isfunction = False
        self._args = []
        self._kw = {}
        return

    def _labelself(self):
        """Provide a label for self, if one is not provided."""
        obj = self.get()
        return obj.__class__.__name__

    def _labelitem(self, idx):
        """Provide a label for adapted items.

        This can be overloaded for list-containers to give names to adapted
        items. This takes prescedent over _labelself, but defers to it when
        returning None (default).
        
        """
        return None

    def get(self):
        """Get the value of the adapted object."""
        if self._value is None:
            self._updateConstraints()
            self._value = self._get()
        return self._value

    def _get(self):
        """Get the value of the adapted object."""
        return self._getter(self._obj)

    def _set(self, val):
        """Set the value of the adapted object."""
        self._setter(self._obj, val)
        # This may have changed our attributes. Abandon the old adapters.
        # They will throw errors if they are used.
        self._abandonRegistered(self)
        return self

    def __call__(self, *args, **kw):
        """Support for callable containers."""
        unbound = self.__getattr__("__call__")
        assert(unbound._container is self)
        adapter = unbound(*args, **kw)
        adapter.name = self.name
        assert(adapter._container is self)
        assert(adapter in self._viewers)
        assert(self in adapter._viewers)
        return adapter

    def __getattr__(self, name):
        """Get an attribute, but wrap it first."""
        if name in self._registeredattrs:
            return self._registeredattrs[name]

        # Check accessors
        if name in self.accessors:
            return self._wrapAccessor(name)

        # Check ignore list
        if name in self.ignore:
            obj = self._get()
            return getattr(obj, name)

        # Write a getter and setter for the object we're trying to retrieve. We
        # don't know anything about that object, so we default to getattr and
        # setattr.
        getter = attrgetter(name)
        setter = attrsetter(name)
        # If we are a function, then we want to make sure that adapters that we
        # spawn are live. If we are live, then we need the same.
        if self._isfunction or isinstance(self._getter, livegetter):
            getter = livegetter(getter)
            setter = nosetter
            obj = self.get
        else:
            obj = self._get()

        adapter = self._addattr(name, obj, getter, setter)

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
        if hasattr(adapter, "_postview"):
            adapter._postview = self._addNode
        elif isViewable(adapter):
            self._registeredattrs[name] = adapter
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
        if self._container is not None:
            self._container._updateConstraints()
        for adapter in self.adapters:
            adapter._updateConstraints()
        self._nlocked = False
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

    def _accessorkey(self, name, args, kw):
        """Name an accessor."""
        key = str(name)
        sargs = map(str, args)
        key += "".join(sargs)
        tkw = kw.items()
        tkw.sort()
        for k, v in tkw:
            key += str(k) + str(id(v))
        return key

    def _abandonRegistered(self, abandoner):
        """Abandon the items in the registry.

        Items in the registry are abandoned when the adapter is 'set'. This
        assures that adapters for the old adapted object are not used.

        This modifies the _get and _set methods to throw AttribueErrors. This
        will keep clients from interacting with abandoned objects.

        """
        for adapter in self.adapters:
            # Stop watching each other
            adapter._removeViewer(self)
            self._removeViewer(adapter)
            # Modify '_getter' and '_setter' of the adapter. We know these are
            # present since the adapter was generated from this container.
            adapter._getter = abandonedgetter(abandoner)
            adapter._setter = abandonedsetter(abandoner)
            adapter._container = None
            if hasattr(adapter, "_abandonRegistered"):
                adapter._abandonRegistered(abandoner)
            # Note that we don't perform any notifications. These will be
            # handled normally when new adapters get made.
        self._registeredattrs = {}
        self._registereditems = {}
        self._registeredaccessors = {}
        return

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

        if key in self._registereditems:
            return self._registereditems[key]

        # Write a getter and setter for the object we're trying to retrieve. We
        # don't know anything about that object, so we default to itemgetter
        # and itemsetter.
        getter = itemgetter(key)
        setter = itemsetter(key)


        # If we are a function, then we want to make sure that adapters that we
        # spawn are live. If we are live, then we need the same.
        if self._isfunction or isinstance(self._getter, livegetter):
            getter = livegetter(getter)
            setter = nosetter
            obj = self.get
        else:
            obj = self._get()

        # Create a suitable name. If the key is a string, then keep it.
        # Otherwise, go for the default name.
        name = key
        if not isinstance(name, str):
            name = self._labelitem(key)

        # Adapt it
        adapter = adapt(obj, name, getter, setter)
        adapter._container = self
        self._registereditems[key] = adapter

        # Watch it
        self._addNode(adapter)

        return adapter

    def _respond(self, msg):
        """Respond to a notification.

        The behavior of _respond is dependent on the message.

        VALUE_CHANGED       --  Set _value to None and notify viewers.
        VARY_CHANGED        --  Notify viewers.

        If the constraint update changes the value of the parameter,
        VALUE_CHANGED will be tacked on to the message.
        
        """
        if self._nlocked: return
        # If we get a VALUE_CHANGED message, then we invalidate our value so it
        # can be recomputed later.
        if (msg & messages.VALUE_CHANGED):
            self._value = None
        self._notify(msg)
        return

    def _identify(self, visitor):
        """Identify self to a visitor."""
        return visitor.onContainer(self)

    def __setstate__(self, state):
        """Set the state after pickling.

        This is provided so that accessors can be re-keyed after a pickle. The
        keys of accessors are based on ids, which change upon pickling.

        """
        self.__dict__.update(state)
        # Re-key all the accessors. The key is based on object ids, which have
        # been modified.
        oldkeys = self._registeredaccessors.keys()
        for key in oldkeys:
            adapter = self._registeredaccessors.pop(key)
            # Getter should be an _AccessorGetter
            keygen = adapter._keygen
            newkey = self._accessorkey(*keygen)
            self._registeredaccessors[newkey] = adapter
        return

    def _show(self):
        """Get a detailed description of the container."""
        pars = getParameters(self)
        spars = [par._show() for par in pars]
        spars.sort()
        out = absName(self) + "\n"
        out += "\n  ".join(spars)
        return out

# End class ContainerAdapter

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

        """
        self.name = name
        self._op = op
        self._registeredoperators = {}
        self._container = None
        return

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
        key = tuple(chain(adargs, (a[1] for a in sorted(adkw.items()))))
        if key in self._registeredoperators:
            op = self._registeredoperators[key]
            return self._registeredoperators[key]

        # Check for errors.
        operation = self._getop()
        if not callable(operation):
            raise TypeError("'%s' object is not callable"%type(operation))
        if isinstance(operation, numpy.ufunc):
            nin = operation.nin
            if adkw: 
                msg = "'%s' does not accept keywords"%self.name
                raise TypeError(msg)
            nargs = len(adargs)
            if nargs != nin:
                msg = "%s() takes exactly %i arguments (%i given)" % \
                    (self.name, nin, nargs)
                raise TypeError(msg)
        elif isfunction(operation) or ismethod(operation):
            callargs = getcallargs(operation, *args, **kw)
        elif isclass(operation):
            callargs = getcallargs(operation.__init__, *args, **kw)
        elif not isbuiltin(operation):
            callargs = getcallargs(operation.__call__, *args, **kw)
        # built-ins, etc. fall through. We can't introspect them.

        # Spawn a ContainerAdapter
        getter = funcgetter()
        op = ContainerAdapter(self.name, self, getter, nosetter)
        getter._container = op
        # Make sure that the ContainerAdapter is viewing the arguments.
        [adarg._addViewer(op) for adarg in adargs]
        [adarg._addViewer(op) for adarg in adkw.values()]
        op._args = adargs
        op._kw = adkw
        op._isfunction = True

        ## Store the container so we don't have to create it more than once.
        self._registeredoperators[key] = op

        # If we have a container, make sure that the new container is properly
        # viewed by it.
        if self._container is not None:
            self._container._addNode(op)
            op._container = self._container

        return op

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
        self._postview = None
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

    def xxx__call__(self, *args, **kw):
        """Call the operation symbolically.

        This sets the arguments for the operation. This is overloaded to call a
        _postview function that may be set by a ContainerAdapter. This function
        tells the Operator to view the container since changes in the
        configuration of a container may be transmitted through its methods.

        *args   --  Node arguments.
        **kw    --  Named node arguments.
        
        """
        op = UnboundOperator.__call__(self, *args, **kw)
        if self._postview is not None:
            self._postview(op)
        return op

# End class MethodAdapter


# Registry of adapters indexed by type of object to be adapted.
registry = {}
registry[types.FunctionType] = UnboundOperator
registry[types.LambdaType] = UnboundOperator
registry[numpy.ufunc] = UnboundOperator
registry[numpy.lib.function_base.vectorize] = UnboundOperator
registry[types.MethodType] = MethodAdapter
#registry[types.MethodType] = UnboundOperator
# Parameters
registry[types.BooleanType] = ParameterAdapter
registry[types.ComplexType] = ParameterAdapter
registry[types.FloatType] = ParameterAdapter
registry[types.IntType] = ParameterAdapter
registry[types.LongType] = ParameterAdapter
registry[types.NoneType] = ParameterAdapter
#registry[types.StringType] = ParameterAdapter
#registry[types.UnicodeType] = ParameterAdapter
# FIXME - Parameters may need to be indexible, but at the same time, they
# should be constrainable as well. Some sort of constraint mechanism is needed
# by ContainerAdapter, and then we can replace ParameterAdapter with
# ContainerAdapter.
registry[numpy.ndarray] = ParameterAdapter
registry[numpy.floating] = ParameterAdapter

__id__ = "$Id$"
