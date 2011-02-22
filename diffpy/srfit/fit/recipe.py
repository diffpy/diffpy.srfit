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
"""Recipe class for creating and organizing fit objects."""

from diffpy.srfit.adapters import adapt
from diffpy.srfit.adapters.nodes import Node
from diffpy.srfit.fit import Parameter, Var
from diffpy.srfit.util import getVaried, isAdapter
from diffpy.srfit.util.cachemanager import CacheManager

__all__ = ["Recipe"]

class Recipe(object):
    """Object for generating equations and managing their parameters.

    The purpose of this class is to generate equations that can be used as part
    of a fit. Various parts of the interface accepts string equations, which
    can simplify the creation of fit complexes. See the methods and examples
    for use.

    Attributes
    Each top-level object passed or extracted during creation is accessible as
    an attribute of the equation class. These objects are adapted, and so have
    the parameter interface.

    """

    def __init__(self, **kw):
        """Initialize the equation with a string.

        kw      --  Named arguments that define the local namespace.

        """
        # Adapt and store keywords
        self._namespace = {}
        # FIXME - rebuild this when pickled.
        self._adapters = {}
        self.register(**kw)
        # FIXME - need objects for accessing constrained parameters. 
        return

    # FIXME - do __setattr__ to set values directly.

    def __getattr__(self, name):
        """Get an attribute from the namespace."""
        try:
            return self._namespace[name]
        except KeyError:
            msg = "'%s' has no attribute '%s'" % (self, name)
            raise AttributeError(msg)

    def register(self, **kw):
        """Register new objects.

        Object adapters will be available as attributes.
        
        """
        for name, val in kw.items():
            adapter = adapt(val, name)
            self._namespace[name] = adapter
            self._adapters[id(val)] = adapter
        return

    def _getAdapter(self, obj):
        """Get an adapter for an object.

        This retrieves a cached adapter or creates a new one and caches it.
        
        """
        key = id(obj)
        adapter = self._adapters.get(key)
        if adapter is None:
            adapter = adapt(obj)
            self._adapters[key] = adapter
        return adapter

    # XXX Should this call register?
    def vars(self, **kw):
        """Create variables and add them to the namespace."""
        for name, val in kw.items():
            v = Var(name, val)
            self._namespace[name] = v
        return

    def var(self, name, val = None):
        """Create and return variables. It gets added to the namespace."""
        v = Var(name, val)
        self._namespace[name] = v
        return v

    def pars(self, **kw):
        """Create parameters and add them to the namespace."""
        for name, val in kw.items():
            v = Parameter(name, val)
            self._namespace[name] = v
        return

    def par(self, name, val = None):
        """Create and return parameter. It gets added to the namespace."""
        p = Parameter(name, val)
        self._namespace[name] = p
        return p

    def constrain(self, target, *args):
        """Get a constraint creator.

        target  --  Target of the constraint.
        args    --  Modfying objects.

        The args and target define a dependency relationship. When the args are
        modified they indicate that the constraint needs to be reevaluated.
        This reevaluation occurs whenever the target is retrieved from the
        recipe (e.g. recipe.target.value evaluates the constraint).

        Returns a callable with which to define the constraint. The constraint
        callable accepts a function of no arguments or a string. If a function
        is passed, it relies on its globals for evaluation. A string constraint
        is evaluated within the namespace defined by the target and the args.

        """
        # Make sure everything is registered
        tadapter = self._getAdapter(target)
        # Check the arguments
        targs = [self._getAdapter(arg) for arg in args]

        nsset = set(self._namespace.values())
        msg = ""
        if tadapter not in nsset:
            msg = "'%s' is not registered" % target
        for targ in targs:
            if targ not in nsset:
                msg = "'%s' is not registered" % target
        if msg:
            raise ValueError(msg)

        cgen = ConstraintBuilder(self, tadapter, targs)
        return cgen

    def build(self, streq):
        """Build an equation from a string.

        The string will be evaluated in the namespace of registered objects,
        parameters and variables.
        
        """
        ns = dict(self._namespace)
        eq = eval(streq, ns)
        eq._makeFunction([], self._namespace)
        return eq

# End class Recipe

class ConstraintBuilder(object):
    """Builder used by Recipes to build constraints."""

    def __init__(self, recipe, target, args):
        self._recipe = recipe
        self._target = target
        self._args = args
        return

    def __call__(self, con):

        if isinstance(con, str):
            # Create a string getter and adapt it
            ns = dict( (arg.name, arg) for arg in self._args )
            ns[self._target.name] = self._target
            getter = StringGetter(con)
            adapter = adapt(ns, "stringpar", getter)
        else:
            getter = ConGetter()
            adapter = adapt(con, "funcpar", getter)

        # Make the adapter a function. This sets up the functional relationship
        # between the adapter and the arguments.
        adapter._makeFunction(self._args, {})

        # Create a proxy parameter for the target. This will hold the
        # constraint and notify the target whenever its value is invalid.
        proxy = Parameter(self._target.name + "_proxy")
        proxy._cache = self._target._cache
        proxy._cache.addNode(proxy)
        proxy.constrain(adapter)

        # FIXME - associate the constraint with the target in the recipe so it
        # can be looked up, etc.

        return

# End class ConstraintBuilder

class ConGetter(object):
    """Getter for function constraints."""

    def __call__(self, func):
        """Get the function value.

        func    --  The function to call. This takes no arguments.

        Returns the value of the function.

        """
        return func()

# End class ConGetter

class StringGetter(object):
    """Getter for string equations."""

    def __init__(self, eqstr):
        self.eqstr = eqstr
        return

    def __call__(self, ns):
        """Get the value.
        
        ns  --  Dictionary holding the namespace of reference objects.

        Returns None, since we don't know what eqstr does.
        
        """
        valspace = dict((name, ref.value) for name, ref in ns.items())
        exec self.eqstr in valspace
        return None

# End class StringGetter

__id__ = "$Id$"
