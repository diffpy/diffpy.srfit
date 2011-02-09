#!/usr/bin/env python

import types

import numpy

__all__ = ["makeName", "absName", "formatEq", "hasNode", "getParameters",
"getVaried", "isAdapter", "isViewable", "isVariable", "isVaried", "isFixed",
"isConstrained", "isFunction", "makeArray"]

def makeName(obj):
    """Make a unique name for an object."""
    return "id" + hex(id(obj))

def absName(adapter):
    """Get an absolute name for an adapter.

    This assumes that the adapter is somehow embeded in a ObjectAdapter, in
    which case it has back-links to the top container.

    """
    path = []
    parent = adapter
    while parent is not None:
        path.append(str(parent.name))
        parent = getContainer(parent)
    absname = ".".join(path[::-1])
    return absname

def getContainer(adapter):
    """Get adapter's container, if any.

    Returns the adapter's container, or None if it is not in a container.

    """
    if hasattr(adapter, "_container"):
        return adapter._container
    return None

def formatEq(eq):
    """Format an equation for printing."""
    from diffpy.srfit.util.visitors import Printer
    printer = Printer()
    eq._identify(printer)
    return printer.output

def hasNode(node, eq):
    """Determine if an equation has given node.

    node    --  The node to identify
    eq      --  The equation to investigate

    """
    from diffpy.srfit.util.visitors import NodeFinder
    finder = NodeFinder(node)
    eq._identify(finder)
    return finder.found

def getParameters(eq):
    """Get all parameters from an equation.

    eq  --  The root node of an equation.

    Returns a set of all parameter object.

    """
    from diffpy.srfit.util.visitors import ParameterGetter
    getter = ParameterGetter()
    eq._identify(getter)
    return getter.vals

def getVaried(eq):
    """Get all varied objects from equation.

    eq  --  The root node of an equation.

    Returns a set of all varied objects.

    """
    from diffpy.srfit.util.visitors import ParameterGetter
    getter = ParameterGetter(isVaried)
    eq._identify(getter)
    return getter.vals

# Functions for duck-typing
def isAdapter(obj):
    """Determine if obj is an adapter

    An adapter has "name", "get" and "set" attributes.

    """
    retval = hasattr(obj, "name") 
    retval &= hasattr(obj, "get")
    retval &= hasattr(obj, "set")
    return retval

def isViewable(obj):
    """Determine if obj is viewable.

    A viewable object has the "_addViewer", "_removeViewer", "_notify" and
    "_respond" attributes.
    
    """
    retval = hasattr(obj, "_addViewer") 
    retval &= hasattr(obj, "_removeViewer") 
    retval &= hasattr(obj, "_notify") 
    retval &= hasattr(obj, "_respond") 
    return retval

def isVariable(obj):
    """Determine if obj can be varied.

    A variable object is viewable, and has the "vary", "fix" and "varied"
    attributes.

    """
    retval = isViewable(obj)
    retval &= hasattr(obj, "vary")
    retval &= hasattr(obj, "fix")
    retval &= hasattr(obj, "isVaried")
    return retval

def isConstrainable(obj):
    """Determine if obj can be constrained.

    A constrainable object has the "constrain", "unconstrain" and
    "isConstrained" methods.

    """
    retval = hasattr(obj, "constrain")
    retval &= hasattr(obj, "unconstrain")
    retval &= hasattr(obj, "isConstrained")
    return retval

def isVaried(obj):
    """Determine if obj is variable and varied."""
    return isVariable(obj) and obj.isVaried()

def isConstrained(obj):
    """Determine if obj has "varied" method and obj.varied() is True"""
    return isConstrainable(obj) and obj.isConstrained()

def isFixed(obj):
    """Determine if obj is not varied and not constrained."""
    return not (isVaried(obj) or isConstrained(obj))

def isFunction(obj):
    """Determine if obj is a function node."""
    return hasattr(obj, "_isfunction") and obj._isfunction

def makeArray(obj):
    """Make an array out of an object."""
    if hasattr(obj, "__iter__"):
        return numpy.asarray(obj)
    else:
        return numpy.array([obj])

