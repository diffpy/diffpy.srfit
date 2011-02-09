#!/usr/bin/env python
# FIXME - need better printing

from diffpy.srfit.util.utilmod import isFunction

class Visitor(object):
    """Abstract class for all visitors to a tree of nodes.

    See implemented visitors for examples of use.
    """

    def onParameter(self, obj):
        """Process a Parameter node."""
        raise NotImplementedError

    def onObject(self, obj):
        """Process a Container node."""
        raise NotImplementedError

    def onFunction(self, obj):
        """Process a node that is being used as a function."""
        raise NotImplementedError

# End class Visitor

class FilterGetter(Visitor):
    """Get nodes according to a filter."""

    def __init__(self, filt):
        self.filt = filt
        self.vals = set()
        self.visited = set()
        return

    def onParameter(self, obj):
        if obj in self.visited: return
        self.visited.add(obj)

        if self.filt(obj):
            self.vals.add(obj)
        if obj.isConstrained():
            obj._constraint._identify(self)
        if isFunction(obj):
            self.onFunction(obj)
        if obj._container is not None:
            obj._container._identify(self)
        return

    def onObject(self, obj):
        if obj in self.visited: return
        self.visited.add(obj)

        if self.filt(obj):
            self.vals.add(obj)
        for arg in obj.adapters:
            arg._identify(self)
        if isFunction(obj):
            self.onFunction(obj)
        if obj.isConstrained():
            obj._constraint._identify(self)
        if obj._container is not None:
            obj._container._identify(self)
        return

    def onFunction(self, obj):
        # Don't check this here since we were redirected from another method
        # where this was already done.
        # if obj in self.visited: return
        # self.visited.add(obj)

        if self.filt(obj):
            self.vals.add(obj)
        for arg in obj._args:
            arg._identify(self)
        for arg in obj._kw.values():
            arg._identify(self)
        return

# End class FilterGetter

class ParameterGetter(Visitor):
    """Get parameters from node network, with an optional filter."""

    def __init__(self, filt = None):
        self.filt = filt
        if self.filt is None:
            self.filt = lambda obj: True
        self.vals = set()
        self.visited = set()
        return

    def onParameter(self, obj):
        if obj in self.visited: return
        self.visited.add(obj)

        if self.filt(obj):
            self.vals.add(obj)
        if obj.isConstrained():
            obj._constraint._identify(self)
        if isFunction(obj):
            self.onFunction(obj)
        if obj._container is not None:
            obj._container._identify(self)
        return

    def onObject(self, obj):
        if obj in self.visited: return
        self.visited.add(obj)

        for arg in obj.adapters:
            arg._identify(self)
        if obj.isConstrained():
            obj._constraint._identify(self)
        if isFunction(obj):
            self.onFunction(obj)
        if obj._container is not None:
            obj._container._identify(self)
        return

    def onFunction(self, obj):
        # Don't check this here since we were redirected from another method
        # where this was already done.
        # if obj in self.visited: return
        # self.visited.add(obj)

        for arg in obj._args:
            arg._identify(self)
        for arg in obj._kw.values():
            arg._identify(self)
        return

# End class ParameterGetter

class NodeFinder(Visitor):
    """Find a node within a tree of nodes.

    After operating on a tree, the 'found' attribute will indicate whether the
    node was found in the tree.
    
    """

    def __init__(self, node):
        """Initialize the finder.

        node    --  The node to find.

        """
        self._node = node
        self.found = False
        return

    def onParameter(self, obj):
        if obj is self._node: self.found = True
        if self.found: return
        if obj.isConstrained():
            obj._constraint._identify(self)
        if isFunction(obj):
            self.onFunction(obj)
        return

    def onObject(self, obj):
        if obj is self._node: self.found = True
        if self.found: return
        for arg in obj.adapters:
            arg._identify(self)
        if obj.isConstrained():
            obj._constraint._identify(self)
        if isFunction(obj):
            self.onFunction(obj)
        return

    def onFunction(self, obj):
        # Don't check this here since we were redirected from another method
        # where this was already done.
        # if obj is self._node: self.found = True
        # if self.found: return
        for arg in obj._args:
            arg._identify(self)
        for arg in obj._kw.values():
            arg._identify(self)
        return

# End class NodeFinder

from diffpy.srfit.fit import functions
class Printer(Visitor):
    """Printer for printing a tree of nodes.

    Attributes:
    output  --  The output generated by the printer.

    """

    _infix = { 
            functions.add : "+",
            functions.subtract : "-",
            functions.multiply : "*",
            functions.divide : "/",
            functions.power : "**",
            }
    _unary = {
            functions.negative : "-",
            }

    def __init__(self):
        """Initialize."""
        self.reset()
        return

    def reset(self):
        """Reset the output string."""
        self.output = ""
        return

    def onParameter(self, obj):
        if isFunction(obj):
            self.onFunction(obj)
            return
        if obj.name is None or obj.name.startswith("_"):
            self.output += str(obj.value)
        else:
            self.output += str(obj.name)
        return

    def onObject(self, obj):
        self.onParameter(obj)
        return

    def onFunction(self, obj):
        instr = self.output
        self.output = ""

        # Deal with infix and unary operations
        if obj._unbound in self.__class__._infix:
            self._onInfix(obj)
        elif obj._unbound in self.__class__._unary:
            self._onUnary(obj)
        else:
            instr += str(obj.name)
            numargs = 0
            for arg in obj._args:
                if numargs != 0: self.output += ", "
                arg._identify(self)
                numargs += 1
            for kw, arg in obj._kw.items():
                if numargs != 0: self.output += ", %s = "%kw
                arg._identify(self)

        if not (self.output.startswith("(") and self.output.endswith(")")):
            self.output = "(%s)"%self.output
        self.output = instr + self.output
        return

    def _onInfix(self, op):
        """Process infix operators."""
        symbol = self.__class__._infix[op._unbound]
        op._args[0]._identify(self)
        self.output += " %s "%symbol
        op._args[1]._identify(self)
        return

    def _onUnary(self, op):
        """Process unary operators."""
        symbol = self.__class__._unary[op._unbound]
        self.output += symbol
        op._args[0]._identify(self)
        return
