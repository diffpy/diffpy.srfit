#!/usr/bin/env python
########################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2008 Trustees of the Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Chris Farrow
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
########################################################################
"""Literal base class used to construct equation trees.

Literals are base pieces of an equation network. The main functionality of a
Literal is to retrieve or calculate its value. A Literal may have a target
other than itself. In this case, the Literal acts as a proxy to its target.

Literals also have an 'identify' method that also acts on behalf of the target.
Vistors perform auxillary operations on an equation network by talking to the
Literals. 

"""

from diffpy.srfit.util.observable import Observable

class Literal(Observable):
    """Abstract class for equation pieces, such as operators and arguments.

    Literal is an Observable. See the Observable class in
    diffpy.srfit.util.observable for what that entails.

    Attributes
    name        --  A name for this Literal. (default None)
    _value      --  A value that the Literal is storing. This is not guaranteed
                    to bevalid.  Use 'getValue' and 'setValue' to retrieve or
                    set a valid value. (default None)
    _target     --  The target of this Literal. (defaults to the Literal)
    _proxies    --  A set of Literals that are acting as proxies to the
                    Literal.
    args       --  List of Literals needed to compute the value of this
                    Literal (other than this Literal, the target, or proxies).

    """ 

    def __init__(self, value = None, name = None):
        """Initialization."""
        Observable.__init__(self)
        self.name = name
        self._value = value
        self._target = self
        self._proxies = set()
        self.args = []
        return

    def setValue(self, val):
        """Set the value of the Literal.

        val --  The value to assign to the target Literal.

        """
        target = self._getDeepTarget()

        notequiv = (val != target._value)
        if notequiv is False:
            return
        if notequiv is True or notequiv.any():
            # Order is important here. We want val = None to notify
            target.notify()
            target._value = val
        # If not notequiv.any() falls through
        return

    def getValue(self):
        """Compute and return the target's value."""
        target = self._getDeepTarget()

        if target._value is None:
            target._updateValue()

        return target._value

    def setTarget(self, target):
        """Set the target of my value.

        target  --  The target literal to proxy. If target is this object, then
                    this object will act on behalf of itself.

        Raises ValueError if the target indirectly refers to this Literal.
        
        """
        self._validateTarget(target)

        # I'm no longer a proxy for my old target
        self._target._proxies.discard(self)
        
        # Check the args of the new target. We want things to stay in a
        # consistent state, so we must restore any changes if we catch an
        # error.
        oldtarget = self._target
        t = target._getDeepTarget()
        self._target = self
        try:
            for arg in t.args:
                self._validateArg(arg)
        except ValueError:
            self._target = oldtarget
            self._target._proxies.add(self)
            raise

        # Change my target
        self._target = target

        # Tell the new target that I'm it's proxy
        if target is not self:
            target._proxies.add(self)
        else:
            self._value = None
        self.notify()
        return

    def notify(self):
        """Notify at self and all observers."""
        Observable.notify(self)
        for p in self._proxies:
            p.notify()
        return

    def _getDeepTarget(self):
        """Get the deep target of this object.

        This recurses through the chain of targets to find the deepest one.
        
        """
        if self._target is self:
            return self
        return self._target._getDeepTarget()

    def flush(self, other):
        """Invalidate my state and notify observers.

        This operates on the target.
        
        """
        target = self._getDeepTarget()
        # We don't have to go any further if we're already flushed
        if target._value is None:
            return
        target._value = None
        target.notify()
        return

    def identify(self, visitor):
        """Identify self to a visitor.

        This operates on the target.

        Do not overload this. Overload _identify.

        """
        target = self._getDeepTarget()
        return target._identify(visitor)

    def _identify(self, visitor):
        """Do the work of identify.

        This should only be called from the target.

        This method must be overloaded by a derived classes.
        
        """
        m = "'%s' must override '_identify'" % self.__class__.__name__
        raise NotImplementedError(m)

    def _updateValue(self):
        """Update my value if possible.

        This should only be called from the target.

        Raises ValueError if the value cannot be updated.

        """
        raise ValueError("I have no value!")

    def _validateTarget(self, target):
        """Validate adding a target to self.

        This does not check the target's args.

        Raises ValueError if we cannot use this target.
        
        """
        for p in self._iterProxies():
            if target is p:
                raise ValueError("Target '%s' causes self-reference"%target)

        return

    def _validateArg(self, arg):
        """Validate an argument.
        
        Raises ValueError if we cannot use this argument.

        """
        # An argument is valid if it does not contain self references when it
        # is added. If a target changes, then all arguments must be retested.

        # Get my target
        target = self._getDeepTarget()

        # Check for my target in the chain of arguments
        atarget = arg._getDeepTarget()
        if atarget is target:
            raise ValueError("'%s' causes self-reference"%arg)
        for a in atarget._iterArgs():
            if a is target:
                raise ValueError("'%s' causes self-reference"%arg)
        return

    def _iterProxies(self):
        """Iterate over Literals that are directly or indirectly my proxies."""
        for p in self._proxies:
            yield p

        for p in self._proxies:
            for l in p._iterProxies():
                yield l

        return

    def _iterArgs(self):
        """Iterate over target literals on which I depend for my value."""
        for a in self.args:
            yield a._getDeepTarget()

        for a in self.args:
            t = a._getDeepTarget()
            for l in t._iterArgs():
                yield l

        return


# version
__id__ = "$Id$"

#
# End of file
