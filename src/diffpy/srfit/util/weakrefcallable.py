#!/usr/bin/env python
##############################################################################
#
# diffpy.srfit      Complex Modeling Initiative
#                   (c) 2016 Brookhaven Science Associates,
#                   Brookhaven National Laboratory.
#                   All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
##############################################################################
"""Picklable storage of callable objects using weak references."""

import types
import weakref


class WeakBoundMethod(object):
    """Callable wrapper to a bound method stored as a weak reference.

    This stores a bound method without keeping its object alive, and can
    call a fallback function once that object has been deleted.

    Attributes
    ----------
    function : FunctionType
        The unbound function extracted from the wrapped bound method.
    fallback : FunctionType or None
        The plain function called when the object holding the bound method
        has been deleted. It is called with this wrapper as the first
        argument, followed by any positional and keyword arguments passed
        for the bound method, and can be used to deregister this wrapper.
    _wref : weakref.ref
        The weak reference to the object the wrapped method is bound to.
    _class : type
        The type of the object the method is bound to. This is only used
        for pickling.
    """

    __slots__ = ("function", "fallback", "_wref", "_class")

    def __init__(self, f, fallback=None):
        """Create a weak reference wrapper to a bound method.

        Parameters
        ----------
        f : MethodType
            The instance-bound method to wrap.
        fallback : FunctionType, optional
            The plain function called instead of the bound method once the
            method's object has been deleted. Default is None.
        """
        # This does not handle builtin methods, but that can be added
        # if necessary.
        self.function = f.__func__
        self.fallback = fallback
        self._class = type(f.__self__)
        self._wref = weakref.ref(f.__self__)
        return

    def __call__(self, *args, **kwargs):
        """Call the wrapped method if its object is still alive.

        If that object has been deleted and a fallback function is set,
        call the fallback function instead.

        Parameters
        ----------
        *args, **kwargs
            The arguments passed to the wrapped bound method.

        Returns
        -------
        object
            The return value of the bound method, or of the fallback
            function once the object has been deleted.

        Raises
        ------
        ReferenceError
            If the method's object has been deleted and no fallback
            function is set.
        """
        mobj = self._wref()
        if mobj is not None:
            return self.function(mobj, *args, **kwargs)
        if self.fallback is not None:
            return self.fallback(self, *args, **kwargs)
        emsg = "Object bound to {} does not exist.".format(self.function)
        raise ReferenceError(emsg)

    # support use of this class in hashed collections

    def __hash__(self):
        return hash((self.function, self._wref))

    def __eq__(self, other):
        rv = self.function == other.function and (
            self._wref == other._wref or None is self._wref() is other._wref()
        )
        return rv

    def __ne__(self, other):
        return not self.__eq__(other)

    # support pickling of this type

    def __getstate__(self):
        """Return state with a resolved weak reference."""
        mobj = self._wref()
        nm = self.function.__name__
        amsg = "Unable to pickle this unbound function by name."
        assert self.function is getattr(self._class, nm), amsg
        state = (self._class, nm, self.fallback, mobj)
        return state

    def __setstate__(self, state):
        """Restore the weak reference in this wrapper upon
        unpickling."""
        self._class, nm, self.fallback, mobj = state
        self.function = getattr(self._class, nm)
        if mobj is None:
            # use a fake weak reference that mimics deallocated object.
            self._wref = self.__mimic_empty_ref
            return
        # Here the referred object exists.
        self._wref = weakref.ref(mobj)
        return

    @staticmethod
    def __mimic_empty_ref():
        return None


# end of class WeakBoundMethod


def weak_ref(f, fallback=None):
    """Create weak-reference wrapper to a bound method.

    Parameters
    ----------
    f : callable
        The object-bound method or plain function to wrap.
    fallback : FunctionType, optional
        The plain function called when the object holding ``f`` has been
        deleted. It is called with the wrapper as the first argument,
        followed by the positional and keyword arguments passed for the
        bound method, and can be used to deregister the wrapper. Default
        is None.

    Returns
    -------
    WeakBoundMethod or callable
        The weak wrapper when ``f`` is a bound method. When ``f`` is a
        plain function or already a WeakBoundMethod, ``f`` itself, and
        ``fallback`` is ignored.
    """
    # NOTE Weak referencing plain functions is probably not needed,
    # because they are already bound to the defining modules.
    rv = f
    if isinstance(f, (types.MethodType, types.BuiltinMethodType)):
        rv = WeakBoundMethod(f, fallback=fallback)
    return rv
