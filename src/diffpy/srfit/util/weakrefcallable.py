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


"""\
Picklable storage of callable objects using weak references.
"""


import weakref
import types


class WeakBoundMethod(object):
    """\
    Callable wrapper to a bound method stored as a weak reference.

    Support storage of bound methods without keeping the associated objects
    alive forever.  Provide facility for a fallback function to be used
    when method-related object is deleted.

    Attributes
    ----------
    function : MethodType
        the unbound function extracted from the wrapped bound method.
    fallback : FunctionType, optional
        plain function to be called when object holding the bound method
        gets deallocated.  The fallback function is called with this
        object as the first argument followed by any positional and
        keyword arguments passed for the bound method.  The fallback
        function can be used to deregister this wrapper.
    _wref : weakref
        weak reference to the object the wrapped method is bound to.
    _class : type
        the type of the object to which the method is bound.
        This is only used for pickling.
    """

    __slots__ = ('function', 'fallback', '_wref', '_class')

    def __init__(self, f, fallback=None):
        """Create a weak reference wrapper to bound method.

        Parameters
        ----------
        f : bound MethodType
            instance-bound method to be wrapped.
        fallback : FunctionType, optional
            plain function to be called instead of the bound method when
            the method associated object gets deallocated.
        """
        # This does not handle builtin methods, but that can be added
        # if necessary.
        self.function = f.__func__.__get__(None, f.im_class)
        self.fallback = fallback
        self._class = type(f.__self__)
        self._wref = weakref.ref(f.__self__)
        return


    def __call__(self, *args, **kwargs):
        """Call the wrapped method if the weak-referenced object is alive.

        If that object does not exist and the fallback function is defined,
        call the fallback function instead.

        Parameters
        ----------
        *args, **kwargs
            same arguments as for the wrapped bound method.

        Raises
        ------
        ReferenceError
            when the method-bound object does not exist and the fallback
            function is not defined.
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
        rv = (self.function == other.function and
              (self._wref == other._wref or
               None is self._wref() is other._wref()))
        return rv


    def __ne__(self, other):
        return not self.__eq__(other)

    # support pickling of this type

    def __getstate__(self):
        """Return state with a resolved weak reference.
        """
        mobj = self._wref()
        nm = self.function.__name__
        assert self.function == getattr(self._class, nm), \
            "Unable to pickle this unbound function by name."
        state = (self._class, nm, self.fallback, mobj)
        return state


    def __setstate__(self, state):
        """Restore the weak reference in this wrapper upon unpickling.
        """
        (self._class, nm, self.fallback, mobj) = state
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

# ----------------------------------------------------------------------------

def weak_ref(f, fallback=None):
    """Create weak-reference wrapper to a bound method.

    Parameters
    ----------
    f : callable
        object-bound method or a plain function.
    fallback : FunctionType, optional
        plain function to be called when object holding ``f`` gets
        deallocated.  The fallback function is called with the
        wrapper object as the first argument followed by positional
        and keyword arguments passed for the bound method.  The
        fallback function can be used to deregister the wrapper.

    Returns
    -------
    WeakBoundMethod
        when `f` is a bound method.  If `f` is a plain function or
        already of WeakBoundMethod type, return `f` and ignore the
        `fallback` argument.
    """
    # NOTE Weak referencing plain functions is probably not needed,
    # because they are already bound to the defining modules.
    rv = f
    if isinstance(f, (types.MethodType, types.BuiltinMethodType)):
        rv = WeakBoundMethod(f, fallback=fallback)
    return rv
