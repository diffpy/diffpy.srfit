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
    alive forever.  Remove this wrapper from an optional holder set when
    the method object is to be finalized.

    Attributes
    ----------
    function : MethodType
        the unbound function extracted from the wrapped bound method.
    holder : set, optional
        set that should drop this wrapper when it becomes inactive.
    _wref : weakref
        weak reference to the object the wrapped method is bound to.
    """

    __slots__ = ('function', 'holder', '_wref')

    def __init__(self, f, holder=None):
        """Create a weak reference wrapper to bound method.

        Parameters
        ----------
        f : bound MethodType
            instance-bound method to be wrapped.
        holder : set, optional
            set which should discard this object when the method
            associated object gets deallocated.
        """
        # This does not handle builtin methods, but that can be added
        # if necessary.
        self.function = f.__func__.__get__(None, f.im_class)
        self.holder = holder
        cb = self.__make_callback(self.holder)
        self._wref = weakref.ref(f.__self__, cb)
        return


    def __call__(self, *args, **kwargs):
        """Call the wrapped method if the weak-referenced object is alive.

        If that object does not exist, remove this wrapper from the holder
        set and return None.

        Parameters
        ----------
        *args, **kwargs
            same arguments as for the wrapped bound method.
        """
        mobj = self._wref()
        if mobj is None:
            if self.holder is not None:
                self.holder.discard(self)
            return None
        return self.function(mobj, *args, **kwargs)


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
        state = (self.function, self.holder, mobj)
        return state


    def __setstate__(self, state):
        """Restore the weak reference in this wrapper upon unpickling.
        """
        (self.function, self.holder, mobj) = state
        if mobj is None:
            # use a fake weak reference that mimics deallocated object.
            self._wref = self.__mimic_empty_ref
            return
        # Here the referred object exists.
        cb = self.__make_callback(self.holder)
        self._wref = weakref.ref(mobj, cb)
        return


    @staticmethod
    def __mimic_empty_ref():
        return None


    @staticmethod
    def __make_callback(holder):
        if holder is None:
            return None
        def cb(wref):
            holder.difference_update([
                m for m in holder
                if isinstance(m, WeakBoundMethod) and m._wref == wref])
        return cb

# end of class WeakBoundMethod

# ----------------------------------------------------------------------------

def weak_ref(f, holder=None):
    """Create weak-reference wrapper to a bound method.

    Parameters
    ----------
    f : callable
        object-bound method or a plain function.
    holder : set, optional
        set that should drop the returned wrapper when the `f`
        gets deallocated.

    Returns
    -------
    WeakBoundMethod
        when `f` is a bound method.  If `f` is a plain function,
        return `f` and ignore the `holder` argument.
    """
    # NOTE Weak referencing plain functions is probably not needed,
    # because they are already bound to the defining modules.
    rv = f
    if isinstance(f, (types.MethodType, types.BuiltinMethodType)):
        rv = WeakBoundMethod(f, holder=holder)
    return rv
