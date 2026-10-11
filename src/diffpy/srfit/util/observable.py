#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                             michael a.g. aïvázis
#                                  orthologue
#                      (c) 1998-2009  all rights reserved
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Derived from pyre-1.0/packages/pyre/patterns/Observable.py
# See pyre-1.0 for full copyright and license information

__all__ = ["Observable"]


from diffpy.srfit.util.weakrefcallable import weak_ref
from diffpy.utils._deprecator import build_deprecation_message, deprecated

observable_base = "diffpy.srfit.util.observable.Observable"
removal_version = "4.0.0"
addObserver_dep_msg = build_deprecation_message(
    observable_base, "addObserver", "add_observer", removal_version
)
removeObserver_dep_msg = build_deprecation_message(
    observable_base, "removeObserver", "remove_observer", removal_version
)
hasObserver_dep_msg = build_deprecation_message(
    observable_base, "hasObserver", "has_observer", removal_version
)


class Observable(object):
    """Provide notification support for classes that maintain dynamic
    associations with multiple clients.

    Observers, i.e. clients of the observable, register event handlers
    that will be invoked to notify them whenever something interesting
    happens to the observable. The nature of what is being observed is
    defined by Observable descendants and their managers. For example,
    instances of pyre.calc.Node are observable by other nodes whose
    value depends on them so that the dependents can be notified about
    value changes and forced to recompute their own value.

    The event handlers are callables that take a tuple starting with the
    observable instance as their single argument. Register them with
    ``add_observer``, remove them with ``remove_observer``, and invoke
    them all with ``notify``.
    """

    def notify(self, other=()):
        """Notify all observers.

        Parameters
        ----------
        other : tuple, optional
            The extra objects passed to each observer after this
            observable. Default is an empty tuple.
        """
        # build a list before notification, just in case the observer's
        # callback behavior involves removing itself from our callback set
        semaphores = (self,) + other
        for callable in tuple(self._observers):
            callable(semaphores)
        return

    # callback management

    def add_observer(self, callable):
        """Add a callable to the set of observers.

        Bound methods are held by weak reference, and are removed
        automatically once their object is deleted.

        Parameters
        ----------
        callable : callable
            The observer to call when this observable notifies.
        """
        f = weak_ref(callable, fallback=_fb_remove_observer)
        self._observers.add(f)
        return

    def remove_observer(self, callable):
        """Remove a callable from the set of observers.

        Parameters
        ----------
        callable : callable
            The observer to remove.

        Raises
        ------
        KeyError
            If ``callable`` is not an observer of this observable.
        """
        f = weak_ref(callable)
        self._observers.remove(f)
        return

    def has_observer(self, callable):
        """Check whether a callable is in the set of observers.

        Parameters
        ----------
        callable : callable
            The observer to look for.

        Returns
        -------
        bool
            The flag that is True if ``callable`` is an observer.
        """
        f = weak_ref(callable)
        rv = f in self._observers
        return rv

    @deprecated(addObserver_dep_msg)
    def addObserver(self, callable):
        """This function has been deprecated and will be removed in
        version 4.0.0.

        Please use diffpy.srfit.util.observable.Observable.add_observer
        instead.
        """
        return self.add_observer(callable)

    @deprecated(removeObserver_dep_msg)
    def removeObserver(self, callable):
        """This function has been deprecated and will be removed in
        version 4.0.0.

        Please use diffpy.srfit.util.observable.Observable.remove_observer
        instead.
        """
        return self.remove_observer(callable)

    @deprecated(hasObserver_dep_msg)
    def hasObserver(self, callable):
        """This function has been deprecated and will be removed in
        version 4.0.0.

        Please use diffpy.srfit.util.observable.Observable.has_observer
        instead.
        """
        return self.has_observer(callable)

    # meta methods

    def __init__(self, **kwds):
        super(Observable, self).__init__(**kwds)
        self._observers = set()
        return


# end of class Observable

# Local helpers --------------------------------------------------------------


def _fb_remove_observer(fobs, semaphores):
    # Remove WeakBoundMethod `fobs` from the observers of notifying object.
    # This is called from Observable.notify when the WeakBoundMethod
    # associated object dies.
    observable = semaphores[0]
    observable.remove_observer(fobs)
    return


# end of file
