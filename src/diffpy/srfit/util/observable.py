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


class Observable(object):
    """
    Provide notification support for classes that maintain dynamic associations with multiple
    clients.

    Observers, i.e. clients of the observable, register event handlers that will be invoked to
    notify them whenever something interesting happens to the observable. The nature of what is
    being observed is defined by Observable descendants and their managers. For example,
    instances of pyre.calc.Node are observable by other nodes whose value depends on them so
    that the dependents can be notified about value changes and forced to recompute their own
    value.

    The event handlers are callables that take the observable instance as their single
    argument.

    interface:
      addObserver: registers its callable argument with the list of handlers to invoke
      removeObserver: remove an event handler from the list of handlers to invoke
      notify: invoke the registered handlers in the order in which they were registered

    """


    def notify(self, other=()):
        """
        Notify all observers
        """
        # build a list before notification, just in case the observer's callback behavior
        # involves removing itself from our callback set
        semaphors = (self,) + other
        for callable in tuple(self._observers):
            callable(semaphors)
        return


    # callback management
    def addObserver(self, callable):
        """
        Add callable to the set of observers
        """
        f = weak_ref(callable, fallback=_fbRemoveObserver)
        self._observers.add(f)
        return


    def removeObserver(self, callable):
        """
        Remove callable from the set of observers
        """
        f = weak_ref(callable)
        self._observers.remove(f)
        return


    def hasObserver(self, callable):
        """
        True if `callable` is present in the set of observers.
        """
        f = weak_ref(callable)
        rv = f in self._observers
        return rv


    # meta methods
    def __init__(self, **kwds):
        super(Observable, self).__init__(**kwds)
        self._observers = set()
        return

# end of class Observable

# Local helpers --------------------------------------------------------------

def _fbRemoveObserver(fobs, semaphors):
    # Remove WeakBoundMethod `fobs` from the observers of notifying object.
    # This is called from Observable.notify when the WeakBoundMethod
    # associated object dies.
    observable = semaphors[0]
    observable.removeObserver(fobs)
    return

# end of file
