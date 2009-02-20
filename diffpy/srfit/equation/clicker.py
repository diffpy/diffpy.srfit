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
"""Class for recording changes in objects.

A Clicker is an object for recording changes of state objects. The main functions are
'click' and '__cmp__', which changes the state of a Clicker and compares its
state to that of other Clickers, respectively. Clickers can be embedded in other
objects and used to record when the state of those objects changes.

A Clicker plays the roll of both observer and subject in the traditional
observer pattern. When clicked, a Clicker will 'click' all of its _observers.
Clickers of the same type share a global state that records the total number of
clicks among that type of object. Clicking a Clicker increments this global
state and sets the local state of the Clicker equal to this value. In this way,
the most recently clicked Clicker will compare at least as large as any other
clicker of its type.  This allows one to keep a reference Clicker to compare
with a pool of other Clickers.  After comparisons are made with the pool, the
reference is be clicked so that it compares larger (more recent) than the
Clickers in the pool.

Note that no effort is made to check for loops in the observation structure.

The Clicker is created with a the clickerFactory method in this module. This
allows one to have several types of clickers with different global state for use
in different types of comparisons.
"""

def clickerFactory():
    """A factory for creating Clicker classes."""

    class Clicker(object):
        """Clicker class for recording state changes."""

        _numclicks = 0

        def __init__(self):
            """Initialize."""
            self._observers = set()
            self._state = 0
            return

        def addObserver(self, other):
            """Add a Clicker that observes this one."""
            self._observers.add(other)
            return

        def addSubject(self, other):
            """Add a clicker to observe."""
            other.addObserver(self)
            return

        def removeObserver(self, other):
            """Remove an observer.

            This has no effect when the passed Clicker is not an observer.
            """
            self._observers.discard(other)
            return

        def removeSubject(self, other):
            """Remove a subject.

            This has no effect when the passed Clicker is not a subject.
            """
            other.removeObserver(self)
            return

        def hasObserver(self, other):
            """Indicate if the passed Clicker is observing this Clicker."""
            return other in self._observers

        def hasSubject(self, other):
            """Indicate if the passed Clicker is a subject of this one."""
            return self in other._observers

        def click(self):
            """Click this Clicker and all of its _observers.

            Observers will be given the same _state as this clicker.
            """
            self.__class__._numclicks += 1
            self._click()
            return

        def _click(self):
            """Increment the local _state without changing the global _state.

            This is used internally. Do not call this method.
            """
            self._state = self.__class__._numclicks
            for clicker in self._observers:
                clicker._click()

        def __cmp__(self, other):
            """Compare the _state of two clickers."""
            return self._state - other._state

    return Clicker


# version
__id__ = "$Id$"

#
# End of file
