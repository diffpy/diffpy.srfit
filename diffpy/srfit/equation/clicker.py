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

A Clicker is an object for recording changes of state in other objects. The main
functions are 'click' and '__cmp__', which changes the counter of a Clicker and
compares its counter to that of other Clickers, respectively. The idea is to
click a clicker whenever the state that it monitors changes. That clicker can
then be compared to other clickers that are used to monitor other objects.

Clickers of the same type share a global counter that records the total number
of independent clicks among that type of Clicker. Clicking a Clicker increments
this global counter and sets the local counter of the Clicker equal to it.  In
this way, the most recently clicked Clicker will compare at least as large as
any other clicker of its type.  To use this information, one can keep a
reference clicker and compare it with other clickers.  When the reference
compares less than another clicker then the state information monitored by that
clicker has changed since the last time it was compared with the monitor, and
should therefore be reprocessed.  After the reference is compared to all other
Clickers, it can be clicked so that it will compare greater than all other
clickers until one of them is clicked. 

Clickers can be composed in a network and can play the roll of observer, subject
or both.  When clicked, a Clicker will impose its counter on all of its
observers.  Thus, observer >= subject is always true. Composition is performed
with 'addObserver' or 'addSubject'. Note that no effort is made to check for
loops in the observation structure.

The Clicker is created with a the clickerFactory method defined in this module.
This allows one to have several types of clickers with different global counters
for use in different types of comparisons. Since the global counter of a clicker
is tied to its class, Clickers can only be compared with clickers with the same
class. Thus, one should not inherit from Clicker and expect the subclasses to
interoperate.
"""

def clickerFactory():
    """A factory for creating Clicker classes.

    Returns Clicker class
    """

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
            self.update()
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
            """Click this Clicker and all of its observers.

            Observers will be given the same state as this clicker.
            """
            self.__class__._numclicks += 1
            self.update()
            return

        def update(self):
            """Update the local state and that of the observers.

            This sets the local state to the global state and updates all
            observers.
            """
            self._state = self.__class__._numclicks
            for clicker in self._observers:
                clicker.update()
            return

        def __cmp__(self, other):
            """Compare the counter of two Clickers.
            
            Raises TypeError if the Clickers are from different classes.
            """
            if self.__class__ is not other.__class__:
                raise TypeError("Cannot compare Clickers of different types")
            return self._state - other._state

        def __str__(self):
            """String representation."""
            return "%i:%i"%(self._state, self.__class__._numclicks)

    return Clicker


# version
__id__ = "$Id$"

#
# End of file
