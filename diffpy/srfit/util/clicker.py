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

A Clicker is an object for recording changes of state in other objects. The
main functions are 'click' and '__cmp__', which changes the counter of a
Clicker and compares its counter to that of other Clickers, respectively. The
way to use a clicker is to click a clicker whenever the state that it monitors
changes. That clicker can then be compared to other clickers that are used to
monitor other objects.

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

Clickers can be composed in a network and can play the roll of observer,
subject or both.  When clicked, a Clicker will impose its counter on all of its
observers.  Thus, observer >= subject is always true. Composition is performed
with 'addObserver' or 'addSubject'. Note that no effort is made to check for
loops in the observation structure.

"""

from _clicker import FlatClicker as Clicker

# version
__id__ = "$Id$"

#
# End of file
