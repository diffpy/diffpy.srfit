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

class Literal(object):
    """Abstract class for equation pieces, such as operators and arguments.

    Attributes
    _value      --  A value that the Literal is storing. This is not guaranteed
                    to bevalid.  Use 'getValue' and 'setValue' to retrieve or
                    set a valid value. (default None)
    """ 

    def __init__(self, value = None):
        """Initialization."""
        self._value = value
        return

    def evaluate(self, node):
        """Evaluate self based on node's dependencies."""
        if self._value is None:
            raise ValueError("I have no value!")
        return self._value

    def setValue(self, val):
        """Set the value of the Literal.

        val --  The value to assign to the target Literal.

        """
        self._value = val
        return
        notequiv = (val != self._value)
        if notequiv is False:
            return
        if notequiv is True or notequiv.any():
            self._value = val
        # If not notequiv.any() falls through
        return

# version
__id__ = "$Id$"

#
# End of file
