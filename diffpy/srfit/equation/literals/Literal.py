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
"""Literal class. 

Literals are base pieces of the equation hierarchy. The have a single method,
'identify' that identifies the Literal to a visitor. Literal has two data
attributes: 'clicker' records changes in the state of the Literal and 'value',
which records the most recent value of the literal. These are used by the
Evaluator visitor.

Note that even though an equation hierarchy is designed according to the visitor
pattern, the primary purpose of the hierarchy is to flexibly and efficiently
evaluate equations. To make this work without complex machinery, some of the
information needed by the Evaluator visitor is stored in the data objects.

The clicker is described in the diffpy.equation.clicker module.
"""

from .. import Clicker

class Literal(object):
    """Abstract class for equation pieces, such as operators and arguments.
    

    Attributes
    name    --  A name for this Literal.
    clicker --  A Clicker instance for recording change in the value
    value   --  The value stored by the Literal.
    """ 

    def __init__(self):
        """Initialization."""
        self.name = None
        self.clicker = Clicker()
        self.value = None
        return

    def identify(self, visitor):
        """Identify self to a visitor."""
        m = "'%s' must override 'identify'" % self.__class__.__name__
        raise NotImplementedError(m)


# version
__id__ = "$Id$"

#
# End of file
