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

Literals are base pieces of the equation hierarchy. The 'identify' method
identifies the Literal to a visitor by calling the identifying method of the
vistior. 

Literals have two data attributes: 
name    --  The name of the literal as it would appear in an equation.
clicker --  Records changes in the state of the Literal.
The clicker is used by the Evaluator visitor. See that class for a more
detailed description.

Note that even though an equation hierarchy is designed according to the
visitor pattern, the primary purpose of the hierarchy is to flexibly and
efficiently evaluate equations. To make this work without overly complex
machinery, some of the information needed by the Evaluator visitor is stored in
the data objects.

The clicker is described in the diffpy.equation.clicker module.

"""

from .. import Clicker

class Literal(object):
    """Abstract class for equation pieces, such as operators and arguments.

    Attributes
    name    --  A name for this Literal (default None).
    clicker --  A Clicker instance for recording change in the value

    """ 

    # Required attributes - used for type checking
    name = None
    clicker = None

    def __init__(self):
        """Initialization."""
        self.name = None
        self.clicker = Clicker()
        return

    def identify(self, visitor):
        """Identify self to a visitor."""
        m = "'%s' must override 'identify'" % self.__class__.__name__
        raise NotImplementedError(m)


# version
__id__ = "$Id$"

#
# End of file
