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
"""Generator class. 

Not all equations can be easily expressed with Arguments, Partitions and
Operators. The Generator class is designed to provide a lower-level interface
for equation builders. It is similar to an Operator, except that a Generator
creates another Literal. This Literal should be defined in the initializer, and
can be updated with the 'generate' method of the Generator. Generators are
treated as leaf nodes (like Arguments and Partitions) of a Literal tree. They
can be dependent on other Literals to update their clicker, but those details
are up to classes that inherit from Generator.

"""

from .abcs import GeneratorABC
from .literal import Literal

import numpy

class Generator(Literal, GeneratorABC):
    """Abstract class for Generator objects.
    
    Attributes
    args    --  A list of Literals that this Generator depends on, modified by
                'addLiteral'. These are ignored by the Evaluator visitor, but
                are used by the ArgFinder visitor, for example.
    name    --  A name for this Generator.
    clicker --  A Clicker instance for recording change in the Generator.
    literal --  The literal modified by the generate method (default None).

    """ 

    # Required attributes - used for type checking
    args = None
    literal = None

    def __init__(self, name = ""):
        """Initialization."""
        Literal.__init__(self)
        self.name = name
        self.literal = None
        self.args = []
        return

    def identify(self, visitor):
        """Identify self to a visitor."""
        visitor.onGenerator(self)
        return

    def addLiteral(self, literal):
        """Add a literal to this generator."""
        self.args.append(literal)
        self.clicker.addSubject(literal.clicker)
        return

    def generate(self, clicker):
        """Generate the Literal.

        clicker --  A Clicker instance for decision making. It is not up to the
                    Evaluator or any other visitor to decide when it can call
                    this method.  The clicker can be used by the Generator to
                    make that decision.

        By default this method does nothing.

        """
        return

    def __str__(self):
        if self.name:
            return "Generator(" + self.name + ")"
        return self.__repr__()

# version
__id__ = "$Id$"

#
# End of file
