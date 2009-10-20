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
"""The Equation class for wrapping an evaluation network as a functor.

Equation is a functor that holds an evaluation network and gives access to the
Argument nodes of that network.  __call__ method evaluates the equation at the
most recent value of its Arguments. The non-constant arguments are accessible
as attributes of the Equation instance.

Example
> # make an evaluation network. Here's a simple one
> add = AdditionOperator()
> a = Argument(name="a") # Don't forget to name these!
> b = Argument(name="b")
> add.addLiteral(a)
> add.addLiteral(b)
> # make an Equation instance and pass the root
> eq = Equation(root = add)
> eq(a=3, b=4) # returns 7
> eq(a=2) # remembers b=4, returns 6
> eq.a.setValue(-3)
> eq.b.setValue(3)
> eq() # uses last assignment of a and b, returns 0

See the class documentation for more information.

"""

# IDEA - Evaluate the branch-weight at every node when the root is added. Use
# this to break the evaluation into smaller problems that can be run in
# parallel.

from .visitors import Validator
from .visitors import ArgFinder
from .literals import Argument
from .literals import Operator
from .literals import Node

class Equation(Operator):
    """Class for giving a function interface to a literal tree.

    This class inherits from Operator. See the Operator documentation.

    The tree is scanned for errors when it is added. Thus, the tree should be
    complete before putting it inside an Equation.

    Attributes
    root        --  The root Literal of the evaluation network.
    args        --  A list of Arguments, which is used by the __call__ method.
                    Note that arguments are extracted depth-first starting at
                    the root.

    """

    def __init__(self, root=None):
        """Initialize.

        root    --  The root node of the Literal tree (optional)

        """
        Operator.__init__(self)
        self.root = None

        # Needed by Operator interface
        self.symbol = None
        self.nout = 1
        self.args = []
        self.operation = self.__call__

        if root is not None:
            self.setRoot(root)
        return

    def setRoot(self, root):
        """Set the root of the Literal tree.

        Raises:
        ValueError if errors are found in the Literal tree.

        """
        validator = Validator()
        root.identify(validator)
        if validator.errors:
            m = "Errors found in Literal tree %s\n"%root
            m += "\n".join(validator.errors)
            raise ValueError(m)

        argfinder = ArgFinder(getconsts=False)
        self.args = root.identify(argfinder)
        self.nin = len(self.args)
        self.root = root
        return

    def __call__(self, *args):
        """Call the equation.
        
        New Argument values are acceped as arguments. The order of accepted
        arguments is given by the args attribute.  The Equation will remember
        values set in this way.

        """
        # Process args
        for idx, val in enumerate(args):
            if idx > len(self.args):
                raise ValueError("Too many arguments")
            arg = self.args[idx]
            arg.value = val

        return self.root.value

# version
__id__ = "$Id$"

#
# End of file
