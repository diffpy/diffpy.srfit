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
"""The Equation class for holding and evaluating an equation.

Equation is a functor that holds a Literal tree that defines an equation. It's
__call__ method evaluates the equation at the most recent value of its
Arguments. The non-constant arguments are accessible as attributes of the
Equation instance.

Example
> # make a Literal tree. Here's a simple one
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

"""

from .visitors import Evaluator
from .visitors import Validator
from .visitors import ArgFinder

class Equation(object):
    """Class for holding and evaluating a Literal tree.

    Instances have attributes that are the non-const Arguments of the tree
    (accessed by name) and a __call__ method that uses an Evaluator to evaluate
    the tree.  It is assumed, but not checked that Arguments have unique names.
    If this is not the case, then you shouldn't try to access the arguments
    through this class.

    The tree is scanned for errors when it is added. Thus, the tree should be
    complete before putting it inside an Equation.
    """

    def __init__(self, root=None):
        """Initialize.

        root    --  The root node of the Literal tree (optional)
        """
        self.evaluator = Evaluator()
        self.root = None
        self.args = {}
        self.arglist = [] # to preserve order
        if root is not None:
            self.setRoot(root)
        return

    def __getattr__(self, name):
        """Gives access to the Arguments as attributes."""
        arg = self.args.get(name)
        if arg is None:
            raise AttributeError("No argument named '%s' here"%name)
        return arg

    def setRoot(self, root):
        """Set the root of the Literal tree.

        Raises:
        ValueError if errors are found in the Literal tree.
        """
        validator = Validator()
        root.identify(validator)
        if( validator.errors ):
            m = "Errors found in Literal tree %s\n"%root
            m += "\n".join(validator.errors)
            raise ValueError(m)

        argfinder = ArgFinder(getconsts=False)
        root.identify(argfinder)
        self.arglist = list(argfinder.args)
        self.args = dict( [(arg.name, arg) for arg in argfinder.args] )
        self.root = root
        return

    def __call__(self, *args, **kw):
        """Call the equation.
        
        New Argument values are acceped as arguments or keyword assignments. The
        order of accepted arguments is given by the arglist attribute.  The
        equation will remember values set in this way.

        Raises
        ValueError when a passed argument cannot be found
        """
        # Process args
        for idx, val in enumerate(args):
            if idx > len(self.arglist):
                raise ValueError("Too many arguments")
            arg = self.arglist[idx]
            arg.setValue(val)

        # Process kw
        for name, val in kw.items():
            arg = self.args.get(name)
            if arg is None:
                raise ValueError("No argument named '%s' here"%name)
            arg.setValue(val)

        # Evaluate the function
        self.root.identify(self.evaluator)
        self.evaluator.clicker.click()
        return self.evaluator.value


# version
__id__ = "$Id$"

#
# End of file
