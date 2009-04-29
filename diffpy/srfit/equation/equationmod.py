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

See the class documentation for more information.
"""

from .visitors import Evaluator
from .visitors import Validator
from .visitors import ArgFinder
from .visitors import PartFinder
from .literals import Generator
from .literals import Argument

class Equation(Generator):
    """Class for holding and evaluating a Literal tree.

    Instances have attributes that are the non-const Arguments of the tree
    (accessed by name) and a __call__ method that uses an Evaluator vistor to
    evaluate the tree.  It is assumed, but not checked that Arguments have
    unique names.  If this is not the case, then one should keep their own list
    of Arguments.

    The tree is scanned for errors when it is added. Thus, the tree should be
    complete before putting it inside an Equation.

    Attributes
    evaluator   --  An Evaluator instance unique to this Equation
    root        --  The root Literal of the equation tree
    argdict     --  A dictionary of Arguments from the root, indexed by
                    name. This is used by the __call__ method.
    args        --  A list of Arguments, used to preserve the order of the
                    Arguments, which is used by the __call__ method.
    name        --  A name for this Equation.
    clicker     --  A Clicker instance for recording change in the Generator.
    literal     --  An Argument to store the value of te
    """

    def __init__(self, root=None):
        """Initialize.

        root    --  The root node of the Literal tree (optional)
        """
        Generator.__init__(self)
        self.evaluator = None
        self.root = None
        self.argdict = {}
        self.literal = Argument()

        if root is not None:
            self.setRoot(root)
        return

    def __getattr__(self, name):
        """Gives access to the Arguments as attributes."""
        arg = self.argdict.get(name)
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
        if validator.errors:
            m = "Errors found in Literal tree %s\n"%root
            m += "\n".join(validator.errors)
            raise ValueError(m)

        argfinder = ArgFinder(getconsts=False)
        root.identify(argfinder)
        self.args = list(argfinder.args)
        self.argdict = dict( [(arg.name, arg) for arg in argfinder.args] )

        partfinder = PartFinder()
        root.identify(partfinder)
        self.evaluator = Evaluator(parts = bool(partfinder.parts))

        self.root = root
        self.clicker.addSubject(root.clicker)
        return

    def __call__(self, *args, **kw):
        """Call the equation.
        
        New Argument values are acceped as arguments or keyword assignments (or
        both).  The order of accepted arguments is given by the args
        attribute.  The Equation will remember values set in this way.

        Raises
        ValueError when a passed argument cannot be found
        """
        # Process args
        for idx, val in enumerate(args):
            if idx > len(self.args):
                raise ValueError("Too many arguments")
            arg = self.args[idx]
            arg.setValue(val)

        # Process kw
        for name, val in kw.items():
            arg = self.argdict.get(name)
            if arg is None:
                raise ValueError("No argument named '%s' here"%name)
            arg.setValue(val)

        # Evaluate the function
        self.root.identify(self.evaluator)
        self.evaluator.click()
        return self.evaluator.value

    # for the Generator interface
    def generate(self, clicker):
        """Generate the Literal.

        clicker --  A Clicker instance for decision making. It is not up to the
                    Evaluator or any other visitor to decide when it can call
                    this method.  The clicker can be used by the Generator to
                    make that decision.

        This stores the value of the equation in the literal attribute.
        """
        self.literal.setValue( self() )
        return


# version
__id__ = "$Id$"

#
# End of file
