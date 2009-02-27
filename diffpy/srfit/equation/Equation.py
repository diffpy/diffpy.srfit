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
"""The Equation class for holding and evaluating an equation."""

from .visitors import Evaluator
from .visitors import Validator
from .visitors import ArgFinder

class Equation(object):
    """Class for holding and evaluating a Literal tree.

    The class has data attributes that are the Arguments of the tree (accessed
    by name) and a __call__ method that uses an Evaluator to evaluate the tree.
    It is assumed, but not checked that Arguments have unique names. If this is
    not the case, then you shouldn't try to access the arguments through this
    class.
    """

    def __init__(self, root=None):
        """Initialize.

        root    --  The root node of the Literal tree (optional)
        """
        self.evaluator = Evaluator()
        self.root = None
        self.args = {}
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

        This will check the tree for errors and extract the Arguments. Thus, the
        tree should be complete before calling this method. This will raise a
        ValueError if errors are found.
        """
        validator = Validator()
        root.identify(validator)
        if( validator.errors ):
            m = "Errors found in Literal tree %r\n"%root
            m += "\n".join(validator.errors)
            raise ValueError(m)

        argfinder = ArgFinder()
        root.identify(argfinder)
        self.args = dict( [(arg.name, arg) for arg in argfinder.args] )
        self.root = root
        return

    def __call__(self, **kw):
        """Call the equation.
        
        New Argument values are acceped as keyword assignments and a ValueError
        is thrown in a passed argument cannot be found.
        """
        for name, val in kw.items():
            arg = self.args.get(name)
            if arg is None:
                raise ValueError("No argument named '%s' here"%name)
            arg.setValue(val)

        self.root.identify(self.evaluator)
        self.evaluator.clicker.click()
        return self.evaluator.value



# version
__id__ = "$Id$"

#
# End of file
